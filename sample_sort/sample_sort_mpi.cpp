#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cassert>

// Function to initialize the array based on the specified input type
void initialize_array(std::vector<int>& arr, const std::string& input_type) {
    if (input_type == "random") {
        srand(static_cast<unsigned int>(arr.size())); 
        for (size_t i = 0; i < arr.size(); i++) {
            arr[i] = rand();
        }
    } else if (input_type == "sorted") {
        for (size_t i = 0; i < arr.size(); i++) {
            arr[i] = i;
        }
    } else if (input_type == "reverse_sorted") {
        for (size_t i = 0; i < arr.size(); i++) {
            arr[i] = arr.size() - i - 1;
        }
    } else if (input_type == "perturbed") {
        for (size_t i = 0; i < arr.size(); i++) {
            arr[i] = i;
        }
        size_t num_perturb = arr.size() / 100;  // Perturb 1% of the elements
        srand(static_cast<unsigned int>(arr.size())); 
        for (size_t i = 0; i < num_perturb; i++) {
            size_t idx1 = rand() % arr.size();
            size_t idx2 = rand() % arr.size();
            std::swap(arr[idx1], arr[idx2]);
        }
    } else {
        std::cerr << "Unknown input_array type: " << input_type << std::endl;
        exit(1);
    }
}

int get_bucket_index(int value, const std::vector<int>& splitters) {
    int bucket = 0;
    while (bucket < splitters.size() && value > splitters[bucket]) {
        bucket++;
    }
    return bucket;
}

int main(int argc, char *argv[]) {
    CALI_CXX_MARK_FUNCTION;

    // MPI variables
    int num_procs, rank, root = 0;
    int i, j, k;
    std::vector<int> input_array, local_array, local_samples, global_samples, buckets, bucket_buffer, local_bucket, sorted_array, output_buffer;

    CALI_MARK_BEGIN("initialization");
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 3) {
        if (rank == root) {
            std::cerr << "Usage: " << argv[0] << " <array_size> <input_type>\n";
        }
        MPI_Finalize();
        exit(1);
    }

    if (rank == root) {
        std::cout << "Array size: " << argv[1] << "\n";
        std::cout << "Number of processes: " << num_procs << "\n";
        std::cout << "Input type: " << argv[2] << "\n";
    }

    int array_size = atoi(argv[1]);
    std::string input_type = argv[2];

    // Validate array size and process count
    if (array_size < num_procs) {
        if (rank == root) {
            std::cerr << "Error: Array size must be >= number of processes\n";
        }
        MPI_Finalize();
        exit(1);
    }
    CALI_MARK_END("initialization");

    // Initialize input array
    if (rank == root) {
        CALI_MARK_BEGIN("data_init_runtime");
        input_array.resize(array_size);
        sorted_array.resize(array_size);
        initialize_array(input_array, input_type);
        CALI_MARK_END("data_init_runtime");
    }

    // Calculate distribution sizes
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_small");
    int base_size = array_size / num_procs;
    int remainder = array_size % num_procs;
    int my_size = base_size + (rank < remainder ? 1 : 0);
    local_array.resize(my_size);
    CALI_MARK_END("comp_small");
    CALI_MARK_END("comp");

    // Distribute data
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_small");
    std::vector<int> sendcounts(num_procs);
    std::vector<int> displs(num_procs);
    if (rank == root) {
        int offset = 0;
        for (i = 0; i < num_procs; i++) {
            sendcounts[i] = base_size + (i < remainder ? 1 : 0);
            displs[i] = offset;
            offset += sendcounts[i];
        }
    }

    MPI_Scatterv(input_array.data(), sendcounts.data(), displs.data(),
                MPI_INT, local_array.data(), my_size, MPI_INT,
                root, MPI_COMM_WORLD);
    CALI_MARK_END("comm_small");
    CALI_MARK_END("comm");

    // Local sort and sample selection
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_small");
    std::sort(local_array.begin(), local_array.end());

    // Improved sample selection
    int samples_per_proc = 3;  // Take multiple samples per process
    local_samples.resize(num_procs - 1);
    
    if (my_size > 0) {
        double stride = static_cast<double>(my_size) / num_procs;
        for (i = 0; i < num_procs - 1; i++) {
            int idx = static_cast<int>(stride * (i + 1));
            idx = std::min(idx, my_size - 1);
            local_samples[i] = local_array[idx];
        }
    } else {
        std::fill(local_samples.begin(), local_samples.end(), INT_MAX);
    }
    CALI_MARK_END("comp_small");
    CALI_MARK_END("comp");

    // Gather and process samples
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_small");
    if (rank == root) {
        global_samples.resize(num_procs * (num_procs - 1));
    }

    MPI_Gather(local_samples.data(), num_procs - 1, MPI_INT,
               global_samples.data(), num_procs - 1, MPI_INT,
               root, MPI_COMM_WORLD);

    if (rank == root) {
        std::sort(global_samples.begin(), global_samples.end());
        // Select evenly spaced splitters
        for (i = 0; i < num_procs - 1; i++) {
            int idx = ((i + 1) * global_samples.size()) / num_procs;
            local_samples[i] = global_samples[idx];
        }
    }

    // Broadcast final splitters
    MPI_Bcast(local_samples.data(), num_procs - 1, MPI_INT, root, MPI_COMM_WORLD);
    CALI_MARK_END("comm_small");
    CALI_MARK_END("comm");

    // Bucket distribution
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_large");
    std::vector<std::vector<int>> temp_buckets(num_procs);
    
    // Distribute elements to buckets
    for (const int& elem : local_array) {
        int bucket = get_bucket_index(elem, local_samples);
        temp_buckets[bucket].push_back(elem);
    }

    // Prepare for all-to-all
    int max_bucket_size = my_size + 1;  // Add 1 for count
    buckets.resize(max_bucket_size * num_procs);
    
    for (i = 0; i < num_procs; i++) {
        buckets[i * max_bucket_size] = temp_buckets[i].size();
        std::copy(temp_buckets[i].begin(), temp_buckets[i].end(),
                 buckets.begin() + i * max_bucket_size + 1);
    }
    CALI_MARK_END("comp_large");
    CALI_MARK_END("comp");

    // All-to-all exchange
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_large");
    bucket_buffer.resize(max_bucket_size * num_procs);
    MPI_Alltoall(buckets.data(), max_bucket_size, MPI_INT,
                 bucket_buffer.data(), max_bucket_size, MPI_INT,
                 MPI_COMM_WORLD);
    CALI_MARK_END("comm_large");
    CALI_MARK_END("comm");

    // Process received data
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_large");
    std::vector<int> merged_data;
    for (i = 0; i < num_procs; i++) {
        int count = bucket_buffer[i * max_bucket_size];
        merged_data.insert(merged_data.end(),
                         bucket_buffer.begin() + i * max_bucket_size + 1,
                         bucket_buffer.begin() + i * max_bucket_size + 1 + count);
    }
    
    std::sort(merged_data.begin(), merged_data.end());
    
    // Prepare for final gather
    local_bucket.resize(2 * base_size + 2);  // Add extra space for safety
    local_bucket[0] = merged_data.size();
    std::copy(merged_data.begin(), merged_data.end(), local_bucket.begin() + 1);
    CALI_MARK_END("comp_large");
    CALI_MARK_END("comp");

    // Final gather
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_large");
    if (rank == root) {
        output_buffer.resize(num_procs * (2 * base_size + 2));
    }
    
    MPI_Gather(local_bucket.data(), 2 * base_size + 2, MPI_INT,
               output_buffer.data(), 2 * base_size + 2, MPI_INT,
               root, MPI_COMM_WORLD);
    CALI_MARK_END("comm_large");
    CALI_MARK_END("comm");

    // Final assembly and verification
    if (rank == root) {
        CALI_MARK_BEGIN("finalization");
        std::vector<int> final_array;
        for (i = 0; i < num_procs; i++) {
            int count = output_buffer[i * (2 * base_size + 2)];
            final_array.insert(final_array.end(),
                             output_buffer.begin() + i * (2 * base_size + 2) + 1,
                             output_buffer.begin() + i * (2 * base_size + 2) + 1 + count);
        }
        
        // Copy to sorted_array
        assert(final_array.size() == array_size);
        std::copy(final_array.begin(), final_array.end(), sorted_array.begin());
        CALI_MARK_END("finalization");

        // Verify sort
        CALI_MARK_BEGIN("correctness_check");
        bool is_sorted = true;
        for (i = 1; i < array_size; i++) {
            if (sorted_array[i - 1] > sorted_array[i]) {
                is_sorted = false;
                std::cout << "Sort failed at position " << i << ": " 
                         << sorted_array[i-1] << " > " << sorted_array[i] << std::endl;
                break;
            }
        }
        
        if (is_sorted) {
            std::cout << "The array has been sorted correctly.\n";
        } else {
            std::cout << "The array is NOT sorted correctly.\n";
        }
        CALI_MARK_END("correctness_check");
    }

    // Metadata collection
    adiak::init(NULL);
    adiak::value("algorithm", "Sample Sort");
    adiak::value("programming_model", "MPI");
    adiak::launchdate();
    adiak::libraries();
    adiak::cmdline();
    adiak::clustername();
    adiak::value("data_type", "int");
    adiak::value("size_of_data_type", 4);
    adiak::value("input_size", array_size);
    adiak::value("input_type", input_type);
    adiak::value("num_procs", num_procs);
    adiak::value("scalability", "strong");
    adiak::value("group_num", 2);
    adiak::value("implementation_source", "online");

    MPI_Finalize();
    return 0;
}