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

int main(int argc, char *argv[]) {
    // Mark the entire main function for Caliper
    CALI_CXX_MARK_FUNCTION;

    // Initialize Adiak metadata collection
    adiak::init(NULL);
    adiak::value("algorithm", "Sample Sort");
    adiak::value("programming_model", "MPI");
    adiak::launchdate();
    adiak::libraries();
    adiak::cmdline();
    adiak::clustername();

    // MPI variables
    int num_procs, rank, root = 0;
    int i, j, k, num_bucket_elements, num_elements, count;
    std::vector<int> input_array, local_array, local_samples, global_samples, buckets, bucket_buffer, local_bucket, sorted_array, output_buffer;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Handle command-line arguments
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

    // Step 1: Initialize and populate the array on the root process
    if (rank == root) {
        CALI_MARK_BEGIN("data_init_runtime");  // Caliper mark for data initialization
        input_array.resize(array_size);  // Allocate space for the input array
        sorted_array.resize(array_size); // Allocate space for the sorted output array
        initialize_array(input_array, input_type);  // Initialize the array
        CALI_MARK_END("data_init_runtime");
    }

    // Broadcast array size to all processes
    MPI_Bcast(&array_size, 1, MPI_INT, root, MPI_COMM_WORLD);
    num_bucket_elements = array_size / num_procs;  // Number of elements for each process
    local_array.resize(num_bucket_elements);  // Allocate space for local array

    // High-level communication region
    CALI_MARK_BEGIN("comm");

    // Step 2: Scatter the input array to all processes
    CALI_MARK_BEGIN("comm_small");  // Caliper mark for small communication (scatter)
    MPI_Scatter(input_array.data(), num_bucket_elements, MPI_INT, local_array.data(), 
                num_bucket_elements, MPI_INT, root, MPI_COMM_WORLD);  // Scatter data
    CALI_MARK_END("comm_small");

    // High-level computation region
    CALI_MARK_BEGIN("comp");

    // Step 3: Sort local data
    CALI_MARK_BEGIN("comp_small");  // Caliper mark for small computation (local sort)
    std::sort(local_array.begin(), local_array.end());
    CALI_MARK_END("comp_small");

    // Step 4: Select samples from the sorted local data
    local_samples.resize(num_procs - 1);
    for (i = 0; i < (num_procs - 1); i++) {
        local_samples[i] = local_array[num_bucket_elements / (num_procs * num_procs) * (i + 1)];
    }

    // Step 5: Gather all samples to the root process
    global_samples.resize(num_procs * (num_procs - 1));
    CALI_MARK_BEGIN("comm_small");  // Caliper mark for small communication (gather)
    MPI_Gather(local_samples.data(), num_procs - 1, MPI_INT, global_samples.data(), num_procs - 1, 
               MPI_INT, root, MPI_COMM_WORLD);  // Gather samples
    CALI_MARK_END("comm_small");

    // Step 6: Root process sorts the samples and selects splitters
    if (rank == root) {
        CALI_MARK_BEGIN("comp_large");  // Caliper mark for large computation (sorting splitters)
        std::sort(global_samples.begin(), global_samples.end());
        for (i = 0; i < num_procs - 1; i++) {
            local_samples[i] = global_samples[(num_procs - 1) * (i + 1)];
        }
        CALI_MARK_END("comp_large");
    }

    // Step 7: Broadcast splitters to all processes
    CALI_MARK_BEGIN("comm_small");  // Caliper mark for small communication (broadcast)
    MPI_Bcast(local_samples.data(), num_procs - 1, MPI_INT, root, MPI_COMM_WORLD);
    CALI_MARK_END("comm_small");

    // Step 8: Redistribute local data based on splitters
    buckets.resize(array_size + num_procs);
    j = 0;
    k = 1;
    for (i = 0; i < num_bucket_elements; i++) {
        if (j < (num_procs - 1)) {
            if (local_array[i] < local_samples[j]) {
                buckets[((num_bucket_elements + 1) * j) + k++] = local_array[i];
            } else {
                buckets[(num_bucket_elements + 1) * j] = k - 1;
                k = 1;
                j++;
                i--;
            }
        } else {
            buckets[((num_bucket_elements + 1) * j) + k++] = local_array[i];
        }
    }
    buckets[(num_bucket_elements + 1) * j] = k - 1;

    // Step 9: Send and receive buckets from other processes
    bucket_buffer.resize(array_size + num_procs);
    CALI_MARK_BEGIN("comm_large");  // Caliper mark for large communication (all-to-all)
    MPI_Alltoall(buckets.data(), num_bucket_elements + 1, MPI_INT, bucket_buffer.data(), 
                 num_bucket_elements + 1, MPI_INT, MPI_COMM_WORLD);  // All-to-all communication
    CALI_MARK_END("comm_large");

    // Step 10: Rearrange the received buckets
    local_bucket.resize(2 * array_size / num_procs);
    count = 1;
    for (j = 0; j < num_procs; j++) {
        k = 1;
        for (i = 0; i < bucket_buffer[(array_size / num_procs + 1) * j]; i++) {
            local_bucket[count++] = bucket_buffer[(array_size / num_procs + 1) * j + k++];
        }
    }
    local_bucket[0] = count - 1;

    // Step 11: Sort local data again after redistribution
    num_elements = local_bucket[0];
    CALI_MARK_BEGIN("comp_large");  // Caliper mark for large computation (final sorting)
    std::sort(local_bucket.begin() + 1, local_bucket.begin() + num_elements + 1);
    CALI_MARK_END("comp_large");

    // End the high-level computation region
    CALI_MARK_END("comp");

    // Step 12: Gather sorted data back to the root process
    output_buffer.resize(2 * array_size);
    if (rank == root) {
        sorted_array.resize(array_size);
    }
    CALI_MARK_BEGIN("comm_large");  // Caliper mark for large communication (gather)
    MPI_Gather(local_bucket.data(), 2 * num_bucket_elements, MPI_INT, output_buffer.data(), 
               2 * num_bucket_elements, MPI_INT, root, MPI_COMM_WORLD);
    CALI_MARK_END("comm_large");

    // End the high-level communication region
    CALI_MARK_END("comm");

    // Step 13: Rearrange gathered sorted data
    if (rank == root) {
        count = 0;
        for (j = 0; j < num_procs; j++) {
            k = 1;
            for (i = 0; i < output_buffer[(2 * array_size / num_procs) * j]; i++) {
                sorted_array[count++] = output_buffer[(2 * array_size / num_procs) * j + k++];
            }
        }

        // Step 14: Correctness check
        CALI_MARK_BEGIN("correctness_check");  // Caliper mark for correctness check
        bool is_sorted = true;
        for (i = 1; i < array_size; i++) {
            if (sorted_array[i - 1] > sorted_array[i]) {
                is_sorted = false;
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

    // Finalize MPI and Caliper
    MPI_Finalize();
    return 0;
}