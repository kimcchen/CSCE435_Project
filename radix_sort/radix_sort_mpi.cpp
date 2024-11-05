#include <vector>
#include <algorithm>
#include <mpi.h>
#include <iostream>
#include <random>
#include <string>
#include <climits>


#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

namespace {
    const char* whole_program = "main";
    const char* data_init_runtime = "data_init_runtime";
    const char* comm = "comm";
    const char* comm_large = "comm_large";
    const char* comm_small = "comm_small";
    const char* comp = "comp";
    const char* comp_small = "comp_small";
    const char* comp_large = "comp_large";
    const char* correctness_check = "correctness_check";
}


int getMax(const std::vector<int>& vec) {
    if (vec.empty()) {
        return 0;  
    }
    return *std::max_element(vec.begin(), vec.end());
}


void countingSort(std::vector<int>& array, const std::vector<int>& rec_buf, 
                 int digit, int num_process, int rank) {
    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_small);
    std::vector<int> local_count(10, 0);
    
    for (int num : rec_buf) {
        local_count[(num / digit) % 10]++;
    }

    CALI_MARK_END(comp_small);
    CALI_MARK_END(comp);

    // Reduce all the sub counts to root process
    if (rank == 0) {
        std::vector<int> count(10, 0);
        CALI_MARK_BEGIN(comm);
        CALI_MARK_BEGIN(comm_small);
        MPI_Reduce(local_count.data(), count.data(), 10, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        CALI_MARK_END(comm_small);
        CALI_MARK_END(comm);

        CALI_MARK_BEGIN(comp);
        CALI_MARK_BEGIN(comp_large);
        // Calculate cumulative count
        for (int i = 1; i < 10; i++) {
            count[i] += count[i - 1];
        }

       

        std::vector<int> temp_array(array.size());
        for (int i = array.size() - 1; i >= 0; i--) {
            temp_array[count[(array[i] / digit) % 10] - 1] = array[i];
            count[(array[i] / digit) % 10]--;
        }
        array = std::move(temp_array);  // Move assignment instead of memcpy
        CALI_MARK_END(comp_large);
        CALI_MARK_END(comp);
    } else {
        CALI_MARK_BEGIN(comm);
        CALI_MARK_BEGIN(comm_small);
        MPI_Reduce(local_count.data(), nullptr, 10, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        CALI_MARK_END(comm_small);
        CALI_MARK_END(comm);
    }
}


void radix_sort(std::vector<int>& array, int num_process, int rank) {
    int n = array.size();
    int rem = n % num_process;  // elements remaining after division among processes
    int dim, displacement;

    if (rank < rem) {
        dim = n/num_process + 1;
        displacement = rank * dim;
    } else {
        dim = n/num_process;
        displacement = rank * dim + rem;
    }

    std::vector<int> rec_buf(dim);
    std::vector<int> sendcounts;
    std::vector<int> displs;

    if (rank == 0) {
        sendcounts.resize(num_process);
        displs.resize(num_process);
    }

    CALI_MARK_BEGIN(comm);
    CALI_MARK_BEGIN(comm_large);

    MPI_Gather(&dim, 1, MPI_INT, sendcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&displacement, 1, MPI_INT, displs.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(array.data(), sendcounts.data(), displs.data(), MPI_INT,
                 rec_buf.data(), dim, MPI_INT, 0, MPI_COMM_WORLD);

    int local_max = getMax(rec_buf);
    int global_max;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    CALI_MARK_END(comm_large);
    CALI_MARK_END(comm);



    for (int digit = 1; global_max / digit > 0; digit *= 10) {
        countingSort(array, rec_buf, digit, num_process, rank);
    }
}

// Function to generate test data based on array type
std::vector<int> generate_data(int size, const std::string& array_type) {
    std::vector<int> data(size);
    std::random_device rd;
    std::mt19937 gen(rd());

    if (array_type == "Random") {
        // Generate random numbers between 0 and size*10
        std::uniform_int_distribution<> dis(0, INT_MAX);
        for (int i = 0; i < size; i++) {
            data[i] = dis(gen);
        }
    } else if (array_type == "Sorted") {
        // Generate sorted array
        for (int i = 0; i < size; i++) {
            data[i] = i;
        }
    } else if (array_type == "Reverse") {
        // Generate reverse sorted array
        for (int i = 0; i < size; i++) {
            data[i] = size - i;
        }
    } else if (array_type == "OnePrecPret") {
        for (int i = 0; i < size; i++) {
            data[i] = i; // Assign sorted values
        }
        int num_swaps = std::max(1, size / 100);
        for(int i = 0; i < num_swaps; i++) {
            int swap_pos = rand() % (size - 1);  
            std::swap(data[swap_pos], data[swap_pos + 1]);
        }
    } else {
        // Default to random if type is not recognized
        std::uniform_int_distribution<> dis(0, INT_MAX);
        for (int i = 0; i < size; i++) {
            data[i] = dis(gen);
        }
    }
    return data;
}

// Function to verify if array is sorted
bool is_sorted(const std::vector<int>& arr) {
    return std::is_sorted(arr.begin(), arr.end());
}

// Function to print array (for debugging with small arrays)
void print_array(const std::vector<int>& arr, const std::string& message) {
    std::cout << message << ": ";
    for (int num : arr) {
        std::cout << num << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    int rank, num_process;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    if (argc != 3) {
        if (rank == 0) {
            std::cout << "usage: ./radix_sort_mpi <array_size> <array_type>\n";
            std::cout << "array_type options: random, sorted, reverse, equal\n";
        }
        MPI_Finalize();
        return 1;
    }
    

    int array_size = std::stoi(argv[1]);
    std::string array_type = argv[2];

    std::vector<int> array;

    cali::ConfigManager mgr;
    mgr.start();

    double start_time = MPI_Wtime();
    CALI_MARK_BEGIN(whole_program);
    CALI_MARK_BEGIN(data_init_runtime);

    // Only root process generates the initial array
    if (rank == 0) {
        array = generate_data(array_size, array_type);
        
        // // Print first few elements if array is large
        // if (array_size <= 20) {
        //     print_array(array, "Initial array");
        // } else {
        //     std::cout << "Array size: " << array_size << ", Type: " << array_type << std::endl;
        //     std::cout << "First 10 elements: ";
        //     for (int i = 0; i < 10 && i < array_size; i++) {
        //         std::cout << array[i] << " ";
        //     }
        //     std::cout << std::endl;
        // }
    } else {
        array.resize(array_size);  // Non-root processes need space for receiving
    }
    CALI_MARK_END(data_init_runtime);

    // Perform radix sort
    radix_sort(array, num_process, rank);

    // End timing
    double end_time = MPI_Wtime();

    // Verify and print results
    if (rank == 0) {
        // if (array_size <= 20) {
        //     print_array(array, "Sorted array");
        // } else {
        //     std::cout << "First 10 elements after sorting: ";
        //     for (int i = 0; i < 10 && i < array_size; i++) {
        //         std::cout << array[i] << " ";
        //     }
        //     std::cout << std::endl;
        CALI_MARK_BEGIN(correctness_check);
        bool sorted = is_sorted(array);
        CALI_MARK_END(correctness_check);
        std::cout << "Array is " << (sorted ? "correctly" : "not") << " sorted" << std::endl;
        std::cout << "Time taken: " << (end_time - start_time) << " seconds" << std::endl;
        std::cout << "Number of processes: " << num_process << std::endl;
        std::cout << "Array_size: " << array_size << std::endl;
        std::cout << "Array_type: " << array_type << std::endl;
    }


    adiak::init(NULL);
    adiak::launchdate();    // launch date of the job
    adiak::libraries();     // Libraries used
    adiak::cmdline();       // Command line used to launch the job
    adiak::clustername();   // Name of the cluster
    adiak::value("algorithm", "radix"); // The name of the algorithm you are using (e.g., "merge", "bitonic")
    adiak::value("programming_model", "mpi"); // e.g. "mpi"
    adiak::value("data_type", "int"); // The datatype of input elements (e.g., double, int, float)
    adiak::value("size_of_data_type", sizeof(int)); // sizeof(datatype) of input elements in bytes (e.g., 1, 2, 4)
    adiak::value("array_size", array_size); // The number of elements in input dataset (1000)
    adiak::value("array_type", array_type); // For sorting, this would be choices: ("Sorted", "ReverseSorted", "Random", "1_perc_perturbed")
    adiak::value("num_procs", num_process); // The number of processors (MPI ranks)
    adiak::value("scalability", "strong"); // The scalability of your algorithm. choices: ("strong", "weak")
    adiak::value("group_num", 2); // The number of your group (integer, e.g., 1, 10)
    adiak::value("implementation_source", "online"); // Where you got the source code of your algorithm. choices: ("online", "ai", "handwritten").

    mgr.stop();
    mgr.flush();

    MPI_Finalize();
    return 0;
}
