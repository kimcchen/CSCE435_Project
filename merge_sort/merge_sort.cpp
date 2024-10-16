#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <vector>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>
#include <iostream>
#include <random>
#include <algorithm>
#include <math.h>

std::vector<int> merge(std::vector<int> arr1, std::vector<int> arr2){
    size_t size_1 = arr1.size();
    size_t size_2 = arr2.size();
    std::vector<int> final_arr(size_1 + size_2);
    int i = 0;
    int j = 0;
    int k = 0;

    while ((i < size_1) && (j < size_2)){
        if (arr1[i] < arr2[j]){
            final_arr.at(k) = arr1[i];
            i += 1;
        }
        else{
            final_arr.at(k) = arr2[j];
            j += 1;
        }
        k += 1;
    }
    
    // copy leftover elements
    while (i < size_1){
        final_arr.at(k) = arr1[i];
        i += 1;
        k += 1;
    }
        
    while (j < size_2){
        final_arr.at(k) = arr2[j];
        j += 1;
        k += 1;
    }
        
    return final_arr;
}


std::vector<int> merge_sort(std::vector<int> arr){
    // base case
    if (arr.size() <= 1){
        return arr;
    }

    // find midpoint
    int mid = arr.size() / 2;

    // generate two separate subarrays
    std::vector<int> left(arr.begin(), arr.begin() + mid);
    std::vector<int> right(arr.begin() + mid, arr.end());

    // merge the sorted halves
    return merge(merge_sort(left), merge_sort(right));
}


// Function to print a vector
void printVector(const std::vector<int>& arr) {
    for (int num : arr) {
        std::cout << num << " ";
    }
    std::cout << std::endl;
}

int main (int argc, char *argv[]) {
    CALI_CXX_MARK_FUNCTION;

    // get cmd arguments
    int arr_size = atoi(argv[1]);
    int arr_type = atoi(argv[2]);

    // initialize MPI
    int rank, num_procs;
    const char *input_type;
    

    // define caliper region names
    // const char* whole_program = "main";
    const char* data_init_runtime = "data_init_runtime";
    const char* comm = "comm";
    const char* comm_large = "comm_large";
    const char* comp = "comp";
    const char* comp_large = "comp_large";
    const char* correctness_check = "correctness_check";

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);

    cali::ConfigManager mgr;
    mgr.start();

    // DATA INITIALIZATION REGION
    CALI_MARK_BEGIN(data_init_runtime);

    int worker_size = arr_size / num_procs;
    std::vector<int> worker_arr(worker_size);

    if (arr_type == 0) {
        int start = rank * worker_size;
        for (int i = 0; i < worker_size; i++) {
            worker_arr[i] = start + i;
        }
        input_type = "Sorted";
    }
    else if (arr_type == 1) {
        int start = rank * worker_size;
        for (int i = 0; i < worker_size; i++) {
            worker_arr[i] = rand() % arr_size;
        }
        input_type = "Random";
    }
    else if (arr_type == 2) {
        int start_value = arr_size - (rank * worker_size);
        // int start_idx = rank * worker_size;
        for (int i = 0; i < worker_size; i++) {
            // arr[start_idx + i] = start_value - i;
            worker_arr[i] = start_value - i;
        }
        input_type = "ReverseSorted";
    }
    else if (arr_type == 3) {
        int start = rank * worker_size;
        for (int i = 0; i < worker_size; i++) {
            worker_arr[i] = start + i;
        }

        int numPertubed = std::ceil(worker_size*0.01);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist(0, worker_size - 1);
        for (int i = 0; i < numPertubed; ++i) {
            int idx_1 = dist(gen);
            int idx_2 = dist(gen);

            while (idx_1 == idx_2) {
                idx_2 = dist(gen);
            }

            std::swap(worker_arr[idx_1], worker_arr[idx_2]);
        }

        input_type = "1_perc_perturbed";
    }

    CALI_MARK_END(data_init_runtime);
    

    std::vector<int> arr(arr_size);

    // sort local copy
    // COMPUTATION REGION
    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_large);
    std::vector<int> sorted_arr = merge_sort(worker_arr);
    CALI_MARK_END(comp_large);
    CALI_MARK_END(comp);

    // merge sorted pieces back together
    CALI_MARK_BEGIN(comm);
    CALI_MARK_BEGIN(comm_large);
    int step = 1;
    while (step < num_procs) {
        if (rank % (2 * step) == 0) {
            if (rank + step < num_procs) {
                // receive sorted array from another process
                int received_size = worker_size;
                std::vector<int> received_array(received_size);
                MPI_Recv(received_array.data(), received_size, MPI_INT, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // merge the received array with the local array
                sorted_arr = merge(sorted_arr, received_array);
                worker_size += received_size;
            }
        } 
        else {
            // send the sorted array to another process
            int target = rank - step;
            MPI_Send(sorted_arr.data(), worker_size, MPI_INT, target, 0, MPI_COMM_WORLD);
            break;
        }
        
        step *= 2;
    }
    CALI_MARK_END(comm_large);
    CALI_MARK_END(comm);

    // gather the final sorted array
    if (rank == 0) {
        printf("Sorted array: ");
        printVector(sorted_arr);
        
        
        CALI_MARK_BEGIN(correctness_check);
        bool correct = true;
        for (int i = 1; i < arr_size; i++) {
            if (arr[i] < arr[i - 1]) {
                correct = false;
                break;
            }
        }
        printf("Final array is sorted: %s\n", correct ? "true" : "false");
        CALI_MARK_END(correctness_check);
    }


    adiak::init(NULL);
    adiak::launchdate();    // launch date of the job
    adiak::libraries();     // Libraries used
    adiak::cmdline();       // Command line used to launch the job
    adiak::clustername();   // Name of the cluster
    adiak::value("algorithm", "merge"); // The name of the algorithm you are using (e.g., "merge", "bitonic")
    adiak::value("programming_model", "mpi"); // e.g. "mpi"
    adiak::value("data_type", "int"); // The datatype of input elements (e.g., double, int, float)
    adiak::value("size_of_data_type", sizeof(int)); // sizeof(datatype) of input elements in bytes (e.g., 1, 2, 4)
    adiak::value("input_size", arr_size); // The number of elements in input dataset (1000)
    adiak::value("input_type", input_type); // For sorting, this would be choices: ("Sorted", "ReverseSorted", "Random", "1_perc_perturbed")
    adiak::value("num_procs", num_procs); // The number of processors (MPI ranks)
    adiak::value("scalability", "strong"); // The scalability of your algorithm. choices: ("strong", "weak")
    adiak::value("group_num", 2); // The number of your group (integer, e.g., 1, 10)
    adiak::value("implementation_source", "online"); // Where you got the source code of your algorithm. choices: ("online", "ai", "handwritten").


    // Flush Caliper output before finalizing MPI
    mgr.stop();
    mgr.flush();

    MPI_Finalize();
    return 0;
}
