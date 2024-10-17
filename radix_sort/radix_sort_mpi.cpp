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


int getMax(const std::vector<int>& arr) {
    return *std::max_element(arr.begin(), arr.end());
}

// Function to perform counting sort
void countingSort(std::vector<int>& input, std::vector<int>& buckets, int exp) {
    int n = input.size();
    std::vector<int> output(n, 0);
    int bucketCounts[10] = {0};
    // store occurences of the 10^exp digit
    for (int i = 0; i < input.size(); i++) {
        buckets[(input[i] / exp) % 10]++;
    }
    for (int i = 0; i < 10; i++) {
        bucketCounts[i] = buckets[i];
    }
    // Prefix sum to represent actual positions
    for (int i = 1; i < 10; i++) {
        bucketCounts[i] += bucketCounts[i - 1];
    }
    // Using calculated positions populate the output array
    for (int i = n - 1; i >= 0; i--) {
        int index = (input[i] / exp) % 10;
        output[bucketCounts[index] - 1] = input[i];
        bucketCounts[index]--;
    }
    for (int i = 0; i < n; i++) {
        input[i] = output[i];
    }
}

bool correctnessCheck(const std::vector<int>& unsorted, const std::vector<int>& sorted) {
    std::vector<int> sortCheck = unsorted;
    // std::cout << "Unsorted Array: ";
    // for (const auto& val : sortCheck) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;
    std::sort(sortCheck.begin(), sortCheck.end());
    // std::cout << "Sorted Copy Array: ";
    // for (const auto& val : sorted) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;
    if (sortCheck.size() != sorted.size()){
        return false;
    }
    for (size_t i = 0; i < sorted.size(); i++) {
        if (sortCheck[i] != sorted[i]) {
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "usage: ./radix_sort_mpi <array_size> <array_type>\n";
        return 0;
    }

    int array_size = atoi(argv[1]);
    int array_type = atoi(argv[2]);

    const char* whole_program = "main";
    const char* data_init_runtime = "data_init_runtime";
    const char* comm = "comm";
    const char* comm_large = "comm_large";
    const char* comm_small = "comm_small";
    const char* comp = "comp";
    const char* comp_small = "comp_small";
    const char* comp_large = "comp_large";
    const char* correctness_check = "correctness_check";

    int rank, numProc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    CALI_MARK_BEGIN(whole_program);
    
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<int> input;
    std::vector<int> adjustedInput;
    int n, 
        adjustedN,
        maxElement, 
        localArrSize, 
        remainderFiller;

    cali::ConfigManager mgr;
    mgr.start();

    CALI_MARK_BEGIN(data_init_runtime);

    if (rank == 0) {
        n = array_size;

        if (array_type == 1) {
            srand(time(NULL));
            for (int i = 0; i < array_size; i++) {
                input.push_back(rand() % array_size); // Random numbers between 0 and 999
            }
        } else if (array_type == 2) {
            for (int i = 0; i < array_size; i++) {
                input.push_back(i); // Sorted array
            }
        } else if (array_type == 3) {
            for (int i = array_size - 1; i >= 0; i--) {
                input.push_back(i); // Reverse sorted array
            }
        } else {
            printf("Invalid array_type. Please use 'random', 'sorted', or 'reverse'.\n");
            MPI_Finalize();
            return 0;
        }
        // std::cout << "Original array: ";
        // for (int num : input)
        //     std::cout << num << " ";
        // std::cout << std::endl;

        localArrSize = std::ceil(n / (double) numProc);
        remainderFiller = numProc * localArrSize - n;

        // rebuild array with adjusted size
        adjustedN = localArrSize * numProc;
        adjustedInput.resize(adjustedN);
        for (int i = 0; i < adjustedN; i++) {
            if (i < remainderFiller) {
                adjustedInput[i] = 0;
            } else {    
                adjustedInput[i] = input[i - remainderFiller];
            }
        }

        maxElement = *std::max_element(input.begin(), input.end());
    }

    CALI_MARK_END(data_init_runtime);

    CALI_MARK_BEGIN(comm);
    CALI_MARK_BEGIN(comm_small);

    // broadcasts shared vars
    MPI_Bcast(&maxElement, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&localArrSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&adjustedN, 1, MPI_INT, 0, MPI_COMM_WORLD);

    CALI_MARK_END(comm_small);

    MPI_Barrier(MPI_COMM_WORLD);

    // allocate local array
    std::vector<int> local_arr(localArrSize, 0);

    CALI_MARK_BEGIN(comm_small);
    MPI_Scatter(adjustedInput.data(), localArrSize, MPI_INTEGER, local_arr.data(), localArrSize, MPI_INTEGER, 0, MPI_COMM_WORLD);
    CALI_MARK_END(comm_small);
    
    std::vector<int> allCounts(10 * numProc, 0);
    std::vector<int> newlocal_arr(localArrSize);

    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_large);

    // call counting sort for each digit in max element
    for (int exp = 1; maxElement / exp > 0; exp *= 10) {
        std::vector<int> count(10, 0);
        std::vector<int> totalCountsSum(10, 0);
        std::vector<int> prefixSum(10, 0);
        std::vector<int> totalCountsSumLeft(10, 0);

        countingSort(local_arr, count, exp);
        MPI_Allgather(count.data(), 10, MPI_INTEGER, allCounts.data(), 10, MPI_INTEGER, MPI_COMM_WORLD);

        // populate aggreated data from all processes
        for (int i = 0; i < 10 * numProc; i++) {
            int lsd = i % 10;
            int currProc = i / 10;
            int val = allCounts[i];
            if (currProc < rank) {
                totalCountsSumLeft[lsd] += val;
            }
            totalCountsSum[lsd] += val;
            prefixSum[lsd] += val;
        }

        for (int i = 1; i < 10; i++) {
            prefixSum[i] += prefixSum[i - 1];
        }

        MPI_Request request;
        MPI_Status status;

        int lsdSent[10] = {0};
        int valIndexPair[2];
        int destIndex, destProcess, localDestIndex;
        // Move keys within each processor to appropriate buckets (local operation) 
        // Send/receive keys between the processors
        CALI_MARK_BEGIN(comm_small);
        for (int i = 0; i < localArrSize; i++) {
            int val = local_arr[i];
            int lsd = (val / exp) % 10;
            destIndex = prefixSum[lsd] - totalCountsSum[lsd] + totalCountsSumLeft[lsd] + lsdSent[lsd];
            lsdSent[lsd]++;
            destProcess = destIndex / localArrSize;

            valIndexPair[0] = val;
            valIndexPair[1] = destIndex;

            MPI_Isend(valIndexPair, 2, MPI_INT, destProcess, 0, MPI_COMM_WORLD, &request);
            MPI_Recv(valIndexPair, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            localDestIndex = valIndexPair[1] % localArrSize;
            val = valIndexPair[0];
            newlocal_arr[localDestIndex] = val;
        }
        CALI_MARK_END(comm_small);

        local_arr = newlocal_arr; 
    }

    CALI_MARK_END(comp_large);
    CALI_MARK_END(comp);


    int * output;
    if (rank == 0) {
        output = new int[adjustedN]{};
    }

    CALI_MARK_BEGIN(comm_large);
    MPI_Gather(local_arr.data(), localArrSize, MPI_INTEGER, &output[rank * localArrSize], localArrSize, MPI_INTEGER, 0, MPI_COMM_WORLD);
    CALI_MARK_END(comm_large);

    CALI_MARK_END(comm);


    CALI_MARK_BEGIN(correctness_check);

    if (rank == 0) {
        output += remainderFiller; 
        // correctness check
        std::vector<int> oput(output, output + adjustedN);
        if (!correctnessCheck(input, oput)) {
            return 1;
        }
        output -= remainderFiller;
        delete[] output;
    }
    CALI_MARK_END(correctness_check);

    CALI_MARK_END(whole_program);

    adiak::init(NULL);
    adiak::launchdate();    // launch date of the job
    adiak::libraries();     // Libraries used
    adiak::cmdline();       // Command line used to launch the job
    adiak::clustername();   // Name of the cluster
    adiak::value("algorithm", "merge"); // The name of the algorithm you are using (e.g., "merge", "bitonic")
    adiak::value("programming_model", "mpi"); // e.g. "mpi"
    adiak::value("data_type", "int"); // The datatype of input elements (e.g., double, int, float)
    adiak::value("size_of_data_type", sizeof(int)); // sizeof(datatype) of input elements in bytes (e.g., 1, 2, 4)
    adiak::value("array_size", array_size); // The number of elements in input dataset (1000)
    adiak::value("array_type", array_type); // For sorting, this would be choices: ("Sorted", "ReverseSorted", "Random", "1_perc_perturbed")
    adiak::value("numProv", numProc); // The number of processors (MPI ranks)
    adiak::value("scalability", "strong"); // The scalability of your algorithm. choices: ("strong", "weak")
    adiak::value("group_num", 2); // The number of your group (integer, e.g., 1, 10)
    adiak::value("implementation_source", "online"); // Where you got the source code of your algorithm. choices: ("online", "ai", "handwritten").

    mgr.stop();
    mgr.flush();


    MPI_Finalize();
}
