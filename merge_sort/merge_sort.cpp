// #include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <vector>

// #include <caliper/cali.h>
// #include <caliper/cali-manager.h>
// #include <adiak.hpp>
#include <iostream>
#include <algorithm>
#include <math.h>

std::vector<int> merge(std::vector<int> arr1, std::vector<int> arr2){
    size_t size_1 = arr1.size();
    size_t size_2 = arr2.size();
    std::vector<int> final_arr;
    int i = 0;
    int j = 0;

    while ((i < size_1) && (j < size_2)){
        if (arr1[i] < arr2[j]){
            final_arr.push_back(arr1[i]);
            i += 1;
        }
        else{
            final_arr.push_back(arr2[j]);
            j += 1;
        }
    }
    
    // copy leftover elements
    while (i < size_1){
        final_arr.push_back(arr1[i]);
        i += 1;
    }
        
    while (j < size_2){
        final_arr.push_back(arr2[j]);
        j += 1;
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

    // sort subarrays
    std::vector<int> left_sorted = merge_sort(left);
    std::vector<int> right_sorted = merge_sort(right);

    // Merge the sorted halves
    return merge(left_sorted, right_sorted);
}


// Function to print a vector
void printVector(const std::vector<int>& arr) {
    for (int num : arr) {
        std::cout << num << " ";
    }
    std::cout << std::endl;
}

int main() {
    std::vector<int> arr = {38, 27, 43, 3, 9, 82, 10};

    std::cout << "Original array: ";
    printVector(arr);

    std::vector<int> sortedArr = merge_sort(arr);

    std::cout << "Sorted array: ";
    printVector(sortedArr);

    return 0;
}
