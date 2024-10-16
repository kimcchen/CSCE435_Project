
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>
#include <vector>
#include <algorithm>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

using namespace std;

// #define N 128
// #define MAX_ITEM (3*N)
#define ROOT_RANK 0
#define MAX_OUTPUT 1000
#define INCREASING 0
#define DECREASING 1

// Comparison function used by qsort to sort integers in ascending order
int compare_ints_up(int arg1, int arg2) 
{
    return arg1 < arg2;
}

// Comparison function used by qsort to sort integers in descending order
int compare_ints_down(int arg1, int arg2) 
{
    return arg1 > arg2;
}

// Sort the array in place.
// If dir is 0, then sort the numbers in ascending order. Otherwise, sort them
// in descending order
void qsort_ints(vector<int>& array, int dir)
{
    if (dir == INCREASING)
        sort(array.begin(), array.end(), compare_ints_up);
    else
        sort(array.begin(), array.end(), compare_ints_down);
}

// Print up to MAX_OUTPUT elements in an array on one line.
void print_array(const vector<int>& array)
{
    for (size_t i = 0; i < array.size(); i++)
    {
        cout << array[i] << " ";
    }
    cout << endl;
}

void compare_arrays(const vector<int>& a1, const vector<int>& a2)
{
    if (a1 == a2)
        cout << "The comparison test succeeded" << endl;
    else
        cout << "The arrays don't match!" << endl;
} // end compare_arrays

vector<int> create_random_int_array(int len)
{
    vector<int> ret(len);
    srand(1);
    for (int i = 0; i < len; i++)
        ret[i] = rand() % 1000;
    return ret;
} // create_random_int_array

// Return the value of the bit_positionth bit in the_int.
// This method assumes bit_position will be between 0 and 30.
// The bits are ordered from right to left, so
// The bit at position 0 of 4 is 0 because its binary representation is 100.
// The bit at position 1 of 4 is 0.
// The bit at position 2 of 4 is 1. 
// All remaining bits of 4 are 0.
int bit(int the_int, int bit_position)
{
  // shift the bits right bit_position positions.
  return (the_int >> bit_position) % 2;
}

// Return the base 2 logarithm of the integer.
int log2i(int the_int)
{
  // just do some type casting
  return (int)log2((double)the_int);
}

// Return an integer that is the_int, but with the bit_positionth bit
// flipped.
// For example, flipbut(12,3) will return 4.
int flipbit(int the_int, int bit_position)
{
  int one = 1;
  //                shift a 1-bit into the correct position
  //    then xor it with the_int
  return the_int ^ (one << bit_position);
}

// Function to merge data in the increasing order
void merge_up(vector<int>& my_array, int rank, const vector<int>& parray, int prank) {
    for (int i = 0; i < my_array.size(); i++) {
        if (rank < prank) {
            if (my_array[i] >= parray[i])
                my_array[i] = parray[i];
        } else {
            if (my_array[i] < parray[i])
                my_array[i] = parray[i];
        }
    }
}

// Function to merge data in the decreasing order
void merge_down(vector<int>& my_array, int rank, const vector<int>& parray, int prank) {
    for (int i = 0; i < my_array.size(); i++) {
        if (rank < prank) {
            if (my_array[i] < parray[i])
                my_array[i] = parray[i];
        } else {
            if (my_array[i] >= parray[i])
                my_array[i] = parray[i];
        }
    }
}

int main(int argc, char **argv)
{
	const int TAG = 1;
    MPI_Status status;
	int  numprocs, rank, rc;
	vector<int> array, my_array, parray, sorted_array;
	int length, length_per_process;
    
    /* Define Caliper region names */
    const char *whole_computation = "whole_computation";  // For the entire computation process
    const char *data_init_runtime = "data_init_runtime";  // For data initialization during runtime
    const char *comm = "comm";                            // For all communication-related operations
    const char *comm_small = "comm_small";                // For small communication (e.g., MPI_Sendrecv)
    const char *comm_large = "comm_large";                // For large communication (e.g., MPI_Scatter/Gather)
    const char *comp = "comp";                            // For all computation-related operations
    const char *comp_small = "comp_small";                // For small data computation (e.g., local sorting)
    const char *comp_large = "comp_large";                // For large data computation (e.g., merging arrays)
    const char *correctness_check = "correctness_check";         // Encompasses the entire sorting process
  
	/******* Initialize *******/
	rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Initialize Caliper
    cali::ConfigManager mgr;
    mgr.start();

    // double whole_computation_start = MPI_Wtime();
    CALI_MARK_BEGIN(whole_computation);

    // Get the matrix size from argv[1], and ensure argv[1] is passed.
    if (argc < 2) {
        if (rank == ROOT_RANK) {
            cout << "Matrix size not provided. Terminating." << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    CALI_MARK_BEGIN(data_init_runtime);  
    int N = atoi(argv[1]); // Convert the first argument to an integer, which is the matrix size
    // cout << "Matrix size: " << N << endl;
  
  int exponent = log2i(numprocs);
  if (numprocs < 2 || pow(2,exponent) != numprocs)
  {
		cout << "The number of processors must be greater than 1 and must be a power of 2" << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
  }
  
	length_per_process = N/numprocs;
  // if the number of items per proc is not even, then reduce it 
  // so it is even
  if (length_per_process%2 != 0)
    length_per_process--;
  // If the number of items is not evenly divisible by the processors,
  // the reduce it so that it is.
  length = length_per_process*numprocs;
  
  // Allocate space for my initial chunk of the array.
	my_array.resize(length_per_process);
	if (rank == ROOT_RANK)
    {
		// cout << "The number of elements in the array is" << length << endl;
        array = create_random_int_array(length);

        // Print the unsorted global array
        // cout << "The unsorted global array (before scattering):" << endl;
        // print_array(array);
	}
    CALI_MARK_END(data_init_runtime);  
  
  /***** Do that sorting thing *******/
    CALI_MARK_BEGIN(comm);
    CALI_MARK_BEGIN(comm_large);
    MPI_Scatter(array.data(), length_per_process, MPI_INT, my_array.data(), length_per_process, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    CALI_MARK_END(comm_large);
    CALI_MARK_END(comm);

    CALI_MARK_BEGIN(comp);

    // Sort locally
    // qsort_ints(my_array, INCREASING);
    CALI_MARK_BEGIN(comp_small);
    qsort_ints(my_array, bit(rank, 0) == 0 ? INCREASING : DECREASING);
    CALI_MARK_END(comp_small);
    // Perform the bitonic merge
    for (int depth = 1; depth <= log2(numprocs); depth++) {
        for (int p = depth - 1; p >= 0; p--) {
            int prank = flipbit(rank, p);
            parray.resize(length_per_process);

            // Exchange data with partner
            CALI_MARK_BEGIN(comm);
            CALI_MARK_BEGIN(comm_small);
            MPI_Sendrecv(my_array.data(), length_per_process, MPI_INT, prank, TAG,
                         parray.data(), length_per_process, MPI_INT, prank, TAG,
                         MPI_COMM_WORLD, &status);
            CALI_MARK_END(comm_small);
            CALI_MARK_END(comm);

            CALI_MARK_BEGIN(comp_large);
            if (bit(rank, depth) == 0) {
                merge_up(my_array, rank, parray, prank);
            } else {
                merge_down(my_array, rank, parray, prank);
            }
            CALI_MARK_END(comp_large);
        }
        // Local sort after each depth, based on bit value of rank at depth
        CALI_MARK_BEGIN(comp_small);
        qsort_ints(my_array, bit(rank, depth) == 0 ? INCREASING : DECREASING);
        CALI_MARK_END(comp_small);  
    }
    CALI_MARK_END(comp);

     // Gather the results
    CALI_MARK_BEGIN(comm);
    CALI_MARK_BEGIN(comm_large);
    if (rank == ROOT_RANK) {
        // Sort locally as a reference
        // qsort_ints(array, INCREASING);

        sorted_array.resize(length);  // Allocate space for sorted array

        // Compare the gathered result with the reference sorted array
        // compare_arrays(sorted_array, array);
    }
    // Gather the results
    MPI_Gather(my_array.data(), length_per_process, MPI_INT, sorted_array.data(), length_per_process, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    CALI_MARK_END(comm_large);
    CALI_MARK_END(comm);
    

    // Test the result at root
    if (rank == ROOT_RANK) {
        CALI_MARK_BEGIN(correctness_check);
        // Sort locally as a reference
        qsort_ints(array, INCREASING);

        // cout << "Final sorted array from bitonic sort:" << endl;
        // print_array(sorted_array);

        // Compare the gathered result with the reference sorted array
        compare_arrays(sorted_array, array);
        CALI_MARK_END(correctness_check);
    }

    CALI_MARK_END(whole_computation);

    adiak::init(NULL);
    adiak::launchdate();                                          // launch date of the job
    adiak::libraries();                                           // Libraries used
    adiak::cmdline();                                             // Command line used to launch the job
    adiak::clustername();                                         // Name of the cluster
    adiak::value("algorithm", "column_sort");                         // The name of the algorithm you are using (e.g., "merge", "bitonic")
    adiak::value("programming_model", "mpi");         // e.g. "mpi"
    adiak::value("data_type", "int");                         // The datatype of input elements (e.g., double, int, float)
    adiak::value("size_of_data_type", sizeof(int));         // sizeof(datatype) of input elements in bytes (e.g., 1, 2, 4)
    adiak::value("input_size", N);                       // The number of elements in input dataset (1000)
    adiak::value("input_type", "Random");                       // For sorting, this would be choices: ("Sorted", "ReverseSorted", "Random", "1_perc_perturbed")
    adiak::value("num_procs", numprocs);                         // The number of processors (MPI ranks)
    adiak::value("scalability", "strong");                     // The scalability of your algorithm. choices: ("strong", "weak")
    adiak::value("group_num", 2);                      // The number of your group (integer, e.g., 1, 10)
    adiak::value("implementation_source", "online"); // Where you got the source code of your algorithm. choices: ("online", "ai", "handwritten").

	/******* Finalize *******/
    // Stop Caliper before finalizing MPI
    mgr.stop();
    mgr.flush(); 
	MPI_Finalize();
    
	return 0;
}
