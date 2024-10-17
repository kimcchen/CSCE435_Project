
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>
#include <vector>
#include <algorithm>

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
  
	/******* Initialize *******/
	rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Get the matrix size from argv[1], and ensure argv[1] is passed.
    if (argc < 2) {
        if (rank == ROOT_RANK) {
            cout << "Matrix size not provided. Terminating." << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    int N = atoi(argv[1]); // Convert the first argument to an integer, which is the matrix size
    cout << "Matrix size: " << N << endl;
  
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
		cout << "The number of elements in the array is" << length << endl;
        array = create_random_int_array(length);

        // Print the unsorted global array
        cout << "The unsorted global array (before scattering):" << endl;
        print_array(array);
	}
  
  /***** Do that sorting thing *******/
  
    MPI_Scatter(array.data(), length_per_process, MPI_INT, my_array.data(), length_per_process, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    // Sort locally
    // qsort_ints(my_array, INCREASING);
    qsort_ints(my_array, bit(rank, 0) == 0 ? INCREASING : DECREASING);
    // Perform the bitonic merge
    for (int depth = 1; depth <= log2(numprocs); depth++) {
        for (int p = depth - 1; p >= 0; p--) {
            int prank = flipbit(rank, p);
            parray.resize(length_per_process);

            // Exchange data with partner
            MPI_Sendrecv(my_array.data(), length_per_process, MPI_INT, prank, TAG,
                         parray.data(), length_per_process, MPI_INT, prank, TAG,
                         MPI_COMM_WORLD, &status);

            if (bit(rank, depth) == 0) {
                merge_up(my_array, rank, parray, prank);
            } else {
                merge_down(my_array, rank, parray, prank);
            }
        }
        // Local sort after each depth, based on bit value of rank at depth
        qsort_ints(my_array, bit(rank, depth) == 0 ? INCREASING : DECREASING);
    }
     // Gather the results
    if (rank == ROOT_RANK) {
        // Sort locally as a reference
        // qsort_ints(array, INCREASING);

        sorted_array.resize(length);  // Allocate space for sorted array

        // Compare the gathered result with the reference sorted array
        // compare_arrays(sorted_array, array);
    }
    // Gather the results
    MPI_Gather(my_array.data(), length_per_process, MPI_INT, sorted_array.data(), length_per_process, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

    // Test the result at root
    if (rank == ROOT_RANK) {
        // Sort locally as a reference
        qsort_ints(array, INCREASING);

        cout << "Final sorted array from bitonic sort:" << endl;
        print_array(sorted_array);

        // Compare the gathered result with the reference sorted array
        compare_arrays(sorted_array, array);
    }

	/******* Finalize *******/
	MPI_Finalize();
    
	return 0;
}
