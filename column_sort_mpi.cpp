#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

#include <algorithm>

#define MASTER 0
#define FROM_MASTER 1 /* setting a message type */
#define FROM_WORKER 2 /* setting a message type */
#define CHECK_SORTED 3

#define POS_INFINITY INT32_MAX
#define NEG_INFINITY INT32_MIN

int main(int argc, char *argv[]) {
    bool doPrint = false;
    if (argc != 4) {
        printf("usage: ./column_sort_mpi <number_of_processes> <number_of_matrix_rows> <number_of_matrix_columns\n");
        return 0;
    }
    // srand(time(NULL));
    int num_procs = atoi(argv[1]);
    int numRows = atoi(argv[2]);
    int numCols = atoi(argv[3]);

    if (numRows < 2 * (numCols - 1) * (numCols - 1)) {
        printf("<number_of_matrix_rows> must be >= 2 * (<number_of_matrix_columns> - 1)^2\n");
        return 0;
    }

    int matrixSize = numRows * numCols;

    int numColsPerWorker = numCols / num_procs;

    int rank;

    int dest, 
        src,
        offset;

    MPI_Status status;

    int matrix[numCols][numRows],
        shiftedMatrix[numCols+1][numRows],
        rowMatrix[numRows][numCols];
    
    int oneD[matrixSize];
    
    int colNum,
        colIdx,
        rowIdx;
    
    int sorted,
        isSorted;

    double transpose_time,
        untranspose_time,
        shift_time,
        unshift_time,
        sort_time,
        matrix_creation_time,
        check_time,
        total_time,
        total_sort_time;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // cali::ConfigManager mgr;
    // mgr.start();
    double total_time_start = MPI_Wtime();

    if (rank == MASTER) {
        // create matrix
        double matrix_create_start = MPI_Wtime();
        if (doPrint) { printf("\n"); }
        for (int i = 0; i < numCols; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i][j] = rand() % matrixSize;
                // matrix[i][j] = i*numRows + j;
                if (doPrint) { printf("%d,", matrix[i][j]); }
            }
            if (doPrint) { printf("\n"); }
        }
        if (doPrint) { printf("\n"); }
        double matrix_create_end = MPI_Wtime();
        matrix_creation_time = matrix_create_end - matrix_create_start;

        printf("Step 0: Matrix Created\n");
        printf("matrix_creation_time: %f \n\n", matrix_creation_time);

        double total_sort_start = MPI_Wtime();
        double sort_start = MPI_Wtime();
        
        offset = numColsPerWorker;
        for (dest = 1; dest < num_procs; ++dest) {
            MPI_Send(&matrix[offset], numColsPerWorker*numRows, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            // printf("sent %d cols to process %d\n", numColsPerWorker, dest);
            offset = offset + numColsPerWorker;
        }
        
        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            std::sort(matrix[colNum], matrix[colNum] + numRows);
        }

        memcpy(rowMatrix, matrix, numColsPerWorker*numRows*sizeof(int));

        offset = numColsPerWorker*numRows;
        for (src = 1; src < num_procs; ++src) {
            MPI_Recv((*rowMatrix+offset), numColsPerWorker*numRows, MPI_INT, src, FROM_WORKER, MPI_COMM_WORLD, &status);
            offset = offset + numColsPerWorker*numRows;
        }

        double sort_end = MPI_Wtime();
        sort_time = sort_end - sort_start;
        
        if (doPrint) {
            printf("\n");
            for (int i = 0; i < numRows; ++i) {
                for (int j = 0; j < numCols; ++j) {
                    printf("%d,", rowMatrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        printf("Step 1: Matrix Sorted\n");
        printf("sort_time: %f \n\n", sort_time);

        // transpose matrix start
        double transpose_time_start = MPI_Wtime();

        for (int i = 0; i < numCols; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i][j] = rowMatrix[j][i];
            }
        }

        double transpose_time_end = MPI_Wtime();
        transpose_time = transpose_time_end - transpose_time_start;
        // transpose matrix end

        if (doPrint) {
            printf("\n");
            for (int i = 0; i < numCols; ++i) {
                for (int j = 0; j < numRows; ++j) {
                    printf("%d,", matrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        printf("Step 2: Matrix Transposed\n");
        printf("transpose_time: %f \n\n", transpose_time);

        sort_start = MPI_Wtime();

        offset = numColsPerWorker;
        for (dest = 1; dest < num_procs; ++dest) {
            MPI_Send(&matrix[offset], numColsPerWorker * numRows, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            // printf("sent %d cols to process %d\n", numColsPerWorker, dest);
            offset = offset + numColsPerWorker;
        }

        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            std::sort(matrix[colNum], matrix[colNum] + numRows);
        }

        // memcpy(transposedMatrix, matrix, numColsPerWorker * numRows * sizeof(int));

        offset = numColsPerWorker;
        for (src = 1; src < num_procs; ++src) {
            MPI_Recv(&matrix[offset], numColsPerWorker * numRows, MPI_INT, src, FROM_WORKER, MPI_COMM_WORLD, &status);
            offset = offset + numColsPerWorker;
        }

        sort_end = MPI_Wtime();
        sort_time = sort_end - sort_start;
        
        if (doPrint) {
            printf("\n");
            for (int i = 0; i < numCols; ++i) {
                for (int j = 0; j < numRows; ++j) {
                    printf("%d,", matrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        printf("Step 3: Matrix Sorted\n");
        printf("sort_time: %f \n\n", sort_time);

        // untranspose matrix start
        double untranspose_start = MPI_Wtime();

        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                rowMatrix[i][j] = matrix[j][i];
            }
        }

        if (doPrint) {
            printf("\n");
            for (int i = 0; i < numRows; ++i) {
                for (int j = 0; j < numCols; ++j) {
                    printf("%d,", rowMatrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        memcpy(matrix, rowMatrix, matrixSize*sizeof(int));

        double untranspose_end = MPI_Wtime();
        untranspose_time = untranspose_end - untranspose_start;
        // untranspose matrix end

        if (doPrint) {
            printf("\n");
            for (int i = 0; i < numCols; ++i) {
                for (int j = 0; j < numRows; ++j) {
                    printf("%d,", matrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        printf("Step 4: Matrix Untransposed\n");
        printf("untranspose_time: %f \n\n", untranspose_time);

        sort_start = MPI_Wtime();
        offset = numColsPerWorker;
        for (dest = 1; dest < num_procs; ++dest) {
            MPI_Send(&matrix[offset], numColsPerWorker * numRows, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            // printf("sent %d cols to process %d\n", numColsPerWorker, dest);
            offset = offset + numColsPerWorker;
        }

        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            std::sort(matrix[colNum], matrix[colNum] + numRows);
        }

        offset = numColsPerWorker;
        for (src = 1; src < num_procs; ++src) {
            MPI_Recv(&matrix[offset], numColsPerWorker * numRows, MPI_INT, src, FROM_WORKER, MPI_COMM_WORLD, &status);
            offset = offset + numColsPerWorker;
        }

        sort_end = MPI_Wtime();
        sort_time = sort_end - sort_start;

        if (doPrint) {
            printf("\n");
            for (int i = 0; i < numCols; ++i) {
                for (int j = 0; j < numRows; ++j) {
                    printf("%d,", matrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        printf("Step 5: Matrix Sorted\n");
        printf("sort_time: %f \n\n", sort_time);

        // begin shift matrix
        double shift_start = MPI_Wtime();

        int beginShiftAmount = numRows / 2;
        int endShiftAmount = numRows / 2 + numRows % 2;
        int k = 0;
        for (int i = 0; i < numCols; ++i) {
            for (int j = 0; j < numRows; ++j) {
                oneD[k++] = matrix[i][j];
            }
        }
        for (int i = 0; i < beginShiftAmount; ++i) {
            shiftedMatrix[0][i] = NEG_INFINITY;
        }
        for (int i = 0; i < endShiftAmount; ++i) {
            shiftedMatrix[numCols][numRows - endShiftAmount + i] = POS_INFINITY;
        }
        k = 0;
        colIdx = 0;
        rowIdx = beginShiftAmount;
        while (k < matrixSize)
        {
            if (rowIdx % numRows == 0)
            {
                ++colIdx;
            }
            shiftedMatrix[colIdx][rowIdx % numRows] = oneD[k];
            ++k; ++rowIdx;
        }

        double shift_end = MPI_Wtime();
        shift_time = shift_end - shift_start;
        // end shift matrix

        if (doPrint) {
            printf("\n");
            for (int i = 0; i < numCols+1; ++i) {
                for (int j = 0; j < numRows; ++j) {
                    printf("%d,", shiftedMatrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        printf("Step 6: Matrix Shifted\n");
        printf("shift_time: %f \n\n", shift_time);

        sort_start = MPI_Wtime();

        offset = numColsPerWorker;
        for (dest = 1; dest < num_procs; ++dest) {
            MPI_Send(&shiftedMatrix[offset], numColsPerWorker * numRows, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            // printf("sent %d cols to process %d\n", numColsPerWorker, dest);
            offset = offset + numColsPerWorker;
        }

        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            std::sort(shiftedMatrix[colNum], shiftedMatrix[colNum] + numRows);
        }

        std::sort(shiftedMatrix[numCols], shiftedMatrix[numCols] + numRows);

        offset = numColsPerWorker;
        for (src = 1; src < num_procs; ++src) {
            MPI_Recv(&shiftedMatrix[offset], numColsPerWorker * numRows, MPI_INT, src, FROM_WORKER, MPI_COMM_WORLD, &status);
            offset = offset + numColsPerWorker;
        }

        sort_end = MPI_Wtime();
        sort_time = sort_end - sort_start;

        if (doPrint) {
            printf("\n");
            for (int i = 0; i < numCols + 1; ++i) {
                for (int j = 0; j < numRows; ++j) {
                    printf("%d,", shiftedMatrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        printf("Step 7: Matrix Sorted\n");
        printf("sort_time: %f \n\n", sort_time);
        
        // begin unshift matrix
        double unshift_start = MPI_Wtime();

        k = 0;
        for (int i = 0; i < numCols + 1; ++i) {
            for (int j = 0; j < numRows; ++j) {
                if (shiftedMatrix[i][j] == NEG_INFINITY || shiftedMatrix[i][j] == POS_INFINITY) {
                    continue;
                }
                oneD[k++] = shiftedMatrix[i][j];
            }
        }
        k = 0;
        colIdx = -1;
        rowIdx = 0;
        while (k < matrixSize) {
            if (rowIdx % numRows == 0) {
                ++colIdx;
            }
            matrix[colIdx][rowIdx % numRows] = oneD[k];
            ++k; ++rowIdx;
        }

        double unshift_end = MPI_Wtime();
        double total_sort_end = MPI_Wtime();
        unshift_time = unshift_end - unshift_start;
        total_sort_time = total_sort_end - total_sort_start;
        // end unshift matrix

        printf("Step 8: Matrix Unshifted\n");
        printf("unshift_time: %f \n\n", unshift_time);

        printf("total_sort_time: %f \n\n", total_sort_time);
    }

    if (rank != MASTER) {
        for (int i = 0; i < 4; ++i) {
            MPI_Recv(&matrix, numColsPerWorker*numRows, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);

            for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
                std::sort(matrix[colNum], matrix[colNum] + numRows);
            }

            MPI_Send(&matrix, numColsPerWorker*numRows, MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
        }
    }

    double check_start = MPI_Wtime();

    if (rank == MASTER) {
        printf("Step 9: Check if Matrix is Sorted\n");
        if (doPrint) {
            printf("\n");
            for (int i = 0; i < numCols; ++i) {
                for (int j = 0; j < numRows; ++j) {
                    printf("%d,", matrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }
        offset = numColsPerWorker;
        for (dest = 1; dest < num_procs; ++dest) {
            MPI_Send(&matrix[offset], numColsPerWorker * numRows, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            // printf("sent %d cols to process %d\n", numColsPerWorker, dest);
            offset = offset + numColsPerWorker;
        }
    }
    if (rank % 2 == 1) {
        MPI_Recv(&matrix, numColsPerWorker * numRows, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
        sorted = 1;
        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            sorted = sorted && std::is_sorted(matrix[colNum], matrix[colNum] + numRows);
        }
        MPI_Send(&matrix[0][0], 1, MPI_INT, rank-1, CHECK_SORTED, MPI_COMM_WORLD);
        MPI_Recv(&isSorted, 1, MPI_INT, rank-1, CHECK_SORTED, MPI_COMM_WORLD, &status);
        sorted = sorted && isSorted;
        MPI_Send(&sorted, 1, MPI_INT, MASTER, CHECK_SORTED, MPI_COMM_WORLD);
    }
    else if (rank % 2 == 0) {
        if (rank != MASTER) {
            MPI_Recv(&matrix, numColsPerWorker * numRows, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
        }
        isSorted = 1;
        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            isSorted = isSorted && std::is_sorted(matrix[colNum], matrix[colNum] + numRows);
        }
        if (rank+1 < num_procs) {
            MPI_Recv(&sorted, 1, MPI_INT, rank+1, CHECK_SORTED, MPI_COMM_WORLD, &status);
            if (matrix[numColsPerWorker-1][numRows-1] > sorted) {
                sorted = 0;
            } else {
                sorted = 1;
            }
            isSorted = isSorted && sorted;
            MPI_Send(&isSorted, 1, MPI_INT, rank+1, CHECK_SORTED, MPI_COMM_WORLD);
        }
        else {
            MPI_Send(&isSorted, 1, MPI_INT, MASTER, CHECK_SORTED, MPI_COMM_WORLD);
        }
    }

    if (rank == MASTER) {
        sorted = 1;
        for (src = 1; src < num_procs; src += 2) {
            MPI_Recv(&isSorted, 1, MPI_INT, src, CHECK_SORTED, MPI_COMM_WORLD, &status);
            sorted = sorted && isSorted;
        }
        if (num_procs % 2 == 1) {
            MPI_Recv(&isSorted, 1, MPI_INT, num_procs-1, CHECK_SORTED, MPI_COMM_WORLD, &status);
            sorted = sorted && isSorted;
        }
        if (sorted) {
            printf("\nmatrix is sorted\n\n");
        }
        else {
            printf("\nmatrix not sorted\n\n");
        }
        double check_end = MPI_Wtime();
        double total_time_end = MPI_Wtime();
        check_time = check_end - check_start;
        total_time = total_time_end - total_time_start;
        printf("check_time: %f \n\n", check_time);
        printf("total_time: %f \n\n", total_time);
    }
    // mgr.stop();
    // mgr.flush();

    MPI_Finalize();
}