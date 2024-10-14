#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

#include <algorithm>
#include <math.h>

#define MASTER 0
#define FROM_MASTER 1 /* setting a message type */
#define FROM_WORKER 2 /* setting a message type */
#define CHECK_SORTED 3

#define POS_INFINITY INT32_MAX
#define NEG_INFINITY INT32_MIN

int main(int argc, char *argv[]) {
    // bool doPrint = false;
    if (argc != 4) {
        printf("usage: ./column_sort_mpi <number_of_matrix_rows> <number_of_matrix_columns <array_type>\n");
        return 0;
    }
    srand(time(NULL));
    int numRows = atoi(argv[1]);
    int numCols = atoi(argv[2]);
    int arrayType = atoi(argv[3]);

    if (numRows < 2 * (numCols - 1) * (numCols - 1) || numRows <= 0 || numCols <= 0) {
        printf("<number_of_matrix_rows> must be >= 2 * (<number_of_matrix_columns> - 1)^2\n");
        printf("<number_of_matrix_rows> must be > 0\n");
        printf("<number_of_matrix_columns> must be > 0\n");
        return 0;
    }
    if (arrayType < 0 || arrayType > 3) {
        printf("<array_type> must be 0, 1, 2, 3\n");
        printf("0 is sorted\n1 is random\n2 is reverse sorted\n 3 is 1\% perturbed\n");
        return 0;
    }

    int matrixSize = numRows * numCols;

    int rank,
        num_procs;

    int dest, 
        src,
        offset;

    MPI_Status status;

    int matrix[numCols][numRows],
        shiftedMatrix[numCols+1][numRows],
        rowMatrix[numRows][numCols];
    
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

    /* Define Caliper region names */
    const char *transpose = "transpose";
    const char *untranspose = "untranspose";
    const char *shift = "shift";
    const char *unshift = "unshift";
    const char *sort1 = "sort1";
    const char *sort2 = "sort2";
    const char *sort3 = "sort3";
    const char *sort4 = "sort4";
    const char *whole_sort = "whole_sort";
    const char *matrix_creation = "matrix_creation";
    const char *whole_computation = "whole_computation";
    const char *sort_check = "sort_check";

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int numColsPerWorker = numCols / num_procs;

    cali::ConfigManager mgr;
    mgr.start();

    double total_time_start = MPI_Wtime();
    CALI_MARK_BEGIN(whole_computation);

    // if (rank == MASTER) {
    // create matrix
    double matrix_create_start = MPI_Wtime();
    CALI_MARK_BEGIN(matrix_creation);

    if (arrayType == 0) {
        int start = rank * numColsPerWorker * numRows;
        for (int i = 0; i < numColsPerWorker; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i][j] = start + (i*numColsPerWorker) + j;
            }
        }
    }
    else if (arrayType == 1) {
        for (int i = 0; i < numColsPerWorker; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i][j] = rand() % matrixSize;
                // matrix[i][j] = i*numRows + j;
            }
        }
    }
    else if (arrayType == 2) {
        int start = matrixSize - (rank * numColsPerWorker * numRows);
        for (int i = 0; i < numColsPerWorker; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i][j] = start - (i*numColsPerWorker) - j;
            }
        }
    }
    else if (arrayType == 3) {
        int start = rank * numColsPerWorker * numRows;
        for (int i = 0; i < numColsPerWorker; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i][j] = start + (i * numColsPerWorker) + j;
            }
        }

        int numPertubed = std::ceil(numColsPerWorker*numRows*0.01);
        for (int i = 0; i < numPertubed; ++i) {
            int colIdx1 = rand() % numColsPerWorker;
            int rowIdx1 = rand() % numRows;
            int colIdx2 = rand() % numColsPerWorker;
            if (numColsPerWorker != 1) {
                while (colIdx1 != colIdx2) { colIdx2 = rand() % numRows; }
            }
            int rowIdx2 = rand() % numRows;
            while (rowIdx1 == rowIdx2) { rowIdx2 = rand() % numRows; }

            std::swap(matrix[colIdx1][rowIdx1], matrix[colIdx2][rowIdx2]);
        }
    }

    // if (doPrint) {
    //     printf("\n");
    //     for (int i = 0; i < numCols; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    CALI_MARK_END(matrix_creation);
    double matrix_create_end = MPI_Wtime();
    matrix_creation_time = matrix_create_end - matrix_create_start;

    printf("Step 0: Matrix Created\n");
    printf("matrix_creation_time: %f \n\n", matrix_creation_time);

    double total_sort_start = MPI_Wtime();
    double sort_start = MPI_Wtime();
    if (rank == MASTER) {
        CALI_MARK_BEGIN(whole_sort);
        CALI_MARK_BEGIN(sort1);
    }
    
    // offset = numColsPerWorker;
    // for (dest = 1; dest < num_procs; ++dest) {
    //     MPI_Send(&matrix[offset], numColsPerWorker*numRows, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
    //     // printf("sent %d cols to process %d\n", numColsPerWorker, dest);
    //     offset = offset + numColsPerWorker;
    // }
    
    for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
        std::sort(matrix[colNum], matrix[colNum] + numRows);
    }

    if (rank == MASTER) {
        memcpy(rowMatrix, matrix, numColsPerWorker * numRows * sizeof(int));
        offset = numColsPerWorker*numRows;
        for (src = 1; src < num_procs; ++src) {
            MPI_Recv((*rowMatrix+offset), numColsPerWorker*numRows, MPI_INT, src, FROM_WORKER, MPI_COMM_WORLD, &status);
            offset = offset + numColsPerWorker*numRows;
        }

        CALI_MARK_END(sort1);
        double sort_end = MPI_Wtime();
        sort_time = sort_end - sort_start;
        
        // if (doPrint) {
        //     printf("\n");
        //     for (int i = 0; i < numRows; ++i) {
        //         for (int j = 0; j < numCols; ++j) {
        //             printf("%d,", rowMatrix[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

        printf("Step 1: Matrix Sorted\n");
        printf("sort_time: %f \n\n", sort_time);

        // transpose matrix start
        double transpose_time_start = MPI_Wtime();
        CALI_MARK_BEGIN(transpose);

        for (int i = 0; i < numCols; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i][j] = rowMatrix[j][i];
            }
        }

        CALI_MARK_END(transpose);
        double transpose_time_end = MPI_Wtime();
        transpose_time = transpose_time_end - transpose_time_start;
        // transpose matrix end

        // if (doPrint) {
        //     printf("\n");
        //     for (int i = 0; i < numCols; ++i) {
        //         for (int j = 0; j < numRows; ++j) {
        //             printf("%d,", matrix[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

        printf("Step 2: Matrix Transposed\n");
        printf("transpose_time: %f \n\n", transpose_time);

        sort_start = MPI_Wtime();
        CALI_MARK_BEGIN(sort2);

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

        CALI_MARK_END(sort2);
        sort_end = MPI_Wtime();
        sort_time = sort_end - sort_start;
        
        // if (doPrint) {
        //     printf("\n");
        //     for (int i = 0; i < numCols; ++i) {
        //         for (int j = 0; j < numRows; ++j) {
        //             printf("%d,", matrix[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

        printf("Step 3: Matrix Sorted\n");
        printf("sort_time: %f \n\n", sort_time);

        // untranspose matrix start
        double untranspose_start = MPI_Wtime();
        CALI_MARK_BEGIN(untranspose);

        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                rowMatrix[i][j] = matrix[j][i];
            }
        }
        
        CALI_MARK_END(untranspose);
        double untranspose_end = MPI_Wtime();
        untranspose_time = untranspose_end - untranspose_start;
        // untranspose matrix end

        // if (doPrint) {
        //     printf("\n");
        //     for (int i = 0; i < numRows; ++i) {
        //         for (int j = 0; j < numCols; ++j) {
        //             printf("%d,", rowMatrix[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

        printf("Step 4: Matrix Untransposed\n");
        printf("untranspose_time: %f \n\n", untranspose_time);

        sort_start = MPI_Wtime();
        CALI_MARK_BEGIN(sort3);

        offset = numColsPerWorker*numRows;
        for (dest = 1; dest < num_procs; ++dest) {
            MPI_Send((*rowMatrix+offset), numColsPerWorker * numRows, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            // printf("sent %d cols to process %d\n", numColsPerWorker, dest);
            offset = offset + numColsPerWorker*numRows;
        }

        memcpy(matrix, rowMatrix, numColsPerWorker*numRows*sizeof(int));

        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            std::sort(matrix[colNum], matrix[colNum] + numRows);
        }

        offset = numColsPerWorker;
        for (src = 1; src < num_procs; ++src) {
            MPI_Recv(&matrix[offset], numColsPerWorker * numRows, MPI_INT, src, FROM_WORKER, MPI_COMM_WORLD, &status);
            offset = offset + numColsPerWorker;
        }

        CALI_MARK_END(sort3);
        sort_end = MPI_Wtime();
        sort_time = sort_end - sort_start;

        // if (doPrint) {
        //     printf("\n");
        //     for (int i = 0; i < numCols; ++i) {
        //         for (int j = 0; j < numRows; ++j) {
        //             printf("%d,", matrix[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

        printf("Step 5: Matrix Sorted\n");
        printf("sort_time: %f \n\n", sort_time);

        // begin shift matrix
        double shift_start = MPI_Wtime();
        CALI_MARK_BEGIN(shift);

        int beginShiftAmount = numRows / 2;
        int endShiftAmount = numRows / 2 + numRows % 2;

        for (int i = 0; i < beginShiftAmount; ++i) {
            shiftedMatrix[0][i] = NEG_INFINITY;
        }
        for (int i = 0; i < endShiftAmount; ++i) {
            shiftedMatrix[numCols][numRows - endShiftAmount + i] = POS_INFINITY;
        }
        memcpy((*shiftedMatrix+beginShiftAmount+1), matrix, matrixSize*sizeof(int));

        CALI_MARK_END(shift);
        double shift_end = MPI_Wtime();
        shift_time = shift_end - shift_start;
        // end shift matrix

        // if (doPrint) {
        //     printf("\n");
        //     for (int i = 0; i < numCols+1; ++i) {
        //         for (int j = 0; j < numRows; ++j) {
        //             printf("%d,", shiftedMatrix[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

        printf("Step 6: Matrix Shifted\n");
        printf("shift_time: %f \n\n", shift_time);

        sort_start = MPI_Wtime();
        CALI_MARK_BEGIN(sort4);

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

        CALI_MARK_END(sort4);
        sort_end = MPI_Wtime();
        sort_time = sort_end - sort_start;

        // if (doPrint) {
        //     printf("\n");
        //     for (int i = 0; i < numCols + 1; ++i) {
        //         for (int j = 0; j < numRows; ++j) {
        //             printf("%d,", shiftedMatrix[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

        printf("Step 7: Matrix Sorted\n");
        printf("sort_time: %f \n\n", sort_time);
        
        // begin unshift matrix
        double unshift_start = MPI_Wtime();
        CALI_MARK_BEGIN(unshift);

        memcpy(matrix, (*shiftedMatrix+beginShiftAmount+1), matrixSize*sizeof(int));

        CALI_MARK_END(unshift);
        CALI_MARK_END(whole_sort);
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
        MPI_Send(&matrix, numColsPerWorker * numRows, MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
        for (int i = 0; i < 3; ++i) {
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

        CALI_MARK_BEGIN(sort_check);

        // if (doPrint) {
        //     printf("\n");
        //     for (int i = 0; i < numCols; ++i) {
        //         for (int j = 0; j < numRows; ++j) {
        //             printf("%d,", matrix[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

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
        // printf("Step 9: Check if Matrix is Sorted\n");

        // double check_start = MPI_Wtime();
        // CALI_MARK_BEGIN(sort_check);
        // if (doPrint) {
        //     printf("\n");
        //     for (int i = 0; i < numCols; ++i) {
        //         for (int j = 0; j < numRows; ++j) {
        //             printf("%d,", matrix[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

        // offset = numColsPerWorker;
        // for (dest = 1; dest < num_procs; ++dest) {
        //     MPI_Send(&matrix[offset], numColsPerWorker * numRows, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
        //     // printf("sent %d cols to process %d\n", numColsPerWorker, dest);
        //     offset = offset + numColsPerWorker;
        // }

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
        CALI_MARK_END(sort_check);
        double check_end = MPI_Wtime();
        double total_time_end = MPI_Wtime();
        check_time = check_end - check_start;
        total_time = total_time_end - total_time_start;
        printf("check_time: %f \n\n", check_time);
        printf("total_time: %f \n\n", total_time);
    }

    CALI_MARK_END(whole_computation);

    adiak::init(NULL);
    adiak::user();
    adiak::launchdate();
    adiak::libraries();
    adiak::cmdline();
    adiak::clustername();
    adiak::value("num_procs", num_procs);
    adiak::value("matrix_size", matrixSize);
    adiak::value("program_name", "column_sort");
    adiak::value("matrix_datatype_size", sizeof(int));

    mgr.stop();
    mgr.flush();

    MPI_Finalize();
}