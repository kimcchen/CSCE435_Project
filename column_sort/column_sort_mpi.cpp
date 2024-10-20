#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

#include <algorithm>
#include <math.h>
#include <random>

#define POS_INFINITY INT32_MAX
#define NEG_INFINITY INT32_MIN

#define TRANSPOSE 0
#define UNTRANSPOSE 1
#define SHIFT 2
#define UNSHIFT 3
#define CHECK_SORTED 4

#define MASTER 0

using std::is_sorted;
using std::rand;
using std::sort;
using std::srand;

int main(int argc, char *argv[]) {
    CALI_CXX_MARK_FUNCTION;
    if (argc != 4) {
        printf("usage: ./column_sort_mpi <number_of_matrix_rows> <number_of_matrix_columns> <array_type>\n");
        return 0;
    }
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

    int rank, num_procs;
    int colNum;
    MPI_Status status;
    int sorted, isSorted;

    int amountToOffset,
        amountToSend,
        offset,
        recvOffset;

    int dest, src;

    double matrix_creation_time,
        matrix_transpose_time,
        matrix_untranspose_time,
        matrix_sort_time,
        matrix_check_time,
        matrix_shift_time,
        matrix_unshift_time,
        matrix_intermediate_sort_time,
        total_time;

    double matrix_creation_time_start, matrix_creation_time_end,
        matrix_transpose_time_start, matrix_tranpose_time_end,
        matrix_untranspose_time_start, matrix_untranspose_time_end,
        matrix_sort_time_start, matrix_sort_time_end,
        matrix_check_time_start, matrix_check_time_end,
        matrix_shift_time_start, matrix_shift_time_end,
        matrix_unshift_time_start, matrix_unshift_time_end,
        matrix_intermediate_sort_time_start, matrix_intermediate_sort_time_end,
        total_time_start, total_time_end;

    /* Define Caliper region names */
    // const char *transpose = "transpose";
    // const char *untranspose = "untranspose";
    // const char *shift = "shift";
    // const char *unshift = "unshift";
    // const char *sort1 = "sort1";
    // const char *sort2 = "sort2";
    // const char *sort3 = "sort3";
    // const char *sort4 = "sort4";
    // const char *whole_sort = "whole_sort";
    // const char *matrix_creation = "matrix_creation";
    // const char *whole_computation = "whole_computation";
    // const char *sort_check = "sort_check";
    // const char *input_type;
    const char *data_init_X = "data_init_X";
    const char *comm = "comm";
    const char *comm_small = "comm_small";
    const char *comm_large = "comm_large";
    const char *comp = "comp";
    const char *comp_small = "comp_small";
    const char *comp_large = "comp_large";
    const char *correctness_check = "correctness_check";
    const char *input_type;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int numColsPerWorker = numCols / num_procs;
    int localMatrixSize = numColsPerWorker * numRows;

    int *matrix = new int[localMatrixSize];
    int *rowMatrix = new int[localMatrixSize];
    int *finalCol = new int[numRows];

    cali::ConfigManager mgr;
    mgr.start();

    if (rank == MASTER) {
        matrix_creation_time_start = MPI_Wtime();
        total_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(whole_computation);
        // CALI_MARK_BEGIN(matrix_creation);
    }

    CALI_MARK_BEGIN(data_init_X);

    if (arrayType == 0) {
        int start = rank * numColsPerWorker * numRows;
        for (int i = 0; i < numColsPerWorker; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i*numRows + j] = start++;
            }
        }
        input_type = "Sorted";
    }
    else if (arrayType == 1) {
        srand(time(NULL));
        for (int i = 0; i < numColsPerWorker; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i*numColsPerWorker + j] = rand() % localMatrixSize;
            }
        }
        input_type = "Random";
    }
    else if (arrayType == 2) {
        int start = localMatrixSize - (rank * numColsPerWorker * numRows);
        for (int i = 0; i < numColsPerWorker; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i*numColsPerWorker + j] = start - (i*numColsPerWorker) - j;
            }
        }
        input_type = "ReverseSorted";
    }
    else if (arrayType == 3) {
        int start = rank * numColsPerWorker * numRows;
        for (int i = 0; i < numColsPerWorker; ++i) {
            for (int j = 0; j < numRows; ++j) {
                matrix[i*numColsPerWorker + j] = start + (i * numColsPerWorker) + j;
            }
        }

        int numPertubed = std::ceil(numColsPerWorker*numRows*0.01)/2;
        for (int i = 0; i < numPertubed; ++i) {
            int colIdx1 = rand() % numColsPerWorker;
            int rowIdx1 = rand() % numRows;
            int colIdx2 = rand() % numColsPerWorker;
            if (numColsPerWorker != 1) {
                while (colIdx1 != colIdx2) { colIdx2 = rand() % numRows; }
            }
            int rowIdx2 = rand() % numRows;
            while (rowIdx1 == rowIdx2) { rowIdx2 = rand() % numRows; }

            std::swap(matrix[colIdx1*numRows + rowIdx1], matrix[colIdx2*numRows + rowIdx2]);
        }

        input_type = "1_perc_perturbed";
    }

    CALI_MARK_END(data_init_X);

    if (rank == MASTER) {
        // CALI_MARK_END(matrix_creation);
        matrix_creation_time_end = MPI_Wtime();
        matrix_creation_time = matrix_check_time_end - matrix_check_time_start;
        printf("Step 0: Created Matrix\n");
        printf("time: %f\n\n", matrix_creation_time);
    }

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    if (rank == MASTER) {
        matrix_sort_time_start = MPI_Wtime();
        matrix_intermediate_sort_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(whole_sort);
        // CALI_MARK_BEGIN(sort1);
    }

    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_large);

    for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
        sort(&matrix[colNum * numRows], &matrix[colNum * numRows + numRows]);
    }

    CALI_MARK_END(comp_large);
    CALI_MARK_END(comp);

    if (rank == MASTER) {
        // CALI_MARK_END(sort1);
        matrix_intermediate_sort_time_end = MPI_Wtime();
        matrix_intermediate_sort_time = matrix_intermediate_sort_time_end - matrix_intermediate_sort_time_start;
        printf("Step 1: Sorted Columns of Matrix\n");
        printf("time: %f\n\n", matrix_intermediate_sort_time);
    }

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    if (rank == MASTER) {
        matrix_transpose_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(transpose);
    }

    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_small);
    // create a numRows/num_procs x numCols matrix in rowMatrix
    offset = 0;
    amountToOffset = numRows / num_procs;
    while (offset < localMatrixSize) {
        memcpy(&rowMatrix[offset], &matrix[offset], amountToOffset * sizeof(int));
        offset = offset + amountToOffset;
    }

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < amountToOffset; ++i) {
    //         for (int j = 0; j < numCols; ++j) {
    //             printf("%d,", rowMatrix[i*numCols + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    // matrix is now a numCols x numRows/num_procs
    for (int i = 0; i < amountToOffset; ++i) {
        for (int j = 0; j < numCols; ++j) {
            matrix[j*amountToOffset + i] = rowMatrix[i*numCols + j];
        }
    }

    CALI_MARK_END(comp_small);
    CALI_MARK_END(comp);

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numCols; ++i) {
    //         for (int j = 0; j < amountToOffset; ++j) {
    //             printf("%d,", matrix[i*amountToOffset + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    CALI_MARK_BEGIN(comm);
    CALI_MARK_BEGIN(comm_large);

    amountToSend = amountToOffset * numCols / num_procs;
    dest = 0;
    offset = 0;
    int currIdx = 0;
    int startRank = numColsPerWorker * rank;
    recvOffset = startRank;
    colNum = 0;
    offset = 0;

    while (offset < localMatrixSize) {
        if (dest == num_procs) {
            dest = 0;
        }
        // for (colNum )
        MPI_Send(&matrix[offset], amountToSend, MPI_INT, dest, TRANSPOSE, MPI_COMM_WORLD);
        ++dest;
        offset = offset + amountToSend;
    }

    src = 0;
    offset = 0;
    while (offset < localMatrixSize) {
        if (src == num_procs) {
            src = 0;
        }
        MPI_Recv(&matrix[offset], amountToSend, MPI_INT, src, TRANSPOSE, MPI_COMM_WORLD, &status);
        ++src;
        offset = offset + amountToSend;
    }

    CALI_MARK_END(comm_large);
    CALI_MARK_END(comm);

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_small);

    offset = 0;
    colNum = 0;
    recvOffset = 0;
    int count = 0;
    while (recvOffset < localMatrixSize) {
        if (offset >= localMatrixSize) {
            ++count;
            offset = amountToOffset * count;
        }
        memcpy(&rowMatrix[recvOffset], &matrix[offset], amountToOffset * sizeof(int));
        recvOffset = recvOffset + amountToOffset;
        offset = offset + amountToSend;
    }

    memcpy(&matrix[0], &rowMatrix[0], localMatrixSize * sizeof(int));

    CALI_MARK_END(comp_small);
    CALI_MARK_END(comp);

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    if (rank == MASTER) {
        // CALI_MARK_END(transpose);
        matrix_tranpose_time_end = MPI_Wtime();
        matrix_transpose_time = matrix_tranpose_time_end - matrix_transpose_time_start;
        printf("Step 2: Transposed Matrix\n");
        printf("time: %f\n\n", matrix_transpose_time);
    }

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    if (rank == MASTER) {
        matrix_intermediate_sort_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(sort2);
    }

    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_large);

    for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
        sort(&matrix[colNum * numRows], &matrix[colNum * numRows + numRows]);
    }

    CALI_MARK_END(comp_large);
    CALI_MARK_END(comp);

    if (rank == MASTER) {
        // CALI_MARK_END(sort2);
        matrix_intermediate_sort_time_end = MPI_Wtime();
        matrix_intermediate_sort_time = matrix_intermediate_sort_time_end - matrix_intermediate_sort_time_start;
        printf("Step 3: Sorted Columns of Matrix\n");
        printf("time: %f\n\n", matrix_intermediate_sort_time);
    }

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    if (rank == MASTER) {
        matrix_untranspose_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(untranspose);
    }

    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_small);

    // untranspose matrix
    for (int i = 0; i < numColsPerWorker; ++i) {
        for (int j = 0; j < numRows; ++j) {
            rowMatrix[j*numColsPerWorker + i] = matrix[i*numRows + j];
        }
    }

    CALI_MARK_END(comp_small);
    CALI_MARK_END(comp);

    CALI_MARK_BEGIN(comm);
    CALI_MARK_BEGIN(comm_large);

    offset = 0;
    dest = 0;
    currIdx = 0;
    amountToSend = numRows / num_procs * numColsPerWorker;
    // if (rank == MASTER) { printf("amountToSend: %d\n", amountToSend); }

    for (dest = 0; offset < localMatrixSize; ++dest) {
        if (dest == num_procs) {
            dest = 0;
        }
        MPI_Send((&rowMatrix[offset]), amountToSend, MPI_INT, dest, UNTRANSPOSE, MPI_COMM_WORLD);
        offset = offset + amountToSend;
    }

    colNum = 0;
    src = 0;
    offset = 0;
    for (src = 0; offset < localMatrixSize; ++src) {
        if (src == num_procs) {
            src = 0;
        }
        MPI_Recv((&matrix[offset]), amountToSend, MPI_INT, src, UNTRANSPOSE, MPI_COMM_WORLD, &status);
        offset = offset + amountToSend;
    }

    CALI_MARK_END(comm_large);
    CALI_MARK_END(comm);

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_small);

    amountToOffset = numRows / numCols * numColsPerWorker;
    offset = 0;
    recvOffset = 0;
    while (recvOffset < numRows) {
        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            memcpy(&rowMatrix[colNum*numRows + recvOffset], &matrix[offset], amountToOffset * sizeof(int));
            offset  = offset + amountToOffset;
        }
        recvOffset = recvOffset + amountToOffset;
    }

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numRows; ++i) {
    //         for (int j = 0; j < numColsPerWorker; ++j) {
    //             printf("%d,", rowMatrix[i*numColsPerWorker + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    offset = 0;
    while (offset < localMatrixSize) {
        memcpy((&matrix[offset]), (&rowMatrix[offset]), numRows * sizeof(int));
        offset = offset + numRows;
    }

    CALI_MARK_END(comp_small);
    CALI_MARK_END(comp);

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    if (rank == MASTER) {
        // CALI_MARK_END(untranspose);
        matrix_untranspose_time_end = MPI_Wtime();
        matrix_untranspose_time = matrix_untranspose_time_end - matrix_untranspose_time_start;
        printf("Step 4: Untransposed Matrix\n");
        printf("time: %f\n\n", matrix_untranspose_time);
    }

    if (rank == MASTER) {
        matrix_intermediate_sort_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(sort3);
    }

    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_large);

    for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
        sort(&matrix[colNum * numRows], &matrix[colNum * numRows + numRows]);
    }

    CALI_MARK_END(comp_large);
    CALI_MARK_END(comp);

    if (rank == MASTER) {
        // CALI_MARK_END(sort3);
        matrix_intermediate_sort_time_end = MPI_Wtime();
        matrix_intermediate_sort_time = matrix_intermediate_sort_time_end - matrix_intermediate_sort_time_start;
        printf("Step 5: Sorted Columns of Matrix\n");
        printf("time: %f\n\n", matrix_intermediate_sort_time);
    }

    if (rank == MASTER) {
        matrix_shift_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(shift);
    }

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    int beginShiftAmount = numRows / 2;
    int amountToShift = numRows - beginShiftAmount;

    if (rank != num_procs-1) {
        CALI_MARK_BEGIN(comm);
        CALI_MARK_BEGIN(comm_small);

        MPI_Send(&matrix[localMatrixSize-beginShiftAmount], beginShiftAmount, MPI_INT, rank+1, SHIFT, MPI_COMM_WORLD);

        CALI_MARK_END(comm_small);
        CALI_MARK_END(comm);

        CALI_MARK_BEGIN(comp);
        CALI_MARK_BEGIN(comp_small);

        // memcpy(&matrix[beginShiftAmount], &matrix[0], (localMatrixSize-beginShiftAmount) * sizeof(int));
        for (colNum = numColsPerWorker - 1; colNum > 0; --colNum) {
            memcpy(&matrix[colNum * numRows + beginShiftAmount], &matrix[colNum * numRows], (amountToShift) * sizeof(int));
            memcpy(&matrix[colNum * numRows], &matrix[(colNum - 1) * numRows + amountToShift], beginShiftAmount * sizeof(int));
        }
        memcpy(&matrix[beginShiftAmount], &matrix[0], amountToShift * sizeof(int));

        CALI_MARK_END(comp_small);
        CALI_MARK_END(comp);
    }
    else { 
        CALI_MARK_BEGIN(comp);
        CALI_MARK_BEGIN(comp_small);

        memcpy(finalCol, &matrix[localMatrixSize-numRows+amountToShift], beginShiftAmount*sizeof(int));
        for (int i = beginShiftAmount; i < numRows; ++i) {
            finalCol[i] = POS_INFINITY;
        }
        for (colNum = numColsPerWorker - 1; colNum > 0; --colNum) {
            memcpy(&matrix[colNum * numRows + beginShiftAmount], &matrix[colNum * numRows], (amountToShift) * sizeof(int));
            memcpy(&matrix[colNum*numRows], &matrix[(colNum-1)*numRows + amountToShift], beginShiftAmount * sizeof(int));
        }
        memcpy(&matrix[beginShiftAmount], &matrix[0], amountToShift * sizeof(int));

        CALI_MARK_END(comp_small);
        CALI_MARK_END(comp);
    }

    if (rank != 0) {
        CALI_MARK_BEGIN(comm);
        CALI_MARK_BEGIN(comm_small);

        MPI_Recv(&matrix[0], beginShiftAmount, MPI_INT, rank - 1, SHIFT, MPI_COMM_WORLD, &status);

        CALI_MARK_END(comm_small);
        CALI_MARK_END(comm);
    }
    else {
        CALI_MARK_BEGIN(comp);
        CALI_MARK_BEGIN(comp_small);

        for (int i = 0; i < beginShiftAmount; ++i) {
            matrix[i] = NEG_INFINITY;
        }

        CALI_MARK_END(comp_small);
        CALI_MARK_END(comp);
    }

    if (rank == MASTER) {
        // CALI_MARK_END(shift);
        matrix_shift_time_end = MPI_Wtime();
        matrix_shift_time = matrix_shift_time_end - matrix_shift_time_start;
        printf("\nStep 6: Shifted Matrix\n");
        printf("time: %f\n\n", matrix_shift_time);
    }

    if (rank == MASTER) {
        matrix_intermediate_sort_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(sort4);
    }

    CALI_MARK_BEGIN(comp);
    CALI_MARK_BEGIN(comp_large);

    for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
        sort(&matrix[colNum * numRows], &matrix[colNum * numRows + numRows]);
    }

    CALI_MARK_END(comp_large);
    CALI_MARK_END(comp);

    if (rank == MASTER) {
        // CALI_MARK_END(sort4);
        matrix_intermediate_sort_time_end = MPI_Wtime();
        matrix_intermediate_sort_time = matrix_intermediate_sort_time_end - matrix_intermediate_sort_time_start;
        printf("Step 7: Sorted Columns of Matrix\n");
        printf("time: %f\n\n", matrix_intermediate_sort_time);
    }

    if (rank == MASTER) {
        matrix_unshift_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(unshift);
    }

    if (rank != 0) {
        CALI_MARK_BEGIN(comm);
        CALI_MARK_BEGIN(comm_small);

        MPI_Send(&matrix[0], beginShiftAmount, MPI_INT, rank-1, UNSHIFT, MPI_COMM_WORLD);

        CALI_MARK_END(comm_small);
        CALI_MARK_END(comm);

        CALI_MARK_BEGIN(comp);
        CALI_MARK_BEGIN(comp_small);

        // memcpy(matrix, (matrix[0] + beginShiftAmount), (localMatrixSize - beginShiftAmount) * sizeof(int));
        for (colNum = 0; colNum < numColsPerWorker - 1; ++colNum) {
            memcpy(&matrix[colNum * numRows], &matrix[colNum * numRows + beginShiftAmount], (amountToShift) * sizeof(int));
            memcpy(&matrix[colNum * numRows + amountToShift], &matrix[(colNum + 1) * numRows], beginShiftAmount * sizeof(int));
        }
        memcpy(&matrix[localMatrixSize - numRows], &matrix[localMatrixSize - numRows + beginShiftAmount], amountToShift * sizeof(int));

        CALI_MARK_END(comp_small);
        CALI_MARK_END(comp);
    }
    else {
        CALI_MARK_BEGIN(comp);
        CALI_MARK_BEGIN(comp_small);

        // memcpy(matrix, (matrix[0]+beginShiftAmount), (localMatrixSize - beginShiftAmount) * sizeof(int));
        for (colNum = 0; colNum < numColsPerWorker - 1; ++colNum) {
            memcpy(&matrix[colNum * numRows], &matrix[colNum * numRows + beginShiftAmount], (amountToShift) * sizeof(int));
            memcpy(&matrix[colNum * numRows + amountToShift], &matrix[(colNum + 1) * numRows], beginShiftAmount * sizeof(int));
        }
        memcpy(&matrix[(numColsPerWorker-1) * numRows], &matrix[(numColsPerWorker-1) * numRows + beginShiftAmount], amountToShift * sizeof(int));

        CALI_MARK_END(comp_small);
        CALI_MARK_END(comp);
    }

    if (rank != num_procs-1) {
        CALI_MARK_BEGIN(comm);
        CALI_MARK_BEGIN(comm_small);

        MPI_Recv(&matrix[localMatrixSize - beginShiftAmount], beginShiftAmount, MPI_INT, rank + 1, UNSHIFT, MPI_COMM_WORLD, &status);

        CALI_MARK_END(comm_small);
        CALI_MARK_END(comm);
    }
    else {
        CALI_MARK_BEGIN(comp);
        CALI_MARK_BEGIN(comp_small);

        memcpy(&matrix[localMatrixSize - beginShiftAmount], &finalCol[0], beginShiftAmount * sizeof(int));

        CALI_MARK_END(comp_small);
        CALI_MARK_END(comp);
    }

    // if (rank == MASTER) {
    //     printf("\nrank %d\n", rank);
    //     for (int i = 0; i < numColsPerWorker; ++i) {
    //         for (int j = 0; j < numRows; ++j) {
    //             printf("%d,", matrix[i*numRows + j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    if (rank == MASTER) {
        // CALI_MARK_END(unshift);
        // CALI_MARK_END(whole_sort);
        matrix_unshift_time_end = MPI_Wtime();
        matrix_sort_time_end = MPI_Wtime();
        matrix_unshift_time = matrix_unshift_time_end - matrix_unshift_time_start;
        matrix_sort_time = matrix_sort_time_end - matrix_sort_time_start;
        printf("Step 8: Unshifted Matrix\n");
        printf("time: %f\n\n", matrix_unshift_time);

        printf("Matrix Sort Time: %f\n\n", matrix_sort_time);
    }

    // // printf("\nmatrix rank %d\n", rank);
    // // for (int i = 0; i < numColsPerWorker; ++i) {
    // //     for (int j = 0; j < numRows; ++j) {
    // //         printf("%d,", matrix[i][j]);
    // //     }
    // //     printf("\n");
    // // }
    // // printf("\n");

    if (rank == MASTER) {
        matrix_check_time_start = MPI_Wtime();
        // CALI_MARK_BEGIN(sort_check);
    }

    CALI_MARK_BEGIN(correctness_check);

    if (rank % 2 == 1) {
        sorted = 1;
        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            int a = is_sorted(&matrix[colNum * numRows], &matrix[colNum * numRows + numRows]);
            // if (rank == MASTER) { printf("%d, %d,", colNum, a); }
            sorted = sorted && a;
        }
        MPI_Send(&matrix[0], 1, MPI_INT, rank - 1, CHECK_SORTED, MPI_COMM_WORLD);
        MPI_Recv(&isSorted, 1, MPI_INT, rank - 1, CHECK_SORTED, MPI_COMM_WORLD, &status);
        sorted = sorted && isSorted;
        MPI_Send(&sorted, 1, MPI_INT, MASTER, CHECK_SORTED, MPI_COMM_WORLD);
    }
    else {
        isSorted = 1;
        for (colNum = 0; colNum < numColsPerWorker; ++colNum) {
            int a = is_sorted(&matrix[colNum * numRows], &matrix[colNum * numRows + numRows]);
            // if (rank == MASTER) { printf("%d, %d | ", colNum, a); }
            isSorted = isSorted && a;
        }
        MPI_Recv(&sorted, 1, MPI_INT, rank + 1, CHECK_SORTED, MPI_COMM_WORLD, &status);
        if (matrix[localMatrixSize-1] > sorted) {
            sorted = 0;
        }
        else {
            sorted = 1;
        }
        isSorted = isSorted && sorted;
        MPI_Send(&isSorted, 1, MPI_INT, rank + 1, CHECK_SORTED, MPI_COMM_WORLD);
    }

    if (rank == MASTER) {
        sorted = 1;
        for (src = 1; src < num_procs; src += 2) {
            MPI_Recv(&isSorted, 1, MPI_INT, src, CHECK_SORTED, MPI_COMM_WORLD, &status);
            sorted = sorted && isSorted;
        }
        if (sorted) {
            printf("\nmatrix is sorted\n\n");
        }
        else {
            printf("\nmatrix not sorted\n\n");
        }
    }

    CALI_MARK_END(correctness_check);

    if (rank == MASTER) {
        // CALI_MARK_END(sort_check);
        // CALI_MARK_END(whole_computation);
        total_time_end = MPI_Wtime();
        matrix_check_time_end = MPI_Wtime();
        matrix_check_time = matrix_check_time_end - matrix_check_time_start;
        total_time = total_time_end - total_time_start;
        printf("Step 9: Check Matrix\n");
        printf("time: %f\n\n", matrix_check_time);
        printf("Total Time: %f\n\n", total_time);
    }

    delete[] matrix;
    delete[] rowMatrix;
    delete[] finalCol;

    adiak::init(NULL);
    adiak::launchdate();                                  // launch date of the job
    adiak::libraries();                                   // Libraries used
    adiak::cmdline();                                     // Command line used to launch the job
    adiak::clustername();                                 // Name of the cluster
    adiak::value("algorithm", "column_sort");             // The name of the algorithm you are using (e.g., "merge", "bitonic")
    adiak::value("programming_model", "mpi");             // e.g. "mpi"
    adiak::value("data_type", "int");                     // The datatype of input elements (e.g., double, int, float)
    adiak::value("size_of_data_type", sizeof(int));       // sizeof(datatype) of input elements in bytes (e.g., 1, 2, 4)
    adiak::value("input_size", numRows * numCols);        // The number of elements in input dataset (1000)
    adiak::value("input_type", input_type);               // For sorting, this would be choices: ("Sorted", "ReverseSorted", "Random", "1_perc_perturbed")
    adiak::value("num_procs", num_procs);                 // The number of processors (MPI ranks)
    adiak::value("scalability", "strong");                // The scalability of your algorithm. choices: ("strong", "weak")
    adiak::value("group_num", 2);                         // The number of your group (integer, e.g., 1, 10)
    adiak::value("implementation_source", "handwritten"); // Where you got the source code of your algorithm. choices: ("online", "ai", "handwritten").

    mgr.stop();
    mgr.flush();

    MPI_Finalize();
    }