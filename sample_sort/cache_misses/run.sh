#!/bin/bash

module load intel/2020b
module load CMake/3.12.1
module load GCCcore/8.3.0
module load PAPI/6.0.0

# For size 2^26 (67108864)
sbatch mpi.grace_job 67108864 16 random
sbatch mpi.grace_job 67108864 32 random
sbatch mpi.grace_job 67108864 16 sorted
sbatch mpi.grace_job 67108864 32 sorted
sbatch mpi.grace_job 67108864 16 reverse_sorted
sbatch mpi.grace_job 67108864 32 reverse_sorted

# For size 2^25 (33554432)
sbatch mpi.grace_job 33554432 16 random
sbatch mpi.grace_job 33554432 32 random
sbatch mpi.grace_job 33554432 16 sorted
sbatch mpi.grace_job 33554432 32 sorted
sbatch mpi.grace_job 33554432 16 reverse_sorted
sbatch mpi.grace_job 33554432 32 reverse_sorted