import copy

# Define only the missing combinations
missing_cases = [
    # p32 2^28 cases
    (32, 28, "random"),
    (32, 28, "sorted"),
    (32, 28, "reverse_sorted"),
    (32, 28, "perturbed"),

    # p64 2^28 cases
    (64, 28, "random"),
    (64, 28, "sorted"),
    (64, 28, "reverse_sorted"),
    (64, 28, "perturbed"),

    # p128 2^28 cases
    (128, 28, "random"),
    (128, 28, "sorted"),
    (128, 28, "reverse_sorted"),
    (128, 28, "perturbed"),

    # p256 2^28 cases
    (256, 28, "random"),
    (256, 28, "sorted"),
    (256, 28, "reverse_sorted"),
    (256, 28, "perturbed"),

    # p512 2^28 cases
    (512, 28, "random"),
    (512, 28, "sorted"),
    (512, 28, "reverse_sorted"),
    (512, 28, "perturbed"),

    # p1024 2^28 cases
    (1024, 28, "random"),
    (1024, 28, "sorted"),
    (1024, 28, "reverse_sorted"),
    (1024, 28, "perturbed")
]

# directory path variables
outDirPath = "./job_files"
inDirPath = "./job_files/templates"

fileNames = []
for numProcs, exp, arrayType in missing_cases:
    with open(f"{inDirPath}/mpi_2_{exp}.grace_job", "r") as baseFile:
        templateContent = baseFile.readlines()
    
    array_size = 2 ** exp
    
    # Make a fresh copy of the template
    listOfStuffToWriteCopy = copy.deepcopy(templateContent)
    
    # Set nodes and tasks per node
    if numProcs > 32:
        numNodes = numProcs // 32
        tasks_per_node = 32
    else:
        numNodes = 1
        tasks_per_node = numProcs
    
    # Update the SBATCH directives with fixed time and memory
    listOfStuffToWriteCopy[7] = f"#SBATCH --time=00:10:00           #Set the wall clock limit to 10 minutes\n"
    listOfStuffToWriteCopy[8] = f"#SBATCH --nodes={numNodes}\n"
    listOfStuffToWriteCopy[9] = f"#SBATCH --ntasks-per-node={tasks_per_node}\n"
    listOfStuffToWriteCopy[10] = f"#SBATCH --mem=64G                 #Request 64GB memory per node\n"

    # Set array_size, processes and input_type
    listOfStuffToWriteCopy[20] = f"array_size={array_size}\n"
    listOfStuffToWriteCopy[21] = f"processes={numProcs}\n"
    listOfStuffToWriteCopy[22] = f"input_type=\"{arrayType}\"\n"
    
    # Modify CALI_CONFIG and mpirun lines
    listOfStuffToWriteCopy[28] = f"CALI_CONFIG=\"spot(output=p{numProcs}-a{array_size}-i{arrayType}.cali, time.variance,profile.mpi)\" \\\n"
    listOfStuffToWriteCopy[29] = f"mpirun -np $processes ./sample_sort_mpi $array_size $input_type\n"

    filename = f"mpi_2_{exp}_{numProcs}_{arrayType}.grace_job"
    with open(f"{outDirPath}/{filename}", "w") as outFile:
        outFile.writelines(listOfStuffToWriteCopy)
    fileNames.append(filename)

# Create run.sh
with open(f"{outDirPath}/run.sh", "w") as bashFile:
    bashFile.write("#!/bin/bash\n\n")
    bashFile.write("module load intel/2020b\n")
    bashFile.write("module load CMake/3.12.1\n")
    bashFile.write("module load GCCcore/8.3.0\n")
    bashFile.write("module load PAPI/6.0.0\n\n")
    for file in fileNames:
        bashFile.write(f"sbatch {file}\n")