import copy

procs = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
exps = [16, 18, 20, 22, 24, 26, 28]
arrayTypes = [0, 1, 2, 3]
fileNames = []

# directory path variables
# directories must be already made
outDirPath = "./grace_files"
inDirPath = "./grace_job_files"

listOfStuffToWrite = []
for exp in exps:
    with open(f"{inDirPath}/mpi_2_{exp}.grace_job", "r") as baseFile:
        listOfStuffToWrite = baseFile.readlines()
    
    # numProcs is bounded by number of columns
    maxProcs = int(listOfStuffToWrite[20].split("=")[-1].strip())
    # print(maxProcs)
    listOfStuffToWriteCopy = copy.deepcopy(listOfStuffToWrite)
    
    # for every proc, generate a corresponding job file
    # cost minimization done here
    for numProcs in procs:
        if (numProcs > maxProcs):
            break
        
        # 48 max cores per node, but numProcs are not divisible by 48
        # divide by 32 instead
        if (numProcs > 32):
            numNodes = int(numProcs / 32)
            listOfStuffToWriteCopy[8] = f"#SBATCH --nodes={numNodes}\n"
            listOfStuffToWriteCopy[9] = f"#SBATCH --ntasks-per-node=32\n"
        else:
            listOfStuffToWriteCopy[8] = f"#SBATCH --nodes=1\n"
            listOfStuffToWriteCopy[9] = f"#SBATCH --ntasks-per-node={numProcs}\n"

        # set processes value and arrayType variable in job file
        # arrayType value will be taken in from command line
        procLine = f"processes={numProcs}\n"
        listOfStuffToWriteCopy[21] = procLine
        listOfStuffToWriteCopy[22] = f"arrayType=$1\n"
        
        with open(f"{outDirPath}/mpi_2_{exp}_{numProcs}.grace_job", "w") as outFile:
            outFile.writelines(listOfStuffToWriteCopy)
        fileNames.append(f"mpi_2_{exp}_{numProcs}.grace_job")

# for every generated file:
# sbatch <filename> <arrayType>
with open(f"{outDirPath}/run.sh", "w") as bashFile:
    bashFile.write("#!/bin/bash\n\n")
    bashFile.write("module load intel/2020b\n")
    bashFile.write("module load CMake/3.12.1\n")
    bashFile.write("module load GCCcore/8.3.0\n")
    bashFile.write("module load PAPI/6.0.0\n\n")
    for file in fileNames:
        # add this line if want to see progress
        # bashFile.write(f'\necho "{file}\n"\n')
        for arrayType in arrayTypes:
            bashFile.write(f"sbatch {file} {arrayType}\n")