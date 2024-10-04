# CSCE 435 Group project

## 0. Group number: 2

## 1. Group members:
1. Kimberly Chen
2. Spencer Le
3. Suhu Lavu
4. Andrew Mao
5. Jeff Ooi

## 2. Project topic (e.g., parallel sorting algorithms)
Parallel Sorting Algorithms

### 2a. Brief project description (what algorithms will you be comparing and on what architectures)
We will be comparing the following algorithms:
- Bitonic Sort - Kimberly Chen
- Sample Sort -
- Merge Sort -
- Radix Sort -
- Column Sort - Jeff Ooi

We will use the Grace cluster on the TAMU HPRC.

### 2b. Pseudocode for each parallel algorithm
- For MPI programs, include MPI calls you will use to coordinate between processes
- Bitonic Sort Pseudocode
    ```
    func bitonic_sort(matrix, lowIndex, count, direction)
        Initialize MPI
        rank <- MPI rank
        num_procs <- MPI size
    
        if count > 1
            k = count / 2  // Divide the array into two halves
    
            // Sort the first half in ascending order
            bitonicSort(matrix, lowIndex, k, 1)
    
            // Sort the second half in descending order
            bitonicSort(matrix, lowIndex + k, k, 0)
    
            // Merge the two halves according to the 'direction'
            bitonicMerge(arr[], lowIndex, count, direction)
 
            for step = 1 to log2(num_procs)
                partner <- rank XOR step
    
                // Exchange sorted halves with partner process and merge
                if direction == 1 // this means it is ascending
                    if rank < partner
                        MPI send local_array to partner
                        MPI receive partner_array from partner
                        matrix <- bitonic_merge_ascending(matrix, partner_array)
                    else
                        MPI send local_array to partner
                        MPI receive partner_array from partner
                        matrix <- bitonic_merge_descending(matrix, partner_array)
                else // this means it is descending
                    if rank < partner
                        MPI send local_array to partner
                        MPI receive partner_array from partner
                        matrix <- bitonic_merge_descending(matrix, partner_array)
                    else
                        MPI send local_array to partner
                        MPI receive partner_array from partner
                        matrix <- bitonic_merge_ascending(matrix, partner_array)
                end if
            end for
        end if
    
        // Gather all sorted data at master
        if rank is master
            for i: 1 -> num_procs
                MPI receive sorted segments from all processes
            end for
            output "Final sorted array"
        end if
    
        Finalize MPI
    end func
    ```
- Column Sort Pseudocode
    ```
    Assumption: numRows >= 2 * (numCols - 1)^2
    func column_sort(matrix, numRows, numCols)
        Initialize MPI
        rank <- MPI rank
        num_procs <- MPI size
        
        colsPerWorker = numCols / num_procs
        
        if rank is master
            for i: 1 -> num_procs
                MPI send matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1 to process i
            end for
        
        localbuffer <- buffer with enough size to fit matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1
        if rank is not master
            MPI receive matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1 into localbuffer
        
        sort every column in localbuffer
        
        if rank is not master
            MPI send sorted localbuffer to master
        if rank is master
            for i: 1 -> num_procs
            MPI receive all localbuffers into matrix
            end for
        
            transpose matrix
        
            for i: 1 -> num_procs
                MPI send matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1 to process i
            end for
            
        if rank is not master
            MPI receive matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1 into localbuffer
        
        sort every column in localbuffer
        
        if rank is not master
            MPI send sorted localbuffer to master
        if rank is master
            for i: 1 -> num_procs
                MPI receive all localbuffers into matrix
            end for
        
            untranspose matrix

            for i: 1 -> num_procs
                MPI send matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1 to process i
            end for
            
        if rank is not master
            MPI receive matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1 into localbuffer
        
        sort every column in localbuffer

        if rank is not master
            MPI send sorted localbuffer to master
        if rank is master
            for i: 1 -> num_procs
                MPI receive all localbuffers into matrix
            end for
            
            shift matrix
            
            for i: 1 -> num_procs
                MPI send matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1 to process i
            end for
            
        if rank is not master
            MPI receive matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1 into localbuffer
        
        sort every column in localbuffer
        
        if rank is not master
            MPI send sorted localbuffer to master
        if rank is master
            sort last column in matrix
            
            unshift matrix // matrix should now be sorted in Column-Major order
            
            for i: 1 -> num_procs
                MPI send matrix columns i * colsPerWorker to (i + 1) * colsPerWorker -  to process i
            end for
        if rank is not master
            MPI receive matrix columns i * colsPerWorker to (i + 1) * colsPerWorker - 1 into localbuffer
        
        check if localbuffer is sorted
        sorted <- 1 if sorted, 0 if not
        if rank is master
            isSorted <- -1
            MPI reduce sorted to master with minimum sorted into isSorted
            if isSorted is 1
                output "matrix sorted"
            else
                output "matrix not sorted"
        
        Finalize MPI
    end func
    ```

### 2c. Evaluation plan - what and how will you measure and compare
- Input sizes, Input types
- Strong scaling (same problem size, increase number of processors/nodes)
- Weak scaling (increase problem size, increase number of processors)
