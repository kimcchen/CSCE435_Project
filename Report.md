# CSCE 435 Group project

## 0. Group number: 2 (Honors)

## 1. Group members:
1. Kimberly Chen
2. Spencer Le
3. Suhu Lavu
4. Andrew Mao
5. Jeff Ooi

We will communicate via a messaging group chat.

## 2. Project topic (e.g., parallel sorting algorithms)
Parallel Sorting Algorithms

### 2a. Brief project description (what algorithms will you be comparing and on what architectures)
We will be comparing the following algorithms:
- Bitonic Sort - Kimberly Chen
  - Bitonic Sort is a parallel sorting algorithm well-suited for distributed systems where the number of processors is a power of 2. The idea behind bitonic sorting is to iteratively      build "bitonic sequences" (sequences that are first increasing and then decreasing) and then merge them to make the sorted array. The algorithm works as follows:
    1. Locally sorting the processor's assigned portion of the array based on the processor's rank. Processors with even ranks sort in increasing order, while those with odd ranks 
       sort in decreasing order.
    2. Bitonic merging where for each pair of processors, one sends its data to the other, and they compare their values. If the rank of the processor is less than its partner's, it 
       keeps the smaller values; otherwise, it keeps the larger values.
    3. Recursive sorting where processors again sort their local arrays, but now based on the next bit of their rank so that the sorting order alternates between ascending and 
       descending for different processors at each level
    4. Iterative merging where at each depth, the bitonic sequences are merged to form longer sorted sequences, with alternating sorting directions.
    5. Final gathering where the sorted subarrays held by each processor are gathered by the root processor.
  - The algorithm requires p = 2^k processors and for the number of elements to be 2^n. The time complexity of bitonic sort is O(n log^2 n), but since it is distributed among p 
    processors, the parallel time complexity becomes O((n log^2 n) / p). 
- Sample Sort - Spencer Le
  - Sample Sort is a highly scalable parallel sorting algorithm that divides the input data among multiple processors and efficiently sorts large datasets using a combination of local sorting, sampling, and redistribution. It is particularly well-suited for distributed-memory systems. The algorithm works as follows:
    1. Local Sorting: Each processor receives an equal-sized portion of the array to sort. The processors independently sort their local portions using a sequential sorting algorithm (e.g., quicksort). This step reduces the problem size locally before redistribution.
    2. Sampling: After the local sort, each processor selects a small sample of its locally sorted data. These samples are used to choose “splitters,” which will define the boundaries between the ranges of values to be assigned to each processor.
    3. Choosing Splitters: The samples from all processors are gathered on the root processor (or another dedicated processor), which sorts them and selects a set of splitters. These splitters divide the global array into buckets, where each bucket represents a range of values that should be sent to a specific processor.
    4. Redistribution (All-to-All Communication): Using the splitters, each processor sends its locally sorted data to the appropriate processor, ensuring that the data in each bucket is sent to the processor responsible for that range of values. This step involves an all-to-all communication to exchange data between processors.
    5. Final Sorting and Gathering: Once each processor has received its portion of data, it performs a final sort on the received data, ensuring that the data within each processor is fully sorted. Finally, the sorted subarrays are gathered by the root processor to form the final globally sorted array.
   - The algorithm is highly scalable because it minimizes the amount of inter-processor communication and focuses on local sorting and efficient data redistribution. The time complexity of Sample Sort is `O(nlog(n)/p)`, where `n` is the number of elements, and `p` is the number of processors
- Merge Sort - Suhu Lavu
  - Merge sort is a parallel sorting algorithm that uses a divide and conquer approach to sort an array by recursively dividing it into subarrays, sorting those subarrays, and then merging them backk together. The algorithm works as follows:
    1. The initial unsorted array is divided into smaller subarrays. Each processor is assigned a portion of this array.
    2. Each processor sorts their portion using a sequential merge sort where the array is divided into subarrays, sorted, and merged back together.
    3. Once the local arrays are sorted, processors merge their sorted subarrays. Merging involves comparing elements between 2 sorted arrays and combining them into a single sorted array. 
    4. Sorted arrays from multiple processors are merged together in pairs, doubling the size of the sorted array at each step until there is a single sorted array with all elements.
    5. The root process gathers the final sorted array.
  - Sequential merge sort has a time complexity of `O(nlog(n))`. However, because this is distributed among `p` processors, the time complexity becomes `O(nlog(n)/p)` in parallel.
- Radix Sort - Andrew Mao
  - Radix sort is a non comparative sorting algorithm that sorts an array through evaluating and counting individual digits. We sort elements into buckets with a helper counting sort function. This is repeated multiple times starting from the least significant digit until the most significant digit. In contrast to serial radix sort, we must shuffle elements between processors' local buckets. The algorithm works as follows:
    1. Divide the whole unsorted array into local subarrays to be sorted by p processors.
    2. For each digit, from the least significant to the most significant, of the largest number in the array, run serial counting sort for the local array.
    4. Gather counts from all processes to calcaulte the prefix sums and total sums to find the position of each digit like in counting sort.
    5. Move elements within each processor to appropriate local buckets
    6. Send and receive elements to correct positions in local arrays of other processes.
    7. Gather all sorted local arrays back to master.
   - Sequential run time of decimal radix sort is `O(n * d)` where d is the max number of digits or `log10(maxVal)` and n is the number of elements in the array. An ideal parallel runtime would divide the time by p processors, which is `O(n * d / p)`
- Column Sort - Jeff Ooi
  - Column Sort is an eight step matrix parallel sorting algorithm, with the eight steps of the algorithm being as follows:
    1. Sort each column of the matrix.
    2. Transpose the matrix by reading the elements in Column-Major and writing back to the matrix in Row-Major, keeping the original shape and dimensions of the matrix
    3. Sort each column of the transposed matrix.
    4. Untranspose the matrix by reading the elements in Row-Major and writing back to the matrix in Column-Major, keeping the original shape and dimensions of the matrix.
    5. Sort each column of the untransposed matrix.
    6. Shift the elements of the matrix column-wise down by the floor of rows/2 and fill in the empty spaces in the first column with negative infinity. This will create a rows x (columns + 1) matrix. Fill in the remaining spaces in the final column with positive infinity.
    7. Sort each column of the shifted matrix.
    8. Unshift the elements of the matrix by deleting the infinities and shifting every element column-wise up by the floor of rows/2. This will return the matrix back to the original rows x columns dimensions and complete the sort.
  - The algorithm has the restriction that `rows >= 2 * (columns - 1) ^ 2` and `rows % cols == 0`. The longest step of the algorithm is the sorting of the columns, which runs in `O(nlog(n))`. However, because those steps are distributed among `p` processes, the runtime of the algorithm is `O(nlog(n)/p)` time, where `n` is the total number of elements in the matrix.

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
- Sample Sort Pseudocode
    ```
    func parallel_sample_sort(local_data):
        // Helper Functions
        function initialize_array(array_size, input_type):
            array = new array[array_size]
            
            if input_type == "random":
                for i = 0 to array_size - 1:
                    array[i] = random_number()
                    
            else if input_type == "sorted":
                for i = 0 to array_size - 1:
                    array[i] = i
                    
            else if input_type == "reverse_sorted":
                for i = 0 to array_size - 1:
                    array[i] = array_size - i - 1
                    
            else if input_type == "perturbed":
                // Create sorted array first
                for i = 0 to array_size - 1:
                    array[i] = i
                // Perturb 1% of elements
                num_perturb = array_size / 100
                for i = 0 to num_perturb - 1:
                    idx1 = random(0, array_size - 1)
                    idx2 = random(0, array_size - 1)
                    swap(array[idx1], array[idx2])
            
            return array

        function get_bucket_index(value, splitters):
            bucket = 0
            while bucket < splitters.length AND value > splitters[bucket]:
                bucket += 1
            return bucket

        // Main Algorithm
        procedure parallel_sample_sort(array_size, input_type):
            // Initialize MPI
            rank = MPI_Get_rank()
            num_procs = MPI_Get_size()
            root = 0
            
            // Validate input
            if array_size < num_procs:
                print "Error: Array size must be >= number of processors"
                return
            
            // Calculate local sizes
            base_size = array_size / num_procs
            remainder = array_size % num_procs
            local_size = base_size + (rank < remainder ? 1 : 0)
            
            // Initialize and distribute data
            if rank == root:
                input_array = initialize_array(array_size, input_type)
                
            local_array = new array[local_size]
            MPI_Scatter(input_array, local_size, local_array)
            
            // Local sort
            sort(local_array)
            
            // Sample selection
            samples_per_proc = 3
            local_samples = new array[num_procs - 1]
            if local_size > 0:
                stride = local_size / num_procs
                for i = 0 to num_procs - 2:
                    idx = min(floor(stride * (i + 1)), local_size - 1)
                    local_samples[i] = local_array[idx]
            else:
                fill(local_samples, MAX_INT)
            
            // Gather samples and select splitters
            if rank == root:
                global_samples = new array[num_procs * (num_procs - 1)]
                
            MPI_Gather(local_samples, num_procs - 1, global_samples)
            
            if rank == root:
                sort(global_samples)
                // Select evenly spaced splitters
                for i = 0 to num_procs - 2:
                    idx = ((i + 1) * global_samples.length) / num_procs
                    local_samples[i] = global_samples[idx]
                    
            // Broadcast splitters to all processes
            MPI_Broadcast(local_samples)
            
            // Distribute elements to buckets
            temp_buckets = array of num_procs empty arrays
            for each element in local_array:
                bucket = get_bucket_index(element, local_samples)
                temp_buckets[bucket].append(element)
            
            // Prepare for all-to-all communication
            max_bucket_size = local_size + 1
            buckets = new array[max_bucket_size * num_procs]
            
            for i = 0 to num_procs - 1:
                buckets[i * max_bucket_size] = temp_buckets[i].length
                copy(temp_buckets[i], buckets[i * max_bucket_size + 1])
            
            // All-to-all exchange
            bucket_buffer = new array[max_bucket_size * num_procs]
            MPI_Alltoall(buckets, max_bucket_size, bucket_buffer)
            
            // Process received data
            merged_data = empty array
            for i = 0 to num_procs - 1:
                count = bucket_buffer[i * max_bucket_size]
                start = i * max_bucket_size + 1
                end = start + count
                merged_data.append(bucket_buffer[start:end])
            
            // Sort merged data
            sort(merged_data)
            
            // Final gather
            if rank == root:
                final_array = new array[array_size]
                
            MPI_Gather(merged_data, merged_data.length, final_array)
            
            // Verify sort (only root)
            if rank == root:
                is_sorted = true
                for i = 1 to array_size - 1:
                    if final_array[i - 1] > final_array[i]:
                        is_sorted = false
                        break
                        
                print is_sorted ? "Array sorted correctly" : "Sort failed"
            
            MPI_Finalize()
    end func
    ```
- Merge Sort Pseudocode
    ```
    func parallel_merge_sort(array, n):
        Initialize MPI
        rank <- MPI rank
        size <- MPI size

        # calculate size of worker array and create local copy
        worker_size = n / size
        worker_arr = new array[worker_size]

        # send real array to all processes
        MPI_Scatter(array, worker_size, MPI_INT, worker_arr, worker_size, MPI_INT, root=0, MPI_COMM_WORLD)

        # sort local copy
        sorted_arr = merge_sort(worker_arr)

        # merge sorted pieces back together
        step = 1
        while step < size:
            if rank % (2 * step) == 0:
                if rank + step < size:
                    # get sorted array from another process
                    received_size = worker_size * step
                    received_array = new array[received_size]
                    MPI_Recv(received_array, received_size, MPI_INT, source=rank + step, tag=0, MPI_COMM_WORLD)

                    # merge received array with local array
                    sorted_arr = merge(sorted_arr, received_array)
                    worker_size += received_size
            else:
                # send sorted array to another process
                target = rank - step
                MPI_Send(sorted_arr, worker_size, MPI_INT, dest=target, tag=0, MPI_COMM_WORLD)
                break out of loop
            step *= 2

        # master process gathers the final sorted array
        if rank == 0:
            MPI_Gather(sorted_arr, worker_size, MPI_INT, array, worker_size, MPI_INT, root=0, MPI_COMM_WORLD)

        Finalize MPI
    end func

    func merge(arr1, arr2):
        size_1 = size(arr1)
        size_2 = size(arr2)
        final = new array[size_1 + size_2]
        i, j, k = 0

        while i < size_1 and j < size_2:
            if arr1[i] < arr2[j]:
                final[k] = arr1[i]
                i += 1
            else:
                final[k] = arr2[j]
                j += 1
            k += 1
        
        # copy leftover elements
        while i < size_1:
            final[k] = arr1[i]
            i += 1
            k += 1

        while j < size_2:
            final[k] = arr2[j]
            j += 1
            k += 1
        
        return final
    end func

    func merge_sort(arr):
        # base case
        if size(arr) <= 1:
            return arr

        # find midpoint
        mid = size(arr) / 2

        # sort left array        
        left_sorted = merge_sort(arr[0:mid])

        # sort right array
        right_sorted = merge_sort(arr[mid + 1:])

        # merge them together
        return merge(left_sorted, right_sorted)
    end func
    ```

- Radix Sort Pseudocode
    ```
    func radix_sort(matrix, lowIndex, count, direction)
        Initialize MPI
        rank <- MPI rank
        num_procs <- MPI size

            keysPerWorker <- numKeys / num_procs
            
            if rank is master: 
                for i: 1 -> num_procs - 1:
                    MPI send arr[i * keysPerWorker to (i + 1) * keysPerWorker - 1] to process i
                end for
            end if
            
            localbuffer <- [keysPerWorker]
            
            if rank is worker:
                MPI receive arr[i * keysPerWorker to (i + 1) * keysPerWorker - 1] into localbuffer
            end if
            
            # For each iteration, process g bits at a time (Radix Sorting)
            max_bits <- get bits max(arr)
            for i: 0 -> (max_bits / g) do:
                local_bucket_counts <- counting_sort(localbuffer, g, i)
        
                # Aggregate the result of counting sort, make global count of keys per bucket
                global_bucket_counts <- MPI all_to_all(local_bucket_counts)  
                prefix_sums <- get prefix sums of global_bucket_counts
                bucketed_keys <- distribute_keys(localbuffer, g, i, prefix_sums)

                # update these results to preserve order
                exchanged_keys <- MPI all_to_all_exchange(bucketed_keys)
                sort_by_g_bits(exchanged_keys, g, i)
                localbuffer <- exchanged_keys
            end for
            
            if rank is worker:
                MPI send localbuffer to master
            end if
            
            if rank is master: 
                for i: 1 -> num_procs - 1 do:
                    MPI receive sorted buffer into arr
                end for
        
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
We will measure and compare the sorting times of:
- Different input sizes and input types
    - Input sizes are 2^16, 2^18, 2^20, 2^22, 2^24, 2^26, and 2^28
    - Input types are sorted, random, reverse sorted, and 1% perturbed
- Strong scaling (same problem size, increase number of processors)
    - Number of processes are 2, 4, 8, 16, 32, 64, 128, 256, 512, and 1024
- Weak scaling (increase problem size, increase number of processors)
    - Number of processes are 2, 4, 8, 16, 32, 64, 128, 256, 512, and 1024

We will collect them using Caliper and compare them using Thicket.

## 3. Caliper Instrumentation
### 3a. Calltrees
- Bitonic Sort Calltree  
    <img width="473" alt="image" src="https://github.com/user-attachments/assets/61c303aa-33b2-42df-8f0d-43bc1dd63c71">


- Sample Sort Calltree
  
  <img width="755" alt="Screenshot 2024-10-23 at 12 16 30 AM" src="https://github.com/user-attachments/assets/6b935be9-6627-4d6e-b87e-0fab49d38ab1">
  
- Merge Sort Calltree  
    ![Merge Sort Calltree](https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/calltree.png?raw=true)  

- Radix Sort Calltree  
  <img width="708" alt="Screenshot 2024-10-22 at 5 59 43 PM" src="https://github.com/user-attachments/assets/91253424-0569-49a0-b2be-0d585fac01b1">
- Column Sort Calltree  
    `data_init_X` is `data_init_runtime`  
    ![Column Sort Calltree](https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/column_sort_calltree.png?raw=true)

### 3b. Metadata
- Bitonic Sort Metadata  
    <img width="827" alt="Screenshot 2024-10-22 at 11 47 22 PM" src="https://github.com/user-attachments/assets/26110a27-cba5-4047-b591-db9ac2cfa95e">

- Sample Sort Metadata
  <img width="991" alt="Screenshot 2024-10-23 at 12 40 27 AM" src="https://github.com/user-attachments/assets/52fb003a-7d18-459d-be3b-d20caa39f85b">


- Merge Sort Metadata
    ![Merge Sort Metadata](https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/metadata.png?raw=true)

- Radix Sort Metadata  
  <img width="912" alt="Screenshot 2024-10-22 at 11 25 50 PM" src="https://github.com/user-attachments/assets/86e0e3b3-5839-44b1-9945-33d454271993">
- Column Sort Metadata  
    <img alt="Column Sort Metadata" src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/column_sort_metadata.png?raw=true">

## 4. Performance Evaluation

- Bitonic Sort:
  <img width="1358" alt="image" src="https://github.com/user-attachments/assets/bf0dec61-25b5-4a4a-8938-961166a13d00">
  The graph indicates that the speedup maxes out around 512 processes, with the speedup become significantly higher at around 128 processes. We can see from this that using 512 processes would be the most efficient for sorting a randomized array size of 2^16.
- Sample Sort:
  ![Screenshot 2024-11-04 at 10 52 20 AM](https://github.com/user-attachments/assets/79f364f6-5c16-4d16-9376-f8e8040046e0)

  - Strong scaling speedup for different input types with problem size 2^26. Sample Sort shows good scaling up to 64 processes with sorted arrays achieving the best speedup. Performance degrades beyond 64 processes due to communication overhead.

- Merge Sort  
    ![Speedup Graph](https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/figs/2^26/ReverseSorted/main_speedup.png?raw=true)  
    Here is an example of how parallelization can improve performance. The graph shows the speedup of the whole program for an input size of 2^26 and a reverse sorted input. We can see that increasing the number of processes allows for drastic improvements of up to 10x in runtime up to 128 processes, but then results in diminishing returns. More analysis and plots are available in the merge_sort folder.

- Radix Sort:
  - Limitations with the algorithm: Innefficiencies with the current implementation resulted in limited data for higher number of processors. Because of the nature of radix sort using counting sort under the hood, there is a lot of potential for uneven load balances between processors. (There can be more of one digits and a skewed distribution of buckets)  
  ![Whole Comp Time vs Processors](https://github.com/user-attachments/assets/a41b43d3-023d-4cc2-855f-c89532e7c3fe)  
  The graph indicates that the current algorithm is strongly scaled to a certain point where we see diminishing returns at 512+ processors. The high runtime can be attributed to the current inefficient communication.  
  ![Large Comm time](https://github.com/user-attachments/assets/7dd51ee5-6c2c-4ed5-9ca7-b499eae5d19a)  
  *Labeled incorrecly as small* 
- Column Sort  
    ![Strong Scaling Graph](https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots/Strong%20Scaling%20of%20main%20for%20sorting%202^28%20Random%20elements.jpg?raw=true)  
    The graph indicates that the algorithm is strongly scaled, with a bound of about 2 seconds due to the sequential runtime.  
    ![Speedup Graph](https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots/Speedup%20of%20main%20for%20sorting%202^28%20Random%20elements.jpg?raw=true)  
    The graph indicates that the speedup maxes out at 512 processes, with about a 40x speedup compared to 2 processes. Because 512 is the maximum number of processes that can be used to sort 2^28 elements due to the restrictions, this means that using 512 processes would be the most efficient.  
    More analysis and plots in the column_sort folder.

## 5. Presentation

- Bitonic Sort<br/>
  - <strong>comm Strong Scaling Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscaling16.png" width=33% alt="Strong Scaling Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscaling18.png" width=33% alt="Strong Scaling Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscaling20.png" width=33% alt="Strong Scaling Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscaling22.png" width=33% alt="Strong Scaling Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscaling24.png" width=33% alt="Strong Scaling Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscaling26.png" width=33% alt="Strong Scaling Graph comm 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscaling28.png" width=33% alt="Strong Scaling Graph comm 2^28"><br/>
        description
        <br/>
        - <strong>comm Speedup Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomm16.png" width=33% alt="Speedup Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomm18.png" width=33% alt="Speedup Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomm20.png" width=33% alt="Speedup Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomm22.png" width=33% alt="Speedup Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomm24.png" width=33% alt="Speedup Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomm26.png" width=33% alt="Speedup Graph comm 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomm28.png" width=33% alt="Speedup Graph comm 2^28"><br/>
        As seen in the Strong Scaling graphs, the smaller input sizes (2^16, 2^18, 2^20) have a decreasing Speedup as the number of processes increases. This is due to the communication overhead overtaking the benefit of using more processes, so to minimize communication time, it would be best to use 2 processes for those input types. For the rest of the graphs, the graphs appear random. This is due to the Strong Scaling graph being essentially constant, so any deviation from the average will cause a speedup or slowdown to be calculated.
        <br/>
        - <strong>comp_large Strong Scaling Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscalingcomplarge16.png" width=33% alt="Strong Scaling Graph comp_large 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscalingcomplarge18.png" width=33% alt="Strong Scaling Graph comp_large 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscalingcomplarge20.png" width=33% alt="Strong Scaling Graph comp_large 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscalingcomplarge22.png" width=33% alt="Strong Scaling Graph comp_large 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscalingcomplarge24.png" width=33% alt="Strong Scaling Graph comp_large 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscalingcomplarge26.png" width=33% alt="Strong Scaling Graph comp_large 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/strongscalingcomplarge28.png" width=33% alt="Strong Scaling Graph comp_large 2^28"><br/>
        All the graphs show an inversely proprotional relationship between the comp_large time and the number of processes for every input type. Interestingly, the sorted input type is faster for every input size. This is due to the algorithm using std::sort, which uses an implementation of quicksort. Since the input is already sorted, the quicksort algorithmn does not need to make as many, if any, swaps, so it runs quicker. Additionally, for the graphs that have this data point, the graphs level off at about 128 processes, indicating that there is not much more that can be parallelized and sequential runtime is now the limiting factor. These graphs also indicate that the comp_large scales strongly with the number of processes.
        <br/>
        - <strong>comp_large Speedup Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomplarge16.png" width=33% alt="Speedup Graph comp_large 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomplarge18.png" width=33% alt="Speedup Graph comp_large 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomplarge20.png" width=33% alt="Speedup Graph comp_large 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomplarge22.png" width=33% alt="Speedup Graph comp_large 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomplarge24.png" width=33% alt="Speedup Graph comp_large 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomplarge26.png" width=33% alt="Speedup Graph comp_large 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/bitonic_sort/plots/speedupcomplarge28.png" width=33% alt="Speedup Graph comp_large 2^28"><br/>
        Since the Strong Scaling graphs have an inversely proportional relationship between the comp_large time and the number of processes, we would expect a proportional relationship between the comp_large Speedup and the number of processes, which is shown for every graph. Due to the extremely strong strong scaling of comp_large, using the maximum number of processes will lead to the largest speedup, which is ideal.
        <br/>

- Sample Sort<br/>

- Merge Sort<br/>
    - <strong>comm Strong Scaling Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^16/comm.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^18/comm.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^20/comm.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^22/comm.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^24/comm.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^26/comm.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^26"><br/>
        For smaller input sizes, we see that the communication times are very variable. This can be attributed to the fact that when there is a smaller volume of data, communication is relatively constant, as seen by the scale of the y-axis. However, as we increase the input array size, we see that as the number of processes increase, the communication time actually decreases. This is a result of my implementation of merge sort. In the communication step of merge sort, different processes are merging together different parts of the array. This merge computation dominates the total communication time and because the merge computation decreases as the number of processes increase, the average communication time per rank decreases as the number of processes increases. Furthermore, in general, the type of input array (random, sorted, etc.) does not have any effect on communication time.
        <br/>
    - <strong>comm Speedup Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^16/comm_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^18/comm_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^20/comm_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^22/comm_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^24/comm_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^26/comm_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^26"><br/>
        Similar to the comm strong scaling plots above, the smaller input sizes have decreasing speedup as we increase the array size. This is a result of communication overhead when the array size is small, reducing the benefit of parallelization. As you increase the array size, the communication speedup tends to increase as the number of processes increase. This is due to the same reasons explained in the strong scaling section. 
        <br/>
    - <strong>comm Weak Scaling Graph</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/weak_scaling/comm.png?raw=true" width=33% alt="Weak Scaling Graph comm 2^16"><br/>
        The ratio used to generate the weak scaling graph is 4:2, meaning that as input size quadruples, the number of processes doubles. Therefore, in theory, we would obtain a line with a slope of 2. The graphs do tend to follow this general trend, with a small deviation at 64 processes. This deviation can be explained by the fact that at 64 processes, Grace requires 2 nodes which would increase communication times non-proportionally. 
        <br/>

    - <strong>comp_large Strong Scaling Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^16/comp_large.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^18/comp_large.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^20/comp_large.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^22/comp_large.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^24/comp_large.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^26/comp_large.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^26"><br/>
        For the computation time, the trend is the same for all input types and array sizes. As we increase the number of processes, the total computation time decreases up until a certain point due to diminishing returns. For all combinations of inputs, the optimal number of processes appears to be 128 processes when the curve flattens out. The reason why the input type of the array does not change the sorting time is because of my implementation of merge sort. Because I used a naive implementation that does not check to see if pieces of the input are already sorted, the time it takes to sort is the same regardless of how sorted the array is.
        <br/>
    - <strong>comp_large Speedup Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^16/comp_large_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^18/comp_large_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^20/comp_large_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^22/comp_large_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^24/comp_large_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^26/comp_large_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^26"><br/>
        The comp_large speedup graphs match what we see in the strong scaling graphs above. As the number of processes increase, we get a larger speedup, demonstrating the benefits of parallelization. This trend is the same across all array sizes and input types. 
        <br/>
    - <strong>comp_large Weak Scaling Graph</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/weak_scaling/comp_large.png?raw=true" width=33% alt="Weak Scaling Graph comm 2^16"><br/>
        The ratio used to generate the weak scaling graph is 4:2, meaning that as input size quadruples, the number of processes doubles. Therefore, in theory, we would obtain a line with a slope of 2. We can see that the lines in general do follow this trend with slight deviations at 64 processes again.
        <br/>

    - <strong>main Strong Scaling Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^16/main.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^18/main.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^20/main.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^22/main.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^24/main.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^26/main.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^26"><br/>
        We observe that for smaller input types, the time it takes to perform the entire sort actually increases as the number of processes increase. This is due to the fact that at smaller array sizes, the communication times dominate the benefits gained from parallelization, causing the runtime to increase rather than improve. However, as we increase the array size, we see that the time it takes to perform the entire algorithm does in fact decrease as the number of processes increase. I did have a few bad runs on Grace for an array size of 2^26 which corresponds to the outliers seen on that plot. Excluding those, we can see that the optimal number of processes appears to be around 128 as that is when the curve flattens out, indicating that adding any extra processors beyong this point provides no additional benefit. Furthermore, all array types have similar trends and runtimes due to my implementation of merge sort not factoring in the "sortedness" of the initial array.
        <br/>
    - <strong>main Speedup Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^16/main_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^18/main_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^20/main_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^22/main_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^24/main_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/2^26/main_speedup.png?raw=true" width=33% alt="Strong Scaling Graph comm 2^26"><br/>
        The speedup plots match what we expect to see based on the strong scaling graphs. For smaller array sizes, the speedup actually decerases as the number of processes increase. However, for large array sizes, the speedup increases as the number of processes increase up until a certain point. As we increase the array size, the max speedup also increases at a higher number of processes. This makes sense as larger amounts of data have more potential for parallelization and require more threads to achieve this max potential. Ignoring the outliers, the different array types also do follow the same trend.
        <br/>
    - <strong>main Weak Scaling Graph</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/pres_figs/weak_scaling/main.png?raw=true" width=33% alt="Weak Scaling Graph comm 2^16"><br/>
        The ratio used to generate the weak scaling graph is 4:2, meaning that as input size quadruples, the number of processes doubles. Therefore, in theory, we would obtain a line with a slope of 2. The lines do follow this trend in general, again with slight deviations around 64 processes which can be attributed to communication times spiking when increasing the number of nodes on Grace.
        <br/>

    

    - <strong>Cache Miss Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/plots_cm/Total%20L1%20misses.jpg?raw=true" width="33%" alt="L1 Misses Graph">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/merge_sort/plots_cm/Total%20L2%20misses.jpg?raw=true" width="33%" alt="L2 Misses Graph"><br/>
        The behavior for the cache miss plots is as expected. We can see two separate sets of parallel lines on both plots that indicate that as the input size increase the number of cache misses also increase. Furthermore, we can also see that as the number of processes increase, the number of cache misses also increases. Both these observations align with expectations as with more data and more threads accessing that data, the number of cache misses will increase. Finally, we also see, based on the scale of the plots, that the number of L2 cache misses is greater than the number of L1 cache misses, which makes sense as the L1 cache is accessed before the L2 cache and has less memory.

- Radix Sort<br/>
- Radix Sort<br/>
 - <strong>comm Strong Scaling Graphs</strong><br/>
       <img width="33%" alt="commStrongScale2_28" src="https://github.com/user-attachments/assets/4b9a6872-4f17-41f9-9405-7a4a50161787">
  <img width="33%" alt="commStrongScale2_26" src="https://github.com/user-attachments/assets/cfb90f65-8435-4af9-8a60-2cdbeba4e126">
  <img width="33%" alt="commStrongScale2_24" src="https://github.com/user-attachments/assets/705db3ad-7d45-4b57-a556-d455eb1b3815">
  <img width="33%" alt="commStrongScale2_22" src="https://github.com/user-attachments/assets/295631ea-2360-4cc4-afcb-a064bdc66f4a">
  <img width="33%" alt="commStrongScale2_20" src="https://github.com/user-attachments/assets/70f0b450-9671-47c3-a7c3-faeb5c44e606">
  <img width="33%" alt="commStrongScale2_18" src="https://github.com/user-attachments/assets/fb4b9fb5-8ae1-4030-b9da-c167041cee06">
  <img width="33%" alt="commStrongScale2_16" src="https://github.com/user-attachments/assets/f7a1e24e-04ee-4c10-99b4-f932f46ab9f0">
      <br/>
     For lower array sizes such as 2^16 and 2^18, we notice a positive trend where as the number of processors increases, the communication time also increases. For array sizes that are higher that 2^18, we notice a plateau at 32 processors followed by a dip in communication time with a minor increase again for the highest number of processors. This dip could be due to the diminishing amount of data that each individual processor manages. Even when the volume of communication increases (more processors), there is less information passed between processors. The randomly sorted array dominates the total communication time, this is likely due to the fact that I had generated random with a normal distribution from (0, INT32_MAX) instead of capping the max value as the array_size. This makes sense as radix sort's runtime is affected mainly by array_size and the number of digits in each number. 
     <br/>
- <strong>comm Speedup Graphs</strong><br/>
      <img width="33%" alt="commSpeedup2_28" src="https://github.com/user-attachments/assets/7fea8def-679e-42df-857d-da55d1f7347e">
<img width="33%" alt="commSpeedup2_26" src="https://github.com/user-attachments/assets/b026b88a-9024-4a45-939a-5dfdf260b9fb">
<img width="33%" alt="commSpeedup2_24" src="https://github.com/user-attachments/assets/0380814e-6643-4a39-9dd1-039518d93c55">
<img width="33%" alt="commSpeedup2_22" src="https://github.com/user-attachments/assets/09640cd7-93d2-4ed9-82d7-e49649a53379">
<img width="33%" alt="commSpeedup2_20" src="https://github.com/user-attachments/assets/acfa7f0e-8f61-43c0-8746-33e9d407fd72">
<img width="33%" alt="commSpeedup2_18" src="https://github.com/user-attachments/assets/e37a1c66-d30c-4de3-bdc2-97ab889a94ce">
<img width="33%" alt="commSpeedup2_16" src="https://github.com/user-attachments/assets/b09f80f0-9cd1-41b6-9fc2-a52ca7de4927">
       <br/>
        For lower array sizes such as 2^16 and 2^18, we notice a negative trend where as the number of processors increases, the speed up decreases. This makes sense as for smaller array size, the communication overhead from more and more processors dominates because smaller sorting problems won't benefit as much from parallelism. For array sizes that are higher that 2^18, we notice a plateau at 32 processors followed by a surge in speed up at 512 processors. This plateau corresponds to the patterns seen in the strong scaling graphs. This speedup postive rate of change at 32 processors and on corresponds to a small turning point in decreasing communication times in the. This increase could be due to the diminishing amount of data that each individual processor manages. Even when the volume of communication increases (more processors), there is less information passed between individual processors which can lead to quicker communication. 
        <br/>
- <strong>comm Weak Scaling Graph</strong><br/>
         <img width="33%" alt="weakScalingComm" src="https://github.com/user-attachments/assets/4910dd4b-1e27-41a8-92cf-ae027b7599df">
<br/>
        The ratio used to generate the weak scaling graph is 4:2, meaning that as input size quadruples, the number of processes doubles. Ideally, we would expect a our weak scaling graph to have a positive slope of 2 to match the ratio of how we scale our array sizes to the number of processors. In actually, our weak scaling graph initially follows this pattern but soon deviates from this ideal especially as we reach higher number of processors for larger array sizes. This deviation at 32 and 64 processors could be due to load imbalances between local processor counting sorts. While we see a dip in comminucation time in previous graphs, the increase in array size dominates the benefits of parallization even at higher number of processors. 
        <br/>

    - <strong>Computation Strong Scaling Graphs</strong><br/>
      <img width="33%" alt="compStrongScale2_28" src="https://github.com/user-attachments/assets/1701f1a0-322e-4a7f-9f0f-742ab22b86f1">
<img width="33%" alt="compStrongScale2_26" src="https://github.com/user-attachments/assets/1a9667d6-082c-4e66-bdf2-35265132dc04">
<img width="33%" alt="compStrongScale2_24" src="https://github.com/user-attachments/assets/1206c733-346c-4fd8-a460-fbb2806e4ae1">
<img width="33%" alt="compStrongScale2_22" src="https://github.com/user-attachments/assets/c36da454-c92e-42ee-a68f-bffbefea8b8f">
<img width="33%" alt="compStrongScale2_20" src="https://github.com/user-attachments/assets/a4eb2e35-0340-44c0-bea3-c082d455b2ad">
<img width="33%" alt="compStrongScale2_18" src="https://github.com/user-attachments/assets/9938de2c-8a3d-4d86-8645-bdd02086ac68">
<img width="33%" alt="compStrongScale2_16" src="https://github.com/user-attachments/assets/ec806213-d2f7-4370-bdc7-e5bc4d174ff9">
      <br/>
       The computation graphs across all array sizes demonstrates a negative relationship between number of processors and total computation time. This is an expected benefit and behavior from increase parallelism. At this project's scale, breaking down the local sorts into smaller chunks and scattering them amongst more processors is beneficial. One importanat note as before is that the type should not affect the computation time too much. While with random array_types, there is less locality of seeing the same digits (102, 57, 93)  vs (999, 998, 997), a major contributing factor of why random takes longer is the fact that I had generated the array with a normal distribution from (0, INT32_MAX) instead of (0, array_size), leading to much larger amounts of digits on average which is a key factor in radix sorts overall runtime. Another important aspect to point out is the plateau in performance increase at around 64 - 128 processors. While we see the benefits of parallelism for radix sort, there are diminishing returns for higher and higher processor counts. 
        <br/>

    - <strong>Computation Speedup Graphs</strong><br/>
      <img width="33%" alt="compSpeedup2_28" src="https://github.com/user-attachments/assets/aec6bded-43fd-425d-8764-b09e3c22fc49">
<img width="33%" alt="compSpeedup2_26" src="https://github.com/user-attachments/assets/167bf322-430c-4632-afec-48549ae1d06b">
<img width="33%" alt="compSpeedup2_24" src="https://github.com/user-attachments/assets/770a7212-a37f-424c-90f6-0dd1dc5e5533">
<img width="33%" alt="compSpeedup2_22" src="https://github.com/user-attachments/assets/afc610df-9739-4541-906d-1eae19965bc8">
<img width="33%" alt="compSpeedup2_20" src="https://github.com/user-attachments/assets/6b6015f1-edc1-43bb-aa18-da6ef94059b9">
<img width="33%" alt="compSpeedup2_18" src="https://github.com/user-attachments/assets/e1bcece5-67b6-4e94-a4c9-2254ea3c49bd">
<img width="33%" alt="compSpeedup2_16" src="https://github.com/user-attachments/assets/510577d5-052c-4148-93a1-b4059d6d51eb">
  <br/>
       The speedup graphs demonstrate how computation time is largly uneffected by the array type and also works to demonstrate how effective radix sort benefits from parallelization. At most, we see a speed up of over 1000x for 1024 processors. As we increase the array sizes, we slowly approach a more ideal speed up where 64 processors trends toward a speedup of 64x and 256 processors trends toward a speed up of 256x. 
        <br/>
    - <strong>Computation Weak Scaling Graph</strong><br/>
       <img width="33%" alt="weakScalingComp" src="https://github.com/user-attachments/assets/5d27a12d-63a2-4c10-9dc3-210edb5f3437">
<br/>
       The ratio used to generate the weak scaling graph is 4:2, meaning that as input size quadruples, the number of processes doubles. Ideally, we would expect a our weak scaling graph to have a positive slope of 2 to match the ratio of how we scale our array sizes to the number of processors. In actually, our weak scaling graph seems to undershoot this expected results. This indicates a better computation performance than expected. When splitting each array into sub arrays, the computation on sorted, reverse, and one percent preturbed could all benefit more from data locality when using buckets in counting sort leading to less cache misses. 
        <br/>

    - <strong>main Strong Scaling Graphs</strong><br/>
       <img width="33%" alt="mainStrongScale2_28" src="https://github.com/user-attachments/assets/f8980601-f4de-4a12-881a-5c2e56fdba66">
<img width="33%" alt="mainStrongScale2_26" src="https://github.com/user-attachments/assets/b1546b9f-09aa-405a-84a2-d19854dc0350">
<img width="33%" alt="mainStrongScale2_24" src="https://github.com/user-attachments/assets/3de835ee-91b0-4b1b-b38e-d318f9e80c7b">
<img width="33%" alt="mainStrongScale2_22" src="https://github.com/user-attachments/assets/c3be9ab7-bc5d-431a-a25a-28eee0459aec">
<img width="33%" alt="mainStrongScale2_20" src="https://github.com/user-attachments/assets/d8c872a9-d86b-4aaa-8aff-bf6f97b180fb">
<img width="33%" alt="mainStrongScale2_18" src="https://github.com/user-attachments/assets/21dfcdf0-40fd-468f-b136-7e95291a635e">
<img width="33%" alt="mainStrongScale2_16" src="https://github.com/user-attachments/assets/994a9480-5652-4ac8-9efa-2ec2f9c68e58">
  <br/>
       ADD GRAPH ANALYSIS
        <br/>

    - <strong>main Speedup Graphs</strong><br/>
        <img width="33%" alt="mainSpeedup2_28" src="https://github.com/user-attachments/assets/e6b4a5a1-a6d3-467f-aed0-b81fc815f05c">
<img width="33%" alt="mainSpeedup2_26" src="https://github.com/user-attachments/assets/4cd68c57-06d2-487c-be9f-cedc57706b0b">
<img width="33%" alt="mainSpeedup2_24" src="https://github.com/user-attachments/assets/2361ddb0-fd36-45a0-bef4-6cebbbad2809">
<img width="33%" alt="mainSpeedup2_22" src="https://github.com/user-attachments/assets/7ed14522-a1fc-4cea-9043-717135499b2b">
<img width="33%" alt="mainSpeedup2_20" src="https://github.com/user-attachments/assets/76751945-e84a-4b24-b50b-e793d243e1a8">
<img width="33%" alt="mainSpeedup2_18" src="https://github.com/user-attachments/assets/15f749d9-2a74-4e3f-96d0-f2ee56a41427">
<img width="33%" alt="mainSpeedup2_16" src="https://github.com/user-attachments/assets/5ce4aee4-1aa1-495a-823a-2d014cb09e19">
<br/>
       ADD GRAPH ANALYSIS
        <br/>
      <strong>main Weak Scaling Graph</strong><br/>
      <img width="33%" alt="weakScalingMain" src="https://github.com/user-attachments/assets/8692d17f-ef2d-4900-b6dd-800ba085adef">
       <br/>
       ADD GRAPH ANALYSIS
        <br/>
      <strong>Cache Misses Graphs</strong><br/>
      <img width="33%" alt="l1CacheMisses" src="https://github.com/user-attachments/assets/590244f5-ea90-4366-89ef-60f0f7e0f136">
<img width="33%" alt="l2CacheMisses" src="https://github.com/user-attachments/assets/c48f2d84-efe5-40ee-b7c2-857472ba33b8">
       <br/>
       ADD GRAPH ANALYSIS
        <br/>
    - <strong>Additional Notes</strong><br/>

- Column Sort<br/>
    - <strong>comm Strong Scaling Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comm%20for%20sorting%202^16%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comm%20for%20sorting%202^18%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comm%20for%20sorting%202^20%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comm%20for%20sorting%202^22%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comm%20for%20sorting%202^24%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comm%20for%20sorting%202^26%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comm 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comm%20for%20sorting%202^28%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comm 2^28"><br/>
        For the smaller input sizes (2^16, 2^18, 2^20), there is an upwards trend for the comm time, indicating that the communication overhead is overtaking the benefit that more parallelization would give. Across all input types, the larger input sizes, and number of processes, the comm time is about the same, with a few outliers in the graphs, which can be explained with hardware variance. Addiionally, the randomness of the graphs can also be explained with hardware variance as the times are quite small. Some runs may have had nodes that are closer to each other or had nodes that ran a bit faster or slower than the others. The reason for the constant time is due to the implementation. There is a while loop that sends data whose number of iterations is on the scale of the number of processes. The less processes there are, the more data is sent per iteration, but less iterations are needed. So, the constant line seen in the graph indicates that the ratio between the data sending time and the number of iterations is constant. This also means that communication does not scale strongly with the number of proceses whatsoever. However, the algorithm does not spend a lot of time communicating for every input size, as indicated by the extremeley small times on the graphs.
        <br/>
    - <strong>comm Speedup Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comm%20for%20sorting%202^16%20elements.jpg?raw=true" width=33% alt="Speedup Graph comm 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comm%20for%20sorting%202^18%20elements.jpg?raw=true" width=33% alt="Speedup Graph comm 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comm%20for%20sorting%202^20%20elements.jpg?raw=true" width=33% alt="Speedup Graph comm 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comm%20for%20sorting%202^22%20elements.jpg?raw=true" width=33% alt="Speedup Graph comm 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comm%20for%20sorting%202^24%20elements.jpg?raw=true" width=33% alt="Speedup Graph comm 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comm%20for%20sorting%202^26%20elements.jpg?raw=true" width=33% alt="Speedup Graph comm 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comm%20for%20sorting%202^28%20elements.jpg?raw=true" width=33% alt="Speedup Graph comm 2^28"><br/>
        As seen in the Strong Scaling graphs, the smaller input sizes (2^16, 2^18, 2^20) have a decreasing Speedup as the number of processes increases. This is due to the communication overhead overtaking the benefit of using more processes, so to minimize communication time, it would be best to use 2 processes for those input types. For the rest of the graphs, the graphs appear random. This is due to the Strong Scaling graph being essentially constant, so any deviation from the average will cause a speedup or slowdown to be calculated.
        <br/>
    - <strong>comm Weak Scaling Graph</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Weak%20Scaling%20of%20comm%20for%20sorting%20elements.jpg?raw=true" width=33% alt="Weak Scaling Graph comm 2^16"><br/>
        The ratio used to produce the graph is 4:2, so every time the input size quadruples, the number of processes doubles. Thus, we would expect a slope of 2 for the graph, indicating that the work per process should double. The graph actually does follow this slope for a bit, until 512 processes, where all but Sorted input types, experience a large increase in work per process. Once again, this can be explain with hardware variability as these times are quite small. The nodes may have been further away or ran slower than the other nodes, leading to a slowdown.
        <br/>

    - <strong>comp_large Strong Scaling Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comp_large%20for%20sorting%202^16%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comp_large 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comp_large%20for%20sorting%202^18%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comp_large 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comp_large%20for%20sorting%202^20%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comp_large 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comp_large%20for%20sorting%202^22%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comp_large 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comp_large%20for%20sorting%202^24%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comp_large 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comp_large%20for%20sorting%202^26%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comp_large 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20comp_large%20for%20sorting%202^28%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph comp_large 2^28"><br/>
        All the graphs show an inversely proprotional relationship between the comp_large time and the number of processes for every input type. Interestingly, the sorted input type is faster for every input size. This is due to the algorithm using std::sort, which uses an implementation of quicksort. Since the input is already sorted, the quicksort algorithmn does not need to make as many, if any, swaps, so it runs quicker. Additionally, for the graphs that have this data point, the graphs level off at about 128 processes, indicating that there is not much more that can be parallelized and sequential runtime is now the limiting factor. These graphs also indicate that the comp_large scales strongly with the number of processes.
        <br/>
    - <strong>comp_large Speedup Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comp_large%20for%20sorting%202^16%20elements.jpg?raw=true" width=33% alt="Speedup Graph comp_large 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comp_large%20for%20sorting%202^18%20elements.jpg?raw=true" width=33% alt="Speedup Graph comp_large 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comp_large%20for%20sorting%202^20%20elements.jpg?raw=true" width=33% alt="Speedup Graph comp_large 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comp_large%20for%20sorting%202^22%20elements.jpg?raw=true" width=33% alt="Speedup Graph comp_large 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comp_large%20for%20sorting%202^24%20elements.jpg?raw=true" width=33% alt="Speedup Graph comp_large 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comp_large%20for%20sorting%202^26%20elements.jpg?raw=true" width=33% alt="Speedup Graph comp_large 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20comp_large%20for%20sorting%202^28%20elements.jpg?raw=true" width=33% alt="Speedup Graph comp_large 2^28"><br/>
        Since the Strong Scaling graphs have an inversely proportional relationship between the comp_large time and the number of processes, we would expect a proportional relationship between the comp_large Speedup and the number of processes, which is shown for every graph. Due to the extremely strong strong scaling of comp_large, using the maximum number of processes will lead to the largest speedup, which is ideal.
        <br/>
    - <strong>comp_large Weak Scaling Graph</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Weak%20Scaling%20of%20comp_large%20for%20sorting%20elements.jpg?raw=true" width=33% alt="Weak Scaling Graph comp_large 2^16"><br/>
        The ratio used to produce the graph is 4:2, so every time the input size quadruples, the number of processes doubles. Thus, we would expect a slope of 2 for the graph, indicating that the work per process should double. The graph does follow this relationship quite closely for all the input types, meaning that comp_large is also weakly scaled as the number of processes increases. Additionally, the Sorted input type is also the fastest among the input types, once again due to the underlying quicksort of std::sort not needing to make as many, if any, swaps.
        <br/>

    - <strong>main Strong Scaling Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20main%20for%20sorting%202^16%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph main 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20main%20for%20sorting%202^18%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph main 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20main%20for%20sorting%202^20%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph main 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20main%20for%20sorting%202^22%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph main 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20main%20for%20sorting%202^24%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph main 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20main%20for%20sorting%202^26%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph main 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Strong%20Scaling%20of%20main%20for%20sorting%202^28%20elements.jpg?raw=true" width=33% alt="Strong Scaling Graph main 2^28"><br/>
        For smaller input sizes (2^16, 2^18, 2^20), the graphs actually show an increase in main time as the number of processes increases. This is due to the communication overhead overtaking the benefits of using more processes, as indicated by the strong scaling graphs of comm. As the input size increases, however, the graphs start to show the expected inversely proportional relationship between the time and the number of processes. Once again, the Sorted input type is the fastest for those graphs, which was indicated in the comp_large strong scaling graphs. The scaling also begins to level off at 128 processes, which was also indicated by the comp_large strong scaling graphs. This means that for the larger input sizes, the comp_large time overtakes the comm time, so the main graphs will be closer to the comp_large graphs. The graphs for the larger input sizes also indicate that main scales strongly with the number of processes.
        <br/>
    - <strong>main Speedup Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20main%20for%20sorting%202^16%20elements.jpg?raw=true" width=33% alt="Speedup Graph main 2^16">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20main%20for%20sorting%202^18%20elements.jpg?raw=true" width=33% alt="Speedup Graph main 2^18">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20main%20for%20sorting%202^20%20elements.jpg?raw=true" width=33% alt="Speedup Graph main 2^20">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20main%20for%20sorting%202^22%20elements.jpg?raw=true" width=33% alt="Speedup Graph main 2^22">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20main%20for%20sorting%202^24%20elements.jpg?raw=true" width=33% alt="Speedup Graph main 2^24">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20main%20for%20sorting%202^26%20elements.jpg?raw=true" width=33% alt="Speedup Graph main 2^26">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Speedup%20of%20main%20for%20sorting%202^28%20elements.jpg?raw=true" width=33% alt="Speedup Graph main 2^28"><br/>
        For smaller input sizes (2^16, 2^18, 2^20), the Speedup of main goes down, which is expected since the Strong Scaling graphs indicate an increase in time as the number of processes increases. As the input size increases, the graphs show the expected proportional relationship between the Speedup and the number of processes. An interesting thing to note is that as the input size increases, the number of processes with the best speedup also increases. This can be explained by the comp_large time overtaking the comm time as the number of processes increases. Additionally, for larger input types, the Speedup begins to level off at 128 processes, in contrast to the comp_large Speedup never leveling off. This shows that the more random comm Speedup is affecting the Speedup of main, and the outliers in the comm Strong Scaling graph have a noticable effect on the runtime of main.
        <br/>
    - <strong>main Weak Scaling Graph</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_scaling/Weak%20Scaling%20of%20main%20for%20sorting%20elements.jpg?raw=true" width=33% alt="Weak Scaling Graph main 2^16"><br/>
        The ratio used to produce the graph is 4:2, so every time the input size quadruples, the number of processes doubles. Thus, we would expect a slope of 2 for the graph, indicating that the work per process should double. The graph follows this relationship until 32 processes, where it sharply increases. This shows that the algorithm does not follow the ideal weak scaling, so the amount of work per processes is not even as both the amount of processes increase and the input size increases.
        <br/>

    - <strong>Cache Misses Graphs</strong><br/>
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_cm/Min%20L1%20misses%20rank.jpg?raw=true" width="33%" alt="L1 Misses Graph">
        <img src="https://github.com/kimcchen/CSCE435_Project/blob/main/column_sort/plots_cm/Min%20L2%20misses%20rank.jpg?raw=true" width="33%" alt="L2 Misses Graph"><br/>
        The graphs show that as the input size doubles, the number of cache misses per process double. This is expected because the amount of data per process also doubles, meaning there is twice the opportunities for cache misses. They also show that doubling the number of processes will decrease the number of cache misses, which is also expected because by doubling the number of processes, the amount of data per process is halved, so there are less opportunities for cache misses to occur. Additionally, the graphs also show that the number of L1 cache misses is greater than the number of L2 caches misses, which is expected because of L1 cache is accessed first before the L2 cache, so a cache miss on L2 translates to a cache miss on L1 as well, but the L2 cache may already have the data that the L1 cache missed on.
        <br/>
    - <strong>Additional Notes</strong><br/>
        The Speedup is calculated using t_1 as the time it takes for 2 processes multipled by 2, which is why all the speedup graphs start at a 2x speedup.<br/>

        For all the Weak Scaling graphs, the ratio used is 4:2, so as the input size quadruples, the number of processes doubles. Thus, a graph with a slope of 2 is expected, where 2 indicates that the amount of work per process should double.<br/>

        Column Sort has two restrictions:
        - `numRows >= (numCols - 1)^2 * 2`
        - `numRows % numCols == 0`
        <br/>
        Thus, the number of process that can be used is bounded by the maximum number of columns in the matrix, which is why data points are missing from some of the graphs and 1024 processes is not ran at all. 1024 processes would require an input size of 2^31 (1024 x 2097152).<br/>

        Additionally, column sort runs quicker on sorted input compared to the other input types, and any other input type would take about the same time. This is due to the use of quicksort to sort the columns of the matrix in steps 1, 3, 5, and 7. Since the Sorted input type is already sorted, the quicksort algorithm will do less swaps, if any, on the data, leading to a quicker runtime.<br/>
