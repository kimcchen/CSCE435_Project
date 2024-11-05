#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append("/scratch/group/csce435-f24/python-3.10.4/lib/python3.10/site-packages")
sys.path.append("/scratch/group/csce435-f24/thicket")
from glob import glob

import matplotlib.pyplot as plt
import pandas as pd

import thicket as th


# Read all files

# 1_trial is a name of a folder containing the cali files, you may create a folder with a different name and replace the folder name here
tk = th.Thicket.from_caliperreader(glob("/*.cali"))


# View Calltree

print(tk.tree(metric_column="Avg time/rank"))


# Group Performance data by `num_procs` and `array_size` in the Thicket metadata table.

tk.metadata_column_to_perfdata("num_procs")
tk.metadata_column_to_perfdata("array_size")

tk.dataframe = tk.dataframe.reset_index().set_index(["node", "num_procs", "array_size"]).sort_index()

tk.dataframe.head(100)


# Define common variables
processes = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
array_sizes = [65536, 262144, 1048576, 4194304, 16777216, 67108864, 268435456]


# Change font size for all plots
plt.rcParams.update({"font.size": 20})


def extract_and_plot_times(array_size, tk):
    # Dictionary of timing sections to analyze
    timing_sections = {
        'whole_program': 'main',
        'initialization': 'data_init_runtime',
        'communication': {
            'main': 'comm',
            'small': 'comm_small',
            'large': 'comm_large'
        },
        'computation': {
            'main': 'comp',
            'small': 'comp_small',
            'large': 'comp_large'
        },
        'validation': 'correctness_check'
    }
    
    def plot_timing_section(section_name, marker_name, processes):
        times = {
            'min': tk.dataframe.groupby(['num_procs', 'array_size', 'name']).min()['Min time/rank'],
            'avg': tk.dataframe.groupby(['num_procs', 'array_size', 'name']).min()['Avg time/rank'],
            'max': tk.dataframe.groupby(['num_procs', 'array_size', 'name']).min()['Max time/rank']
        }
        
        min_times = []
        avg_times = []
        max_times = []
        
        for num_proc in processes:
            try:
                min_times.append(times['min'].loc[(num_proc, array_size, marker_name)])
                avg_times.append(times['avg'].loc[(num_proc, array_size, marker_name)])
                max_times.append(times['max'].loc[(num_proc, array_size, marker_name)])
            except KeyError as e:
                print(f"Missing data for process count {num_proc} in {marker_name}")
                continue
            
        plt.figure(figsize=(12, 8))
        plt.plot(processes[:len(min_times)], min_times, marker='o', linestyle='-', color='r', label='Min Time')
        plt.plot(processes[:len(avg_times)], avg_times, marker='o', linestyle='-', color='b', label='Avg Time')
        plt.plot(processes[:len(max_times)], max_times, marker='o', linestyle='-', color='y', label='Max Time')
        
        plt.title(f'{section_name} Time vs. Processes\nArray Size: {array_size}')
        plt.xlabel('Number of Processes')
        plt.ylabel('Time (seconds)')
        plt.legend()
        plt.grid(True)
        plt.show()
        
        return min_times, avg_times, max_times

    # Plot main program timing
    main_times = plot_timing_section('Whole Program', timing_sections['whole_program'], processes)
    
    # Plot initialization timing
    init_times = plot_timing_section('Data Initialization', timing_sections['initialization'], processes)
    
    # Plot communication timings
    comm_times = {
        'main': plot_timing_section('Total Communication', timing_sections['communication']['main'], processes),
        'small': plot_timing_section('Small Communication', timing_sections['communication']['small'], processes),
        'large': plot_timing_section('Large Communication', timing_sections['communication']['large'], processes)
    }
    
    # Plot computation timings
    comp_times = {
        'main': plot_timing_section('Total Computation', timing_sections['computation']['main'], processes),
        'small': plot_timing_section('Small Computation', timing_sections['computation']['small'], processes),
        'large': plot_timing_section('Large Computation', timing_sections['computation']['large'], processes)
    }
    
    # Plot validation timing
    validation_times = plot_timing_section('Correctness Check', timing_sections['validation'], processes)
    
    # Create a comparison plot of major components
    plt.figure(figsize=(12, 8))
    
    plt.plot(processes, comm_times['main'][1], marker='o', label='Communication')
    plt.plot(processes, comp_times['main'][1], marker='s', label='Computation')
    plt.plot(processes, init_times[1], marker='^', label='Initialization')
    plt.plot(processes, validation_times[1], marker='v', label='Validation')
    
    plt.title(f'Component Time Comparison vs. Processes\nArray Size: {array_size}')
    plt.xlabel('Number of Processes')
    plt.ylabel('Average Time (seconds)')
    plt.legend()
    plt.grid(True)
    plt.show()


# Call the function to extract and plot data for a given array size
extract_and_plot_times(65536, tk)
