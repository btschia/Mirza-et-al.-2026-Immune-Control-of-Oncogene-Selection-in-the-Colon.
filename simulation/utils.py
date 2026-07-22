"""
Utility functions for cell dynamics simulation, data processing, and file I/O.
"""
import numpy as np
import csv
from collections import defaultdict
from datetime import datetime

# --- Data processing helpers ---
def compute_relative_fraction(cumulative_monoclonal, av_size):
    denom = av_size[:, 1]
    valid_mask = (denom != 0) & ~np.isnan(denom)
    return np.where(valid_mask, cumulative_monoclonal / denom, 0)

def compute_monoclonal_fraction(cumulative_monoclonal):
    cum = np.array(cumulative_monoclonal)
    last_value = cum[-1]
    if last_value == 0:
        return np.zeros_like(cum)
    valid_mask = (cum != 0.0) & (~np.isnan(cum)) & (cum <= last_value)
    return np.where(valid_mask, cum / last_value, 0)

def transpose_dict(data):
    transposed_data = defaultdict(list)
    for key, values in data.items():
        for index, value in enumerate(values):
            transposed_data[index].append(value)
    return dict(transposed_data)

# --- File I/O helpers ---
def generate_file_name(file_name_start):
    now = datetime.now()
    timestring = now.strftime("%Y%m%d-%H%M%S")
    file_time_to_monoclonality = file_name_start + '_time_to_monoclonality_' + timestring + '.csv'
    return file_time_to_monoclonality

def write_nested_list_to_file(file_name, data):
    flattened = []
    for sublist in data:
        row = []
        for lst in sublist:
            row.append(','.join(map(str, lst)))
        flattened.append(row)
    with open(file_name, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(flattened)

def write_used_arguments_to_file(file_name, given_arguments):
    ns_dict = vars(given_arguments)
    with open(file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(ns_dict.keys())
        writer.writerow(ns_dict.values())

def save_to_csv(transposed_data, file_name, header_lines=None):
    if not any(transposed_data.values()):
        print(f"No entries found in transposed_data. CSV file '{file_name}' will not be written.")
        return
    with open(file_name, mode="w", newline="") as file:
        # Only write argument summary/comment lines
        if header_lines:
            for line in header_lines:
                file.write(line + "\n")
        writer = csv.writer(file)
        # Always write the CSV header line once
        writer.writerow(["Timepoint"] + [f"Rep_{i}" for i in range(len(next(iter(transposed_data.values()))) )])
        for index, values in transposed_data.items():
            writer.writerow([index] + values)

def write_nested_list_to_file(data, file_name, header_lines=None):
    if not data:
        print(f"No entries found in store_time. CSV file '{file_name}' will not be written.")
        return
    flattened = []
    for sublist in data:
        row = []
        for lst in sublist:
            row.append(','.join(map(str, lst)))
        flattened.append(row)

    # Write to CSV file
    with open(file_name, 'a', newline='') as file:
        # Only write argument summary/comment lines
        if header_lines:
            for line in header_lines:
                file.write(line + "\n")
        writer = csv.writer(file)
        writer.writerows(flattened)

def is_two_left_of_three_row0(grid):
    for x in range(4):  # avoid index out of bounds
        if grid[0, x] == 2 and grid[0, x + 1] == 3:
            return True
    return False
