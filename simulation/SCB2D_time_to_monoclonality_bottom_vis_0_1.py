'''
Main simulation script for 2D cell dynamics and monoclonality analysis.
'''

# Standard library imports
import argparse
import os
from collections import Counter, defaultdict
from datetime import datetime

# Third-party imports
import numpy as np

# Local module imports
import cell_dynamics_2D
from utils import compute_relative_fraction, compute_monoclonal_fraction, transpose_dict, save_to_csv, write_nested_list_to_file, write_used_arguments_to_file


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("kr_a", type=float, help="relocation rate for cell population a")
    parser.add_argument("kd_a", type=float, help="division rate for cell population a")
    parser.add_argument("rep", type=int, help="number of times to repeat simulation")
    parser.add_argument("it", type=int, help="number of times to repeat each individual simulation run (should match number of crypts in exp.)") # crypts to simulate "in parallel" per simulation
    parser.add_argument("-kr", "--relocation", type=float, help="relocation rate for cell population b")
    parser.add_argument("-kd", "--division", type=float, help="division rate for cell population b")
    parser.add_argument("-r", "--rows", type=int, help="specify number of rows") # height of crypt
    parser.add_argument("-c", "--cols", type=int, help="specify number of columns") # width of crypt
    parser.add_argument("-stat", "--status", type=int, help="print repetition and iteration count (0=no status)")
    parser.add_argument("-id", "--file_ID", type=int, help="add ID to distinguish files to write to")
    return parser.parse_args()


def main():
    args = parse_arguments()
    cell_populations = 2 if args.relocation and args.division else 1
    rows_total = args.rows if args.rows else 10
    columns_total = args.cols if args.cols else 6
    cells_total = rows_total * columns_total
    repetitions = args.rep
    iteration_times = args.it
    kr_a = args.kr_a
    kd_a = args.kd_a
    kr_b = args.relocation if args.relocation else None
    kd_b = args.division if args.division else None
    file_id = args.file_ID
    file_name_start = f"{file_id}_m_SCB2Dmono_until_row_0"
    now = datetime.now()
    timestring = now.strftime("%Y%m%d-%H%M%S")
    # Only keep bottom file outputs
    bottom_write_to = f"{timestring}{file_name_start}_percent_fixed.csv"
    pos_and_time_write_to = f"{timestring}{file_name_start}_starting_position_time_to_monoclonality.csv"
    # Prepare header lines for CSV: arguments summary and timepoints
    arg_summary = '# Arguments: ' + ', '.join(f'{k}={v}' for k, v in vars(args).items())
    # We'll generate the timepoint header after the data is ready
    grouped_data_relative_fraction_bottom = defaultdict(list)
    store_time = [[[] for _ in range(iteration_times)] for _ in range(repetitions)]

    rng = np.random.default_rng()
    swap_count = division_count = swap_count_mutated = division_count_mutaed = 0
    time_bins = 100
    days = 1000
    defined_day = cells_total * time_bins
    for rep in range(repetitions):
        iteration_count = 0
        monoclonal_bottom = np.zeros(days)
        av_size_bottom = np.zeros((days, 2))
        for iteration in range(iteration_times):
            time_count = 0
            if args.status:
                print("rep", rep, "kr_kd", kr_a/kd_a, "iteration", iteration_count+1, "of", iteration_times)
            current_time = 0
            cell_position = np.zeros((rows_total, columns_total), dtype=int)
            # Place the initial labeled cell in row 0, random column from 0 to 4
            initial_row = 0
            initial_col = rng.integers(0, columns_total) # since always at row 0, is also initial starting position
            cell_position[initial_row, initial_col] = 1
            if cell_populations == 2:
                mutated_cell = cell_position[initial_row, initial_col]
                mutated_cell_pool = [mutated_cell]
            while current_time < days:
                time_count += 1
                for i in range(defined_day):
                    x = rng.integers(0, columns_total)
                    y = rng.integers(0, rows_total)
                    lineage_nr = cell_position[y][x]
                    if cell_populations == 2 and lineage_nr in mutated_cell_pool:
                        rand1 = rng.uniform()
                        if rand1 < kd_b / time_bins:
                            rand_div = rng.uniform()
                            rand_div_left_right = rng.random()
                            cell_position, division_count_mutaed = cell_dynamics_2D.division_2D(x, y, cell_position, division_count_mutaed, rows_total, columns_total, rand_div, rand_div_left_right)
                        rand1 = rng.uniform()
                        if rand1 < kr_b / time_bins:
                            rand2 = rng.uniform()
                            cell_position, swap_count_mutated = cell_dynamics_2D.relocation_2D(x, y, cell_position, swap_count_mutated, rows_total, rand2)
                    else:
                        rand1 = rng.uniform()
                        if rand1 < kd_a / time_bins:
                            rand_div = rng.uniform()
                            rand_div_left_right = rng.random()
                            cell_position, division_count = cell_dynamics_2D.division_2D(x, y, cell_position, division_count, rows_total, columns_total, rand_div, rand_div_left_right)
                        rand1 = rng.uniform()
                        if rand1 < kr_a / time_bins:
                            rand2 = rng.uniform()
                            cell_position, swap_count = cell_dynamics_2D.relocation_2D(x, y, cell_position, swap_count, rows_total, rand2)
                rows_0_to_3 = cell_position[0:2] #check only two most bottom
                row_0 = rows_0_to_3[0]
                is_row_0_monoclonal = np.all(row_0 == 1)
                is_row_0_all_zeros = np.all(row_0 == 0)
                no_row_has_1 = not np.any(rows_0_to_3 == 1)
                #print(f"rows 0 to 1: {rows_0_to_3}")
                #print(f"any rows 0 to 3 no 1: {no_row_has_1}")

            
                if is_row_0_monoclonal:
                    # Row 0 has become monoclonal by labeled cell (=1)
                    av_size_bottom[current_time:, 1] += 1
                    monoclonal_bottom[current_time:] += 1
                
                    store_time[rep][iteration_count].extend([current_time, initial_col])
                    break
                
                if is_row_0_all_zeros and no_row_has_1:
                    current_lineage = 999
                    store_time[rep][iteration_count].extend([current_time, current_lineage])
                    # Row 0 has become monoclonal by unlabeled cell (=0)
                    break

                elif not is_row_0_monoclonal:
                    # Row 0 is not yet monoclonal labeled cell (=1), including all-zeros case
                    av_size_bottom[current_time][1] += 1

                current_time += 1
            iteration_count += 1
        cumulative_monoclonal_bottom = monoclonal_bottom
        if np.sum(monoclonal_bottom != 0):
            cumulative_monoclonal_bottom = np.array(cumulative_monoclonal_bottom)
            rel_frac_bottom = compute_relative_fraction(cumulative_monoclonal_bottom, av_size_bottom)
            grouped_data_relative_fraction_bottom[rep] = rel_frac_bottom.tolist()
    # Now, after all loops, generate and write the header only once
    grouped_data_relative_fraction_transpose_bottom = transpose_dict(grouped_data_relative_fraction_bottom)
    if grouped_data_relative_fraction_transpose_bottom:
        # Only pass the argument summary as header_lines
        header_lines = [arg_summary]
    else:
        header_lines = [arg_summary]
    save_to_csv(grouped_data_relative_fraction_transpose_bottom, bottom_write_to, header_lines=header_lines)
    if store_time:
        # Only pass the argument summary as header_lines
        header_lines = [arg_summary]
    else:
        header_lines = [arg_summary]
    write_nested_list_to_file(store_time, pos_and_time_write_to, header_lines=header_lines)


if __name__ == "__main__":
    main()








