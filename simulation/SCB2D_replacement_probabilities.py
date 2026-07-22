import numpy as np
import argparse
from collections import Counter
from datetime import datetime
import cell_dynamics_2D
import utils

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("kr_a", type=float, help="relocation rate for cell population a")
    parser.add_argument("kd_a", type=float, help="division rate for cell population a")
    parser.add_argument("rep", type=int, help="number of times to repeat simulation")
    parser.add_argument("it", type=int, help="number of times to repeat each indiviual simulation run")
    parser.add_argument("-kr", "--relocation", type=float, help="relocation rate for cell population b")
    parser.add_argument("-kd", "--division", type=float, help="division rate for cell population b")
    parser.add_argument("-r", "--rows", type=int, help="specify number of rows")
    parser.add_argument("-c", "--cols", type=int, help="specify number of columns")
    parser.add_argument("-stat", "--status", type=int, help="print repetition and iteration count (0=no status)")
    parser.add_argument("-id", "--file_ID", type=int, help="add ID to distinguish files to write to")
    parser.add_argument("-dt", "--delta_t", type=float, help="delta t, 10^x")
    parser.add_argument("--min_replacements", type=int, default=None, help="Minimum number of valid replacements required per simulation (optional)")
    args = parser.parse_args()

    cell_populations = 2 if args.relocation and args.division else 1
    rows_total = args.rows if args.rows else 10
    columns_total = args.cols if args.cols else 6
    cells_total = rows_total * columns_total
    repetitions = args.rep
    iteration_times = args.it

    # Sanity check: if min_replacements is set and greater than iteration_times, warn and exit
    if args.min_replacements is not None and args.min_replacements > iteration_times:
        print("min_replacements cannot be greater than the number of iterations per repetition (it). Reduce --min_replacements or increase --it.")
        return

    u = 10
    u **= -args.delta_t if args.delta_t else -2

    kr_a = args.kr_a
    kd_a = args.kd_a
    kr_b = args.relocation if args.relocation else None
    kd_b = args.division if args.division else None

    file_id = args.file_ID
    if cell_populations == 2:
        file_name_start = f"{file_id}_m_SCB2Dmono_until_row_4_{u}"
    else:
        file_name_start = f"{file_id}_SCB_2Dmono_until_row_4_{u}"
    now = datetime.now()
    timestring = now.strftime("%Y%m%d-%H%M%S")
    file_time_to_monoclonality = f"{timestring}_{file_name_start}_replacement_counts.csv"
    arg_summary = '# Arguments: ' + ', '.join(f'{k}={v}' for k, v in vars(args).items())
    rng = np.random.default_rng()
    count_cell_action = np.zeros((rows_total, columns_total))
    store_replacements_wt_mut = [[[] for _ in range(iteration_times)] for _ in range(repetitions)]
    time_bins = 100
    current_time = 0

    for rep in range(repetitions):
        iteration_count = 0
        # If --min_replacements is set, only simulations with at least that many total replacements (mut_replaces_wt + wt_replaces_mut)
        # are counted as valid and stored. The script will keep running new simulations until the required number of valid ones
        # (iteration_times) is reached for each repetition. If --min_replacements is not set, all simulations are stored as usual.
        if args.min_replacements is not None:
            # Only count valid iterations, ignore 'it' argument for total attempts
            while iteration_count < iteration_times:
                if args.status:
                    print(f"rep {rep} kr_kd {kr_a/kd_a} valid iteration {iteration_count+1} of {iteration_times}")
                cell_position = np.zeros((rows_total, columns_total), dtype=int)
                division_count = 0
                division_count_mutated = 0
                mut_replaces_wt = 0
                wt_replaces_mut = 0
                swap_count = 0
                swap_count_mutated = 0
                cell_position[0, 2] = 2
                mutated_cell = cell_position[0, 2]
                mutated_cell_pool = [mutated_cell]
                cell_position[0, 3] = 3
                while utils.is_two_left_of_three_row0(cell_position):
                    x = rng.integers(0, columns_total-1, endpoint=True)
                    y = rng.integers(0, rows_total-1, endpoint=True)
                    lineage_nr = cell_position[y][x]
                    count_cell_action[y][x] += 1
                    if cell_populations == 2 and lineage_nr in mutated_cell_pool:
                        rand1 = rng.uniform()
                        if rand1 < kd_b / time_bins:
                            rand_div = rng.uniform()
                            rand_div_left_right = rng.random()
                            cell_position, division_count_mutated, mut_replaces_wt, wt_replaces_mut = cell_dynamics_2D.division_2D_count_replacements(
                                x, y, cell_position, division_count_mutated, rows_total, columns_total, rand_div, rand_div_left_right, mut_replaces_wt,
                                wt_replaces_mut, mutated_cell_pool)
                        rand1 = rng.uniform()
                        if rand1 < kr_b / time_bins:
                            rand2 = rng.uniform()
                            cell_position, swap_count_mutated = cell_dynamics_2D.relocation_2D(
                                x, y, cell_position, swap_count_mutated, rows_total, rand2)
                    else:
                        rand1 = rng.uniform()
                        if rand1 < kd_a / time_bins:
                            rand_div = rng.uniform()
                            rand_div_left_right = rng.random()
                            cell_position, division_count, mut_replaces_wt, wt_replaces_mut = cell_dynamics_2D.division_2D_count_replacements(
                                x, y, cell_position, division_count, rows_total, columns_total, rand_div, rand_div_left_right, mut_replaces_wt,
                                wt_replaces_mut, mutated_cell_pool)
                        rand1 = rng.uniform()
                        if rand1 < kr_a / time_bins:
                            rand2 = rng.uniform()
                            cell_position, swap_count = cell_dynamics_2D.relocation_2D(
                                x, y, cell_position, swap_count, rows_total, rand2)
                # Optional: enforce minimum replacements if set
                if (mut_replaces_wt + wt_replaces_mut) >= args.min_replacements:
                    store_replacements_wt_mut[rep][iteration_count].extend([mut_replaces_wt, wt_replaces_mut])
                    iteration_count += 1
                    current_time += 1
        else:
            for iteration in range(iteration_times):
                if args.status:
                    print(f"rep {rep} kr_kd {kr_a/kd_a} iteration {iteration+1} of {iteration_times}")
                cell_position = np.zeros((rows_total, columns_total), dtype=int)
                division_count = 0
                division_count_mutated = 0
                mut_replaces_wt = 0
                wt_replaces_mut = 0
                swap_count = 0
                swap_count_mutated = 0
                cell_position[0, 2] = 2
                mutated_cell = cell_position[0, 2]
                mutated_cell_pool = [mutated_cell]
                cell_position[0, 3] = 3
                while utils.is_two_left_of_three_row0(cell_position):
                    x = rng.integers(0, columns_total-1, endpoint=True)
                    y = rng.integers(0, rows_total-1, endpoint=True)
                    lineage_nr = cell_position[y][x]
                    count_cell_action[y][x] += 1
                    if cell_populations == 2 and lineage_nr in mutated_cell_pool:
                        rand1 = rng.uniform()
                        if rand1 < kd_b / time_bins:
                            rand_div = rng.uniform()
                            rand_div_left_right = rng.random()
                            cell_position, division_count_mutated, mut_replaces_wt, wt_replaces_mut = cell_dynamics_2D.division_2D_count_replacements(
                                x, y, cell_position, division_count_mutated, rows_total, columns_total, rand_div, rand_div_left_right, mut_replaces_wt,
                                wt_replaces_mut, mutated_cell_pool)
                        rand1 = rng.uniform()
                        if rand1 < kr_b / time_bins:
                            rand2 = rng.uniform()
                            cell_position, swap_count_mutated = cell_dynamics_2D.relocation_2D(
                                x, y, cell_position, swap_count_mutated, rows_total, rand2)
                    else:
                        rand1 = rng.uniform()
                        if rand1 < kd_a / time_bins:
                            rand_div = rng.uniform()
                            rand_div_left_right = rng.random()
                            cell_position, division_count, mut_replaces_wt, wt_replaces_mut = cell_dynamics_2D.division_2D_count_replacements(
                                x, y, cell_position, division_count, rows_total, columns_total, rand_div, rand_div_left_right, mut_replaces_wt,
                                wt_replaces_mut, mutated_cell_pool)
                        rand1 = rng.uniform()
                        if rand1 < kr_a / time_bins:
                            rand2 = rng.uniform()
                            cell_position, swap_count = cell_dynamics_2D.relocation_2D(
                                x, y, cell_position, swap_count, rows_total, rand2)
                store_replacements_wt_mut[rep][iteration].extend([mut_replaces_wt, wt_replaces_mut])
                current_time += 1
    cell_position = np.zeros((rows_total, columns_total), dtype=int)
    iteration_count += 1
    current_time = 0
    if store_replacements_wt_mut:
        # Only pass the argument summary as header_lines
        header_lines = [arg_summary]
    else:
        header_lines = [arg_summary]
    utils.write_nested_list_to_file(store_replacements_wt_mut, file_time_to_monoclonality, header_lines=header_lines )

if __name__ == "__main__":
    main()








