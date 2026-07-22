#############################
### 2D functions

def relocation_2D(x, y, cell_position, swap_count, rows_total, rand2):
    """
    Relocate a cell in a 2D arrangement by swapping it with a neighbor above or below.
    Args:
        x (int): Column index.
        y (int): Row index.
        cell_position (np.ndarray): 2D array of cell positions.
        swap_count (int): Number of swaps performed so far.
        rows_total (int): Total number of rows.
        rand2 (float): Random value to determine direction.
    Returns:
        tuple: Updated cell_position and swap_count.
    """
    if (rand2 > 0.5 or y == 0) and y < rows_total - 1:
        # Swap with the cell above
        cell_position[y, x], cell_position[y + 1, x] = cell_position[y + 1, x], cell_position[y, x]
    elif y > 0:
        # Swap with the cell below
        cell_position[y, x], cell_position[y - 1, x] = cell_position[y - 1, x], cell_position[y, x]
    swap_count += 1
    return cell_position, swap_count


def division_2D(x, y, cell_position, division_count, rows_total, columns_total, rand_div, rand_div_left_right):
    """
    Perform cell division in a 2D arrangement, duplicating a cell to a neighbor (left, right, or above).
    Args:
        x (int): Column index.
        y (int): Row index.
        cell_position (np.ndarray): 2D array of cell positions.
        division_count (int): Number of divisions performed so far.
        rows_total (int): Total number of rows.
        columns_total (int): Total number of columns.
        rand_div (float): Random value to determine division direction.
        rand_div_left_right (float): Random value to determine left/right division.
    Returns:
        tuple: Updated cell_position and division_count.
    """
    tmp_val = cell_position[y, x]  # Cell that duplicates
    if rand_div > 0.5:  # Lateral division
        if rand_div_left_right > 0.5:  # Right division
            if x < columns_total - 1:  # Not at the last column
                if y < rows_total - 1:
                    cell_position[y + 1:rows_total, x + 1] = cell_position[y:rows_total - 1, x + 1]
                cell_position[y, x + 1] = tmp_val  # Duplicate cell to the right
            else:  # At the last column
                if y < rows_total - 1:
                    cell_position[y + 1:rows_total, 0] = cell_position[y:rows_total - 1, 0]
                cell_position[y, 0] = tmp_val  # Duplicate cell to the first column
        else:  # Left division
            if x > 0:  # Not at the first column
                if y < rows_total - 1:
                    cell_position[y + 1:rows_total, x - 1] = cell_position[y:rows_total - 1, x - 1]
                cell_position[y, x - 1] = tmp_val  # Duplicate cell to the left
            else:  # At the first column
                if y < rows_total - 1:
                    cell_position[y + 1:rows_total, columns_total - 1] = cell_position[y:rows_total - 1, columns_total - 1]
                cell_position[y, columns_total - 1] = tmp_val  # Duplicate cell to the last column
    else:  # Upward division
        if y < rows_total - 1:
            cell_position[y + 1:rows_total, x] = cell_position[y:rows_total - 1, x]  # Shift cells down
        cell_position[y, x] = tmp_val  # Keep the original cell
    division_count += 1
    return cell_position, division_count


def division_2D_count_replacements(x, y, cell_position, division_count, rows_total, columns_total, rand_div, rand_div_left_right, mut_replaces_wt, wt_replaces_mut, mutated_cell_pool):
    """
    Perform cell division in a 2D arrangement, duplicating a cell to a neighbor (left, right, or above),
    with additional functionality for counting replacements.
    Args:
        x (int): Column index.
        y (int): Row index.
        cell_position (np.ndarray): 2D array of cell positions.
        division_count (int): Number of divisions performed so far.
        rows_total (int): Total number of rows.
        columns_total (int): Total number of columns.
        rand_div (float): Random value to determine division direction.
        rand_div_left_right (float): Random value to determine left/right division.
        mut_replaces_wt (int): Weight for mutation replacements.
        wt_replaces_mut (int): Weight for weight replacements.
        mutated_cell_pool (list): Pool of mutated cells.
    Returns:
        tuple: Updated cell_position, division_count, mut_replaces_wt, and wt_replaces_mut.
    """
    tmp_val = cell_position[y, x]  # Cell that duplicates
    if rand_div > 0.5:  # Lateral division
        if rand_div_left_right > 0.5:  # Right division
            if x < columns_total - 1: 
                target_val = cell_position[y, x+1]  # Not at the last column
                if tmp_val == 2 and target_val == 3:
                     mut_replaces_wt += 1
                if y < rows_total - 1:
                    cell_position[y + 1:rows_total, x + 1] = cell_position[y:rows_total - 1, x + 1]
                cell_position[y, x + 1] = tmp_val  # Duplicate cell to the right
            else:  # At the last column
                if y < rows_total - 1:
                    cell_position[y + 1:rows_total, 0] = cell_position[y:rows_total - 1, 0]
                cell_position[y, 0] = tmp_val  # Duplicate cell to the first column
        else:  # Left division
            if x > 0:  # Not at the first column
                target_val = cell_position[y,x -1] 
                if tmp_val == 3 and target_val == 2:
                     wt_replaces_mut += 1
                if y < rows_total - 1:
                    cell_position[y + 1:rows_total, x - 1] = cell_position[y:rows_total - 1, x - 1]
                cell_position[y, x - 1] = tmp_val  # Duplicate cell to the left
            else:  # At the first column
                if y < rows_total - 1:
                    cell_position[y + 1:rows_total, columns_total - 1] = cell_position[y:rows_total - 1, columns_total - 1]
                cell_position[y, columns_total - 1] = tmp_val  # Duplicate cell to the last column
    else:  # Upward division
        if y < rows_total - 1:
            cell_position[y + 1:rows_total, x] = cell_position[y:rows_total - 1, x]  # Shift cells down
        cell_position[y, x] = tmp_val  # Keep the original cell
    division_count += 1
    return cell_position, division_count, mut_replaces_wt, wt_replaces_mut



