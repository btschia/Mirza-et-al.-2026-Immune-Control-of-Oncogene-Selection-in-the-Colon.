"""
Loader for the experimental (mouse) data the simulations are compared against.

The Excel sheet has a fixed layout: the conditions are arranged in horizontal
blocks of `block_width` columns separated by `gap` empty columns, and within a
block the rows "Avg" and "SD" hold the measured percentage of fixed crypts at
the sampling days.
"""

import os

import pandas as pd

# Sampling days of the experiment, in days after induction.
DEFAULT_SAMPLING_DAYS = [6, 9, 15, 22]

# Rows of the sheet on which a condition name is expected.
DEFAULT_NAME_ROWS = [4, 14]

# Row labels that are not condition names.
NON_CONDITION_LABELS = {"nan", "avg", "sd", "se", ""}


def load_experimental_data_if_available(file_path, **loader_arguments):
    """
    Load the experimental data if the sheet is present on this machine.

    The measured data are not part of the repository. If the file is absent, an
    empty dict is returned, so that the notebooks can be executed on the
    simulation output alone; the figures are then drawn without the measured
    data points.

    Parameters
    ----------
    file_path : str
        Path to the Excel file.
    loader_arguments
        Passed on to `load_experimental_data_from_excel`.

    Returns
    -------
    dict
        The loaded data, or an empty dict if the file does not exist.
    """
    if not os.path.exists(file_path):
        return {}

    return load_experimental_data_from_excel(file_path, **loader_arguments)


def load_experimental_data_from_excel(file_path,
                                      name_rows=DEFAULT_NAME_ROWS,
                                      block_width=5,
                                      gap=1,
                                      default_days=DEFAULT_SAMPLING_DAYS):
    """
    Load the experimental data from the Excel sheet with the fixed block layout.

    Parameters
    ----------
    file_path : str
        Path to the Excel file.
    name_rows : list of int
        Row indices on which condition names are located.
    block_width : int
        Number of columns per data block.
    gap : int
        Number of empty columns between two blocks.
    default_days : list of int
        Sampling days assigned to the values of a condition. A condition with a
        single value is assumed to be measured on day 9.

    Returns
    -------
    dict
        "block<N>_<condition name>" -> {"Avg": [...], "SD": [...], "Days": [...]}
    """
    sheet = pd.read_excel(file_path, header=None)

    experimental_data = {}
    num_columns = sheet.shape[1]
    block_step = block_width + gap

    for block_start in range(0, num_columns, block_step):
        block_end = block_start + block_width

        # Incomplete or entirely empty blocks are skipped.
        if block_end > num_columns or sheet.iloc[:, block_start:block_end].isna().all().all():
            continue

        for name_row in name_rows:
            condition_name = str(sheet.iloc[name_row, block_start]).strip()
            if condition_name.lower() in NON_CONDITION_LABELS:
                continue

            # The Avg and SD rows are located in the 14 rows below the condition name.
            block = sheet.iloc[name_row + 1:name_row + 15, block_start:block_end]
            first_column = block.iloc[:, 0].astype(str)
            average_row = block[first_column.str.fullmatch("Avg", case=False, na=False)]
            std_row = block[first_column.str.fullmatch("SD", case=False, na=False)]

            averages = (pd.to_numeric(average_row.iloc[0, 1:], errors="coerce")
                        .dropna().tolist()) if not average_row.empty else []
            standard_deviations = (pd.to_numeric(std_row.iloc[0, 1:], errors="coerce")
                                   .dropna().tolist()) if not std_row.empty else []

            days = [9] if len(averages) == 1 else default_days[:len(averages)]
            block_id = f"block{block_start // block_step + 1}"

            experimental_data[f"{block_id}_{condition_name}"] = {
                "Avg": averages,
                "SD": standard_deviations,
                "Days": days,
            }

    return experimental_data
