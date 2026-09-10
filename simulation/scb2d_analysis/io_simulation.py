"""
Readers for the raw output files of the SCB 2D crypt simulation.

Every reader takes a directory and returns plain Python containers, so that the
statistics modules remain independent of the file format.

Three output formats are handled:

`read_visible_monoclonal_probabilities_from_directory`
    "percent fixed" files: one line per simulation day, holding
    P(Monoclonal | Visible) of every repetition.

`read_fixation_events_from_directory`
    "pos and time" files: for every crypt, the day it became monoclonal and the
    founder lineage that took it over.

`read_replacement_counts_from_directory`
    "replacement probability" files: for every crypt, how often a labeled and
    an unlabeled cell replaced the resident stem cell.
"""

import csv
import os
from collections import defaultdict

# Lineage id written by the simulation when a crypt was colonized by an
# *unlabeled* cell. Labeled cells carry ids 0-5, i.e. one of
# the bottom stem cell positions.
UNLABELED_LINEAGE_ID = 999

# Status string written for crypts that never became monoclonal within the
# simulated period.
NOT_MONOCLONAL_MARKER = "not monoclonal"

# File extensions produced by the simulation.
SIMULATION_FILE_EXTENSIONS = (".csv", ".txt")

# The simulation writes its parameters into the first line of every file, as a
# comment of the form
#     # Arguments: kr_a=0.0625, kd_a=0.25, rep=100, it=100, ...
COMMENT_PREFIX = "#"


def _data_lines(file):
    """
    Yield the lines of a simulation file, without the header comment.
    """
    return (line for line in file if not line.startswith(COMMENT_PREFIX))


def _parse_numeric_token(token):
    """Convert one CSV token to int or float; other tokens are returned unchanged."""
    try:
        return float(token) if "." in token else int(token)
    except ValueError:
        return token  # for example a status string such as "not monoclonal"


def _parse_integer_token(token):
    """Convert one CSV token to int; tokens that fail are returned as stripped text."""
    try:
        return int(token)
    except ValueError:
        return token.strip()


def _simulation_files(directory):
    """Yield the full path of every simulation output file in `directory`."""
    for file_name in sorted(os.listdir(directory)):
        if file_name.endswith(SIMULATION_FILE_EXTENSIONS):
            yield os.path.join(directory, file_name)


def read_visible_monoclonal_probabilities_from_directory(directory):
    """
    Read P(Monoclonal | Visible) per simulation day from a "percent fixed" folder.

    File layout: a comment line with the simulation parameters, a column header
    ("Timepoint, Rep_0, Rep_1, ..."), then one line per simulation day,
        day, value_repetition_1, value_repetition_2, ...
    Lines that cannot be parsed, the header line among them, is skipped.

    Returns
    -------
    dict
        simulation day -> list of P(Monoclonal | Visible), pooled over all
        repetitions found in the folder.
    """
    probabilities_by_day = defaultdict(list)

    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)

        # Only files are processed; subdirectories are ignored.
        if not os.path.isfile(file_path):
            continue

        with open(file_path, "r", encoding="utf-8") as file:
            for line in file:
                fields = line.strip().split(",")
                try:
                    day = int(fields[0])
                    probabilities = [float(value) for value in fields[1:]]
                except (ValueError, IndexError):
                    # Rows with invalid or missing values are ignored.
                    continue
                probabilities_by_day[day].extend(probabilities)

    return probabilities_by_day


def read_fixation_events_from_directory(directory, verbose=False):
    """
    Read the monoclonal conversion events from a "pos and time" folder.

    File layout:
      - one header line with the simulation parameters, which is skipped
      - one row per simulation run; each run simulates 100 crypts in parallel
      - each cell of a row holds one comma-separated "day, lineage_id" pair:
        the day the crypt became monoclonal, and the lineage that colonized
          (UNLABELED_LINEAGE_ID = unlabeled founder cell)

    Returns
    -------
    list
        simulation_runs[run][crypt] -> [day, lineage_id]
    """
    simulation_runs = []

    for file_path in _simulation_files(directory):
        if verbose:
            print(f"reading {os.path.basename(file_path)}")

        try:
            with open(file_path, "r", newline="") as file:
                reader = csv.reader(_data_lines(file))

                for row in reader:
                    # Skip rows where all cells are empty (trailing blank lines).
                    if all(cell.strip() == "" for cell in row):
                        continue

                    day_lineage_pairs = []
                    for cell in row:
                        if cell.strip() == "":
                            continue

                        # Pairs may be quoted ("3, 12") or bare (3, 12).
                        cell_content = cell.strip().strip('"')
                        day_lineage_pairs.append([
                            _parse_numeric_token(token)
                            for token in cell_content.split(",") if token.strip()
                        ])

                    simulation_runs.append(day_lineage_pairs)

        except FileNotFoundError:
            print(f"File not found: {file_path}")
        except Exception as error:
            print(f"An error occurred while reading '{file_path}': {error}")

    return simulation_runs


def read_replacement_counts_from_directory(directory):
    """
    Read the stem cell replacement counts from a "replacement probability" folder.

    File layout:
      - one header line with the simulation parameters, which is skipped
      - one row per simulation run
      - each cell of a row holds one record per crypt, either
          "labeled_replacements, unlabeled_replacements"                  (crypt became monoclonal)
        or
          "not monoclonal, labeled_replacements, unlabeled_replacements"  (crypt became not monoclonal)

    Returns
    -------
    list
        simulation_runs[run][crypt] -> [labeled, unlabeled]
        or [NOT_MONOCLONAL_MARKER, labeled, unlabeled]
    """
    simulation_runs = []

    for file_path in _simulation_files(directory):
        try:
            with open(file_path, "r", newline="") as file:
                reader = csv.reader(_data_lines(file))

                for row in reader:
                    # Skip rows where all cells are empty (trailing blank lines).
                    if all(cell.strip() == "" for cell in row):
                        continue

                    replacement_records = []
                    for cell in row:
                        if cell.strip() == "":
                            continue

                        cell_content = cell.strip().strip('"')  # remove outer quotes
                        tokens = [token.strip()
                                  for token in cell_content.split(",") if token.strip()]

                        if NOT_MONOCLONAL_MARKER in tokens:
                            # Keep the marker in front and convert the counts behind it.
                            record = [NOT_MONOCLONAL_MARKER] + [
                                _parse_integer_token(token)
                                for token in tokens if token != NOT_MONOCLONAL_MARKER
                            ]
                        else:
                            record = [_parse_integer_token(token) for token in tokens]

                        replacement_records.append(record)

                    simulation_runs.append(replacement_records)

        except FileNotFoundError:
            print(f"File not found: {file_path}")
        except Exception as error:
            print(f"Error reading '{file_path}': {error}")

    return simulation_runs
