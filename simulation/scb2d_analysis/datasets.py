"""
Location of the input data and assignment of the simulations to the conditions.

The simulation output is shipped with the package, in the `data` directory next
to this module, so that a fresh download reproduces the figures without any
configuration. To analyse a different set of runs, point `SIMULATION_ROOT` at
another directory of the same layout, either below or through the environment
variable `SCB2D_SIMULATION_ROOT`; no other module requires modification.

Layout of the data directory:

    wt/, braf/                 runs of the main figures
    replacement/<run>/         replacement probability output
    kr/<condition>/<variant>/  runs of the kr sensitivity analysis
    analytical/<run>/          runs compared against the analytical solution

Every run folder except those below `replacement` holds two subfolders,
`pos_and_time` and `percent_fixed`, one per simulation output format; see
`pos_and_time_folder` and `percent_fixed_folder`. The folder names are shortened
compared to the names the simulation writes, which keeps the paths below the
Windows limit of 260 characters. The parameters of a run are documented where it
is configured below.

Parameter naming used in the notes:
    p1, p2  division probabilities
    kr      replacement rate
    kd      death rate of the labeled (mutant) cell
"""

import os

# --- Roots -------------------------------------------------------------------
# The simulation output shipped with the package.
PACKAGED_DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# Overridden by the environment variables of the same name, if these are set.
SIMULATION_ROOT = os.environ.get("SCB2D_SIMULATION_ROOT", PACKAGED_DATA_ROOT)

# The measured data are not part of the package; this is where they are looked
# for. See `experimental_data_available`.
DOWNLOAD_ROOT = os.environ.get(
    "SCB2D_DOWNLOAD_ROOT",
    os.path.expanduser(os.path.join("~", "Downloads")),
)

_REPLACEMENT_ROOT = os.path.join(SIMULATION_ROOT, "replacement")
_KR_ROOT = os.path.join(SIMULATION_ROOT, "kr")
_ANALYTICAL_ROOT = os.path.join(SIMULATION_ROOT, "analytical")


# --- Experimental data -------------------------------------------------------
EXPERIMENTAL_DATA_FILE = os.path.join(DOWNLOAD_ROOT, "Final_Fixed Values(1).xlsx")

# Sampling days of the experiment.
EXPERIMENT_DAYS = [6, 9, 15, 22]


# --- Output formats of a run -------------------------------------------------
# Subfolders a run folder is divided into, one per simulation output format.
POS_AND_TIME_SUBFOLDER = "pos_and_time"
PERCENT_FIXED_SUBFOLDER = "percent_fixed"


def pos_and_time_folder(run_folder):
    """Return the subfolder of `run_folder` with the "pos and time" output."""
    return os.path.join(run_folder, POS_AND_TIME_SUBFOLDER)


def percent_fixed_folder(run_folder):
    """Return the subfolder of `run_folder` with the "percent fixed" output."""
    return os.path.join(run_folder, PERCENT_FIXED_SUBFOLDER)


# --- Runs of the main figures ------------------------------------------------
# Condition -> run folder. WT is the wild type of the proximal colon, BRAF the
# a mutant (kd = 0.5).
MAIN_RUNS = {
    "WT": os.path.join(SIMULATION_ROOT, "wt"),
    "BRAF": os.path.join(SIMULATION_ROOT, "braf"),
}

# Condition -> folder with the "percent fixed" output, for P(Monoclonal | Visible).
PERCENT_FIXED_FOLDERS = {
    condition: percent_fixed_folder(run_folder)
    for condition, run_folder in MAIN_RUNS.items()
}

# Condition -> folder with the "pos and time" output, for fixation and extinction.
FIXATION_FOLDERS = {
    condition: pos_and_time_folder(run_folder)
    for condition, run_folder in MAIN_RUNS.items()
}


# --- Comparison with the analytical solution ---------------------------------
# Runs that are compared against the analytical clone size distribution of the 
# one-dimensional stem cell ring (see `analytical.py`). The name of a run is 
# its figure label.
ANALYTICAL_COMPARISON_RUNS = {
    # kr = 0, kd = 0.21: without replacement between the rows the simulation
    # reduces to the one-dimensional ring of the analytical solution.
    "WT, kr = 0, kd = 0.21": os.path.join(_ANALYTICAL_ROOT, "wt_kr0"),
    # p1 = 0.0525, p2 = 0.21
    "WT, p1 = 0.0525, p2 = 0.21": os.path.join(_ANALYTICAL_ROOT, "wt"),
    # p1 = 0.0625, p2 = 0.25, kr = 0.0625, kd = 0.25, labeling row 0 / visibility
    # row 0 and 1. The labeled cells carry no advantage, kd equals p2, so this is
    # a wild type as well and is the comparison for the mutant runs below.
    "WT, LI 0-1, 0.0625 / 0.25": os.path.join(_ANALYTICAL_ROOT, "li_0_1"),
    # as above, with a mutant kd of 0.5, 0.67 and 0.83
    "LI 0-1, 0.0625 / 0.25, mutant kd = 0.5": os.path.join(_ANALYTICAL_ROOT, "li_0_1_kd0_5"),
    "LI 0-1, 0.0625 / 0.25, mutant kd = 0.67": os.path.join(_ANALYTICAL_ROOT, "li_0_1_kd0_67"),
    "LI 0-1, 0.0625 / 0.25, mutant kd = 0.83": os.path.join(_ANALYTICAL_ROOT, "li_0_1_kd0_83"),
}


def analytical_comparison_missing_paths():
    """Return every folder of `ANALYTICAL_COMPARISON_RUNS` that is missing."""
    return [folder
            for run_folder in ANALYTICAL_COMPARISON_RUNS.values()
            for folder in (pos_and_time_folder(run_folder),
                           percent_fixed_folder(run_folder))
            if not os.path.exists(folder)]


# --- Sensitivity runs --------------------------------------------------------
# A sensitivity analysis compares runs that differ in a single simulation
# parameter (here kr), for one or more conditions. It is represented by a nested dict
#
#     {condition: {variant: run folder}}
#
# where the run folder holds the two output subfolders described above.
#
# The structure is not specific to a particular parameter. A different one is
# analysed by building a dict of the same shape with `build_sensitivity_runs`
# and assigning it in `plot_influence_of_kr_publication.ipynb`, from which the
# figures are derived.

# Name of the run the other variants of a condition are compared against.
REFERENCE_VARIANT = "reference"


def build_sensitivity_runs(root, folder_names_by_condition):
    """
    Build a sensitivity analysis from folder *names* relative to `root`.

    Parameters
    ----------
    root : str
        Directory the runs live in.
    folder_names_by_condition : dict
        {condition: {variant: folder name}}. The condition and variant keys are
        free text and serve as figure labels.

    Returns
    -------
    dict
        {condition: {variant: absolute run folder}}
    """
    return {
        condition: {variant: os.path.join(root, folder_name)
                    for variant, folder_name in folder_names.items()}
        for condition, folder_names in folder_names_by_condition.items()
    }


def sensitivity_missing_paths(sensitivity_runs, include_percent_fixed=True):
    """
    Return every folder of `sensitivity_runs` that is missing on this machine.

    Serves the same purpose as `missing_paths`, for runs that are not part of
    the default figure set.
    """
    missing = []
    for variants in sensitivity_runs.values():
        for run_folder in variants.values():
            expected = [pos_and_time_folder(run_folder)]
            if include_percent_fixed:
                expected.append(percent_fixed_folder(run_folder))
            missing += [folder for folder in expected if not os.path.exists(folder)]
    return missing


# The analysis shown in the manuscript: the replacement rate of the mutant cell,
# kr, doubled and halved around its fitted value of 0.0625, with all remaining
# parameters held fixed. All runs use labeling row 0 / visibility row 0 and 1.
KR_SENSITIVITY = build_sensitivity_runs(_KR_ROOT, {
    # kd = 0.83, p1 = 0.0625, p2 = 0.25
    "BRAF; MHCII fl/fl": {
        REFERENCE_VARIANT: os.path.join("mhcii_fl_fl", "reference"),
        "kr x 2": os.path.join("mhcii_fl_fl", "kr_x2"),
        "kr / 2": os.path.join("mhcii_fl_fl", "kr_div2"),
    },
    # kd = 0.5, p1 = 0.0625, p2 = 0.25
    "BRAF; MHCII fl/+": {
        REFERENCE_VARIANT: os.path.join("mhcii_fl_het", "reference"),
        "kr x 2": os.path.join("mhcii_fl_het", "kr_x2"),
        "kr / 2": os.path.join("mhcii_fl_het", "kr_div2"),
    },
})


# --- Replacement probability -------------------------------------------------
# Condition -> {"path": folder with the replacement output, "parameters": note}.
# The order of the entries is the order in which the conditions are reported.
#
# The three mutant runs are the counterparts of the runs of the analytical
# comparison above: same p1, p2 and kr, mutant kd of 0.5, 0.67 and 0.83.
REPLACEMENT_DATASETS = {
    "WT": {
        "path": os.path.join(_REPLACEMENT_ROOT, "wt"),
        "parameters": "kr = 0.0625, p2 = 0.25",
    },
    "LI 0-1, 0.0625 / 0.25, mutant kd = 0.5": {
        "path": os.path.join(_REPLACEMENT_ROOT, "kd_0_5"),
        "parameters": "kd = 0.5, labeling row 0 / visibility row 0 and 1",
    },
    "LI 0-1, 0.0625 / 0.25, mutant kd = 0.67": {
        "path": os.path.join(_REPLACEMENT_ROOT, "kd_0_67"),
        "parameters": "kd = 0.67, labeling row 0 / visibility row 0 and 1",
    },
    "LI 0-1, 0.0625 / 0.25, mutant kd = 0.83": {
        "path": os.path.join(_REPLACEMENT_ROOT, "kd_0_83"),
        "parameters": "kd = 0.83, labeling row 0 / visibility row 0 and 1",
    },
}

# Conditions included in the figures, in reporting order.
REPORTED_REPLACEMENT_CONDITIONS = list(REPLACEMENT_DATASETS)


def experimental_data_available():
    """
    Return True if the Excel sheet with the measured data is present.

    The measured data are optional: without them the notebooks draw the
    simulated curves alone, which is why the sheet is not part of
    `missing_paths`.
    """
    return os.path.exists(EXPERIMENTAL_DATA_FILE)


def missing_paths():
    """
    Return every configured simulation folder that does not exist on this
    machine.

    Intended to be evaluated in the first cell of a notebook, so that an
    unreachable input is reported before the analysis is started rather than in
    the middle of it. The optional experimental data are reported separately by
    `experimental_data_available`.
    """
    configured = list(PERCENT_FIXED_FOLDERS.values())
    configured += list(FIXATION_FOLDERS.values())
    configured += [dataset["path"] for dataset in REPLACEMENT_DATASETS.values()]
    return [path for path in configured if not os.path.exists(path)]
