"""
Reproduce the figures of the manuscript without Jupyter.

The script executes the code cells of the figure notebooks in order and writes
every panel into the `figures` directory next to this file. It is the
non-interactive counterpart of opening the notebooks and running them from top
to bottom; the notebooks remain the place where the figures are defined, so the
two cannot disagree.

Usage
-----
    python run_notebooks.py                 all notebooks
    python run_notebooks.py <notebook> ...  the given notebooks only

Only matplotlib and the packages the analysis itself needs are required; Jupyter
is not. With Jupyter installed, the equivalent command is

    jupyter nbconvert --to notebook --execute --inplace <notebook>
"""
import json
import os
import sys
import warnings

import matplotlib

# The figures are written to file, so no interactive backend is needed.
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (after the backend is chosen)

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

# The notebooks locate the package by walking up from the working directory,
# which fails when the script is started from somewhere else. Putting the
# directory that contains the package on the import path here makes the script
# independent of where it is run from.
_PACKAGE_PARENT = os.path.dirname(PACKAGE_DIR)
if _PACKAGE_PARENT not in sys.path:
    sys.path.insert(0, _PACKAGE_PARENT)

# Executed in this order; the notebooks are independent of each other.
NOTEBOOKS = (
    "plot_simulation_results_publication.ipynb",
    "plot_influence_of_kr_publication.ipynb",
    "plot_analytical_solution_publication.ipynb",
)

FIGURE_DIR = os.path.join(PACKAGE_DIR, "figures")


def run_notebook(path):
    """
    Execute the code cells of one notebook in a namespace of their own.

    Markdown cells are ignored. An exception in a cell aborts the notebook, as
    it would when running the cells by hand.
    """
    with open(path, encoding="utf-8") as handle:
        notebook = json.load(handle)

    namespace = {"__name__": "__main__"}
    file_name = os.path.basename(path)

    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if not source.strip():
            continue
        exec(compile(source, f"{file_name}:cell{index}", "exec"), namespace)


def main(arguments):
    # Read by the notebooks; without it they compute the figures but write none.
    os.environ["SCB2D_SAVE_FIGURES"] = "1"

    # `plt.show()` has nothing to show on a non-interactive backend.
    warnings.filterwarnings("ignore", message="FigureCanvasAgg is non-interactive")

    notebooks = arguments or NOTEBOOKS
    for name in notebooks:
        path = name if os.path.isabs(name) else os.path.join(PACKAGE_DIR, name)
        print(f"--- {os.path.basename(path)}")
        run_notebook(path)
        plt.close("all")

    written = sorted(os.listdir(FIGURE_DIR)) if os.path.isdir(FIGURE_DIR) else []
    print(f"\n{len(written)} figures in {FIGURE_DIR}")
    for name in written:
        print(f"  {name}")


if __name__ == "__main__":
    main(sys.argv[1:])
