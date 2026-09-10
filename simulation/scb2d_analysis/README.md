# scb2d_analysis

Analysis and figure code for the SCB 2D crypt simulation, as used in the
manuscript. The package reads the raw output of the simulation, computes the
reported statistics and assembles the figure panels. It is self-contained: the
modules, the three figure notebooks and the simulation output they read are all
located in this directory, so that the figures can be reproduced without any
configuration.

## Contents

| File | Contents |
| --- | --- |
| `data/` | the simulation output the notebooks read, about 150 MB |
| `datasets.py` | input paths and the assignment of simulations to conditions; the only module containing paths |
| `io_simulation.py` | readers for the three raw simulation output formats |
| `monoclonality.py` | statistics for P(Monoclonal \| Visible) |
| `fixation.py` | fixation and extinction probability curves over simulation time |
| `replacement.py` | replacement probability of a labeled stem cell |
| `experimental_data.py` | loader for the measured mouse data the simulations are compared against |
| `analytical.py` | analytical clone size distribution of the one-dimensional stem cell ring |
| `plot_style.py` | palette, line styles and the shared matplotlib configuration |
| `plots.py` | plotting functions, each drawing into the current axes |
| `plot_simulation_results_publication.ipynb` | the main figures of the manuscript |
| `plot_influence_of_kr_publication.ipynb` | sensitivity of the results to the replacement rate `kr` |
| `plot_analytical_solution_publication.ipynb` | simulation and measurement against the analytical solution |
| `run_notebooks.py` | executes the notebooks without Jupyter and exports the panels |

## Requirements

Python 3.12 with `numpy`, `scipy`, `pandas`, `matplotlib`, `cycler` and
`openpyxl` (required by `pandas` to read the Excel sheet). The notebooks were
last executed with Python 3.12.4, numpy 2.3.4, scipy 1.15.3, pandas 2.2.3,
matplotlib 3.9.2, cycler 0.12.1 and openpyxl 3.1.5.

The figures use the Arial font family. On a system without Arial, matplotlib
substitutes its default font; the panels remain correct but differ typographically
from the published version.

## Input data

The simulation runs the figures are based on are shipped in `data/`, so nothing
has to be configured after downloading the repository:

```
data/
    wt/, braf/                 runs of the main figures
    replacement/<run>/         replacement probability output
    kr/<condition>/<variant>/  runs of the kr sensitivity analysis
    analytical/<run>/          runs compared against the analytical solution
```

Every run folder except those below `replacement` is divided into the
subfolders `pos_and_time` and `percent_fixed`, one per simulation output format.
The folder names are shortened compared to the names the simulation writes, which
keeps the paths below the Windows limit of 260 characters; which run a folder
holds is documented where it is configured in `datasets.py`.

To analyse **different runs**, point `SIMULATION_ROOT` at another directory of
the same layout, either at the top of `datasets.py` or through the environment
variable `SCB2D_SIMULATION_ROOT`. All remaining paths are derived from it, so no
other module has to be adjusted.

`datasets.missing_paths()` returns every configured simulation folder that is not
present. It is evaluated in the first cell of the notebooks, so that an
unreachable input is reported before the analysis is started.

The **measured data** are not part of the repository. They are looked for in
`DOWNLOAD_ROOT` (environment variable `SCB2D_DOWNLOAD_ROOT`), by default as
`Final_Fixed Values(1).xlsx`.

The measured data are optional. If the Excel sheet is not available, the
notebooks report

```
no experimental data available; the figures are drawn without the measured data points
```

and every figure is produced from the simulation output alone. Whether the sheet
is present can be queried with `datasets.experimental_data_available()`; the
loader `experimental_data.load_experimental_data_if_available` returns an empty
dict in that case.

Three output formats of the simulation are read (see `io_simulation.py`):

- *percent fixed*: one line per simulation day, holding P(Monoclonal | Visible)
  of every repetition.
- *pos and time*: for every crypt, the day it became monoclonal and the founder
  lineage that took it over.
- *replacement probability*: for every crypt, how often a labeled and an
  unlabeled cell replaced the resident stem cell.

The parameters appearing in the folder names are `p1` and `p2` (division
probabilities), `kr` (replacement rate) and `kd` (death rate of the labeled,
mutant cell).

## Running the notebooks

Start Jupyter or an editor with Jupyter support in this directory, or in the
directory above it, and execute the notebooks from top to bottom. The first cell
of each notebook adds the directory containing the package to `sys.path`, so
that `import scb2d_analysis` resolves in both cases.

`plot_simulation_results_publication.ipynb` produces P(Monoclonal | Visible)
over time, the fixation and the extinction probability, and the distribution of
the replacement probability. `plot_influence_of_kr_publication.ipynb` produces
the same three probability panels for the runs with a modified `kr`.
`plot_analytical_solution_publication.ipynb` compares the simulated and the
measured curves with the analytical solution of the stem cell ring, for equally
likely replacement directions and for a bias between them.

Reading the simulation folders is the time-consuming step; it is done once per
notebook and the result reused by the subsequent figures.

## Exporting the panels

Every figure cell ends with

```python
if SAVE_FIGURES:
    save_figure("<panel name>", FIGURE_DIR)
```

`SAVE_FIGURES` is defined in the first cell of each notebook and is off by
default, so running a notebook interactively writes no files. Set it to `True`
there to export while working, or reproduce every panel at once with

```powershell
python run_notebooks.py
```

which executes the three notebooks in order and writes the twelve panels as PDF
into `figures/`. The script needs matplotlib only; with Jupyter installed the
equivalent for a single notebook is

```powershell
jupyter nbconvert --to notebook --execute --inplace <notebook>
```

The panels are defined in the notebooks alone, so the script and an interactive
run cannot produce different figures. To change the file format, pass
`formats=("pdf", "svg")` to `save_figure`, or edit `FIGURE_FORMATS` in
`plot_style.py`.

## Applying the sensitivity analysis to another parameter

`plot_influence_of_kr_publication.ipynb` is not specific to `kr`. The runs it
analyses are described by a nested dictionary

```python
{condition: {variant: run folder}}
```

which is built with `datasets.build_sensitivity_runs` and assigned to
`SENSITIVITY_RUNS` in the first cell of the notebook. The variant named
`datasets.REFERENCE_VARIANT` is the run the remaining variants are compared
against. Every figure, color and legend entry is derived from that dictionary,
so a different parameter requires no further modification. `KR_SENSITIVITY` in
`datasets.py` is the definition used in the manuscript and can serve as a
template.

## Conventions

- Every plotting function draws into the *current* matplotlib axes, so that the
  notebook retains control over figure size, layout and export.
  `plot_style.setup_plot_style` is called after the figure has been created and
  before plotting into it.
- The simulation counts days from 0 and the experiment from 1; the modules shift
  the simulated days accordingly, so that both are reported on the same axis.
- The independent unit is the simulation run, not the individual crypt. Crypts
  are therefore averaged within a run first, and the reported statistics are
  computed over the run means.
