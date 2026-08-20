# Immune Control of Oncogene Selection in the Colon.

This repository contains the Python scripts used for the analyses presented in:

> **Mirza H. et al. (2026). _Immune Control of Oncogene Selection in the Colon._**

The simulations are based on the **Stochastic Conveyor Belt (SCB)** framework introduced by Corominas-Murtra *et al.* [1] and extended to two-dimensional intestinal crypt dynamics by Azkanaz *et al.* [2].

## Overview

This repository contains the scripts used to simulate stochastic cell dynamics in the colonic crypt and to analyze monoclonality.

## Repository structure

### 2D crypt simulations 

**`SCB2D_time_to_monoclonality_bottom_vis_0_1.py`**

Main simulation script for 2D cell dynamics and monoclonality analysis.

**`SCB2D_replacement_probabilities.py`**

Script to count cell replacements of labeled cells versus unlabeled cell.

These script imports:

- **`cell_dynamics_2D.py`** – Functions implementing the two-dimensional Stochastic Conveyor Belt (SCB) model.
- **`utils.py`** – Utility functions for simulation, data processing, and file input/output.

## Running the simulations

```
usage:
SCB2D_time_to_monoclonality_bottom_vis_0_1.py
[-h]
[-kr RELOCATION]
[-kd DIVISION]
[-r ROWS]
[-c COLS]
[-stat STATUS]
[-id FILE_ID]
kr_a kd_a rep it
```

### Positional arguments

| Argument | Description |
|----------|-------------|
| `kr_a` | Relocation rate for cell population A |
| `kd_a` | Division rate for cell population A |
| `rep` | Number of independent simulation repeats |
| `it` | Number of iterations per simulation |

### Optional arguments

| Argument | Description |
|----------|-------------|
| `-kr`, `--relocation` | Relocation rate for cell population B |
| `-kd`, `--division` | Division rate for cell population B |
| `-r`, `--rows` | Number of crypt rows |
| `-c`, `--cols` | Number of crypt columns |
| `-stat`, `--status` | Print repetition and iteration progress (`0` = off) |
| `-id`, `--file_ID` | Identifier appended to output files |

### Example

```bash
python SCB2D_time_to_monoclonality_bottom_vis_0_1.py \
0.0625 0.25 100 100 \
-r 20 \
-c 5
```

## Visualizing simulation results 

The folder scb2d_analysis stores all scripts and notebooks to create the plots presented in the manuscript:

python run_notebooks.py plot_simulation_results_publication.ipynb
  - P(Monoclonal | Visible)
  - Extinction Probability
  - Fixation Probability

python run_notebooks.py plot_influence_of_kr_publication.ipynb
for changing mutants' kr
  - P(Monoclonal | Visible)
  - Extinction Probability
  - Fixation Probability

python run_notebooks.py plot_analytical_solution_publication.ipynb
  - Analytical solution of 1D ring model compared to SCB in 2D

The folder contains also precomupted simulation results in data.
The experimental data is not available.

To plot all figures at once.
python run_notebooks.py

the figures will be stored in \scb2d_analysis\figures


## Analytical model

The analytical solution is implemented in:

- *(add filenames here)*

## Requirements

This code requires Python 3.x together with the packages listed in `requirements.txt`.

## License

Specify the license under which this code is distributed (e.g. MIT License).

## References

[1] Corominas-Murtra, B. et al. *Stem cell lineage survival as a noisy competition for niche access.* Proc. Natl. Acad. Sci. USA **117**, 16969–16975 (2020).

[2] Azkanaz, M. et al. *Retrograde movements determine effective stem cell numbers in the intestine.* Nature **607**, 548–554 (2022).

[3] Mirza, H. et al. *Immune Control of Oncogene Selection in the Colon.* Nature (2026).
