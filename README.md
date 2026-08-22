# Immune Control of Oncogene Selection in the Colon.

This repository contains the Python scripts used for the analyses presented in:

> **Mirza H. et al. (2026). _Immune Control of Oncogene Selection in the Colon._ Research Square preprint, [doi:10.21203/rs.3.rs-8673399/v1](https://doi.org/10.21203/rs.3.rs-8673399/v1)**

The simulations are based on the **Stochastic Conveyor Belt (SCB)** framework introduced by Corominas-Murtra *et al.* [1], applied to two-dimensional intestinal crypt dynamics by Azkanaz *et al.* [2] and further extended in this work to incorporate biased competition.

## Overview

This repository provides the scripts used to simulate stochastic cell dynamics in the colonic crypt and to analyze monoclonality, based on the SCB framework [1,2].

---

## Repository Structure

- **2D crypt simulations**
  - `SCB2D_time_to_monoclonality_bottom_vis_0_1.py`: Main simulation script for 2D cell dynamics and monoclonality analysis.
  - `SCB2D_replacement_probabilities.py`: Script to count cell replacements of labeled vs. unlabeled cells.
- **Modules**
  - `cell_dynamics_2D.py`: Functions implementing the two-dimensional SCB model.
  - `utils.py`: Utility functions for simulation, data processing, and file I/O.
- **Analysis**
  - `scb2d_analysis/`: Scripts and notebooks for data processing and visualization.
  - `data/`: Precomputed simulation results.
  - `figures/`: Output figures.

---


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
| `it` | Number of iterations per simulation (= nr. of crypts labelled in parallel)|

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
python SCB2D_time_to_monoclonality_bottom_vis_0_1.py 0.0625 0.25 100 100 -r 20 -c 5 -stat 1
```

## Visualizing simulation results 
All code for data processing and visualization was refactored and optimized using the AI-assited tool (Claude Code).

Absolutely! Here is the **suggested text in Markdown**, ready to be pasted into your README:

---

## Running Simulations

To run the main simulation script, use the following command:

```bash
python SCB2D_time_to_monoclonality_bottom_vis_0_1.py [options] kr_a kd_a rep it
```

### Positional Arguments

| Argument | Description |
|----------|-------------|
| `kr_a`   | Relocation rate for cell population A |
| `kd_a`   | Division rate for cell population A   |
| `rep`    | Number of independent simulation repeats |
| `it`     | Number of iterations per simulation (number of crypts labelled in parallel) |

### Optional Arguments

| Argument | Description |
|----------|-------------|
| `-kr`, `--relocation` | Relocation rate for cell population B |
| `-kd`, `--division`   | Division rate for cell population B   |
| `-r`, `--rows`        | Number of crypt rows                  |
| `-c`, `--cols`        | Number of crypt columns               |
| `-stat`, `--status`   | Print repetition and iteration progress (`0` = off) |
| `-id`, `--file_ID`    | Identifier appended to output files   |
| `-h`, `--help`        | Show help message and exit            |

### Example Usage

Run a simulation with specific parameters:

```bash
python SCB2D_time_to_monoclonality_bottom_vis_0_1.py 0.0625 0.25 100 100 -r 20 -c 5 -stat 1
```

This command runs 100 repeats of 100 iterations each, with 20 rows and 5 columns in the crypt, and prints progress.

---

## Visualizing Simulation Results

All scripts and notebooks for data processing and visualization are located in the `scb2d_analysis` folder.

### To generate all figures at once:

```bash
python run_notebooks.py
```

### To run individual notebooks:

- **Simulation results:**  
  ```bash
  python run_notebooks.py plot_simulation_results_publication.ipynb
  ```
  - Generates: P(Monoclonal | Visible), Extinction Probability, Fixation Probability

- **Influence of mutant relocation rate:**  
  ```bash
  python run_notebooks.py plot_influence_of_kr_publication.ipynb
  ```
  - Generates: P(Monoclonal | Visible), Extinction Probability, Fixation Probability for varying mutant `kr`

- **Analytical solution comparison:**  
  ```bash
  python run_notebooks.py plot_analytical_solution_publication.ipynb
  ```
  - Compares the analytical solution of the 1D ring model to the 2D SCB simulation

**Note:**  
- Precomputed simulation results are available in the `data/` folder.
- Experimental data is **not** included.
- All generated figures will be saved in `scb2d_analysis/figures/`.

---

## Analytical Solution

The analytical solution is implemented in:

- `scb2d_analysis/plot_analytical_solution_publication.ipynb`  
  (Run with `python run_notebooks.py plot_analytical_solution_publication.ipynb`)

This notebook compares the analytical solution of the 1D ring model to the 2D SCB simulation results.

---

<!-- 

## Requirements

This code requires Python 3.x together with the packages listed in `requirements.txt`.

## License

Specify the license under which this code is distributed (e.g. MIT License).
-->
## References

[1] Corominas-Murtra, B. et al. *Stem cell lineage survival as a noisy competition for niche access.* Proc. Natl. Acad. Sci. USA **117**, 16969–16975 (2020).

[2] Azkanaz, M. et al. *Retrograde movements determine effective stem cell numbers in the intestine.* Nature **607**, 548–554 (2022).

[3] Mirza, H. et al. *Immune Control of Oncogene Selection in the Colon.* Research Square preprint, https://doi.org/10.21203/rs.3.rs-8673399/v1 (2026).
