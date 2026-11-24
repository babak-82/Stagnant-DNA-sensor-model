Stagnant Electrochemical DNA Sensor Model

Python code for: Mass-transport limits of stagnant electrochemical DNA sensors: a reaction-diffusion analysis of response time and detection limit

This repository contains the Python implementation of the reaction-diffusion model and computational analyses associated with the manuscript.

The code is intended to support computational reproducibility by providing the model implementation, numerical analyses, supplementary checks, and figure-generation scripts.

Repository Contents

File| Description
"stagnant_cell_model.py"| Finite-volume solver for the reaction-diffusion model, including stiff BDF integration, well-mixed and equilibrium references, surface-state calculations, and time-to-coverage analysis.
"lod_model.py"| Baseline parameters, unit conversions, calibration curves, detection-limit calculations, and depletion bounds.
"reproduce_all.py"| Executes the computational workflow and generates numerical outputs, claims, and figures.
"edge_effect.py"| Radial edge-effect correction for disk electrodes using the Shoup and Szabo expression.
"make_figures.py"| Generates publication figures from the computational outputs.
"figstyle.py"| Shared figure styling, typography, and export settings for PNG, PDF, and SVG formats.
"extra_checks.py"| Additional numerical and physical consistency checks, including solver tolerance, finite-domain Laplace testing, beta-table analysis, association-rate scenarios, optimum probe density, and thermal estimates.
"requirements.txt"| Python package dependencies.
"claims.md"| Numerical claims and values associated with the manuscript.
"LICENSE"| MIT License.
"CITATION.cff"| Citation metadata for this research software.

Installation

Python and the required dependencies must be installed before running the scripts.

Clone or download this repository, navigate to its directory, and install the dependencies:

pip install -r requirements.txt

Reproducibility Workflow

Run the scripts in the following order.

1. Full computational workflow

python reproduce_all.py

This script executes the computational stages and generates the numerical outputs required by the subsequent analyses and figure-generation workflow.

2. Additional numerical and physical checks

python extra_checks.py

This script performs supplementary checks of the model assumptions, numerical behavior, and selected physical estimates.

3. Figure generation

python make_figures.py

The figure-generation script uses the numerical outputs produced by "reproduce_all.py".

To generate selected figures, the script also supports figure-specific arguments:

python make_figures.py 3 5

Optional execution modes

For a faster, coarse-grid pipeline test:

python reproduce_all.py --fast

To execute selected computational stages:

python reproduce_all.py --stages lod,noise

Available stages are documented in the "reproduce_all.py" source code.

Computational Outputs

The computational workflow generates numerical outputs and intermediate files under the "results/" directory.

These include:

- Numerical JSON outputs.
- Updated computational claims.
- First-version figures.
- Final figure exports.

The figure-generation workflow produces PNG, PDF, and SVG files.

The exact output structure depends on the execution stage and script configuration.

Numerical and Physical Checks

The supplementary analysis in "extra_checks.py" covers:

- Solver-tolerance sensitivity.
- Finite-domain Laplace-equation testing.
- Beta-parameter table analysis.
- Association-rate ("k_on") scenarios.
- Optimal probe-density analysis.
- Thermal-estimation checks.

These checks are intended to provide additional scrutiny of the numerical implementation and selected physical assumptions.

Implementation Notes

- The implementation internally uses centimeters (cm) and moles (mol) for its calculations.
- Unit conversions from the SI parameters are handled through the model implementation.
- The computational workflow uses fixed random seeds for reproducibility, including the Sobol sequence, bootstrap analysis, and figure-related randomization.
- The Sobol evaluation stage supports saving intermediate evaluations and resuming interrupted computations.
- Figure typography uses Inter when available, with alternative system fonts as fallbacks.

Reproducibility and Claims

The "claims.md" file documents numerical values associated with the manuscript.

For reproducibility, the claims should be updated from the outputs of the latest successful full computational run.

Users should independently execute the scripts and compare the resulting numerical values and figures with the corresponding manuscript statements.

License

This software is distributed under the MIT License. See the "LICENSE" file for the complete license terms.

Citation

If you use this software in academic work, please cite the repository and its archived version.

Citation metadata are provided in "CITATION.cff".

A version-specific Zenodo DOI will be added after the first archived release.

Author

Babak Sharifi Taskuh

Islamic Azad University, Sari Branch

ORCID: "0009-0002-5449-1696" (https://orcid.org/0009-0002-5449-1696)

Repository

GitHub: https://github.com/babak-82/Stagnant-DNA-sensor-model
