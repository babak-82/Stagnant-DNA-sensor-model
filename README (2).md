# Code for "Mass-transport limits of stagnant electrochemical DNA sensors: a reaction-diffusion analysis of response time and detection limit"

Reproduces every number, table and figure of the paper. `claims.md` lists each number quoted in the text, as produced by the last full run.

## Files (in order of importance)
| file | content |
|---|---|
| `stagnant_cell_model.py` | finite-volume solver of the reaction-diffusion model (stiff BDF), well-mixed and equilibrium references, surface-state solver, time to a coverage fraction |
| `lod_model.py` | baseline parameters (Table 1), unit conversions, calibration curves, detection limit and depletion bound |
| `reproduce_all.py` | runs all 20 stages from scratch and writes `results/*.json`, `results/claims.md` and the first-version figures |
| `edge_effect.py` | radial (edge) correction for a disk electrode, Shoup and Szabo expression |
| `make_figures.py` | publication figures 1 to 9 (new design) from `results/*.json` |
| `figstyle.py` | shared figure theme (fonts, palette, export of PNG 600 dpi, PDF, SVG) |
| `requirements.txt` | package versions used |
| `claims.md` | every number quoted in the paper, from the last full run |
| `LICENSE`, `CITATION.cff` | MIT licence and citation metadata |

## Run
```
pip install -r requirements.txt
python reproduce_all.py                    # all stages from scratch, about 9 minutes on one CPU core
python reproduce_all.py --fast             # coarse grids, pipeline test
python reproduce_all.py --stages lod,noise # selected stages (list in the docstring of reproduce_all.py)
python make_figures.py                     # figures 1 to 9 -> results/figures_v2/Figure_N.png, .pdf, .svg
python make_figures.py 3 5                 # selected figures
```
`make_figures.py` needs the `results/*.json` files, so run `reproduce_all.py` first.
Outputs: `results/*.json`, `results/claims.md`, `results/figures/*` (first-version figures), `results/figures_v2/*` (final figures).

## Notes
- Internally the code works in cm and mol; the conversions from the SI values of Table 1 are in `lod_model.derived`.
- Random seeds are fixed (Sobol sequence seed 7, bootstrap seed 11, figure jitter and schematic seed 1).
- The Sobol stage saves its evaluations to `results/sobol_evals.npy` every 50 runs and resumes from there if interrupted.
- The figures use the font Inter if it is installed and fall back to Arial, Helvetica or DejaVu Sans.
- Tested with Python 3.12.3 and the versions in `requirements.txt`.

## Licence and citation
MIT licence (see `LICENSE`). Citation metadata are in `CITATION.cff`; add the Zenodo DOI and the paper reference after publication.
