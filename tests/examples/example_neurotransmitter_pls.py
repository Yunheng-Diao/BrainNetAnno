"""
Example: Neurotransmitter PLS analysis

Functionality
-------------
Link FC weights/deviations with per-connection neurotransmitter contributions
using PLS, CV-based component selection, and permutation tests of weights.

Pipeline steps
--------------
1) Read FC matrix (square CSV, no header) and keep upper-triangle non-zero entries
2) Read per-connection contributions (CSV with 'Region_Pair')
3) Align data by 'Region_Pair' and construct X (features) and Y (FC values)
4) Select components via CV (MSE); fit PLS
5) Compute explained variance and permutation p-values of weights
6) Save weights/p-values and optional plots

Parameters
----------
- fc_matrix_path (str): Path to square FC CSV (no header).
- nt_contrib_path (str): Path to per-connection contribution CSV.
- top_k (int): Select top-k connections by absolute FC to reduce noise.
- mse_plot_path (Optional[str]): Path to save CV MSE plot; paired *_data.csv is written.
- explained_variance_plot_path (Optional[str]): Path to save explained-variance plot; paired *_data.csv is written.
- output_weights_path (Optional[str]): Path to save weights+p-values CSV.
- max_components (int), cv_splits (int), random_state (int), n_permutations (int): Analysis settings.

Input file formats
------------------
- FC: square CSV without header; non-zero upper-triangle entries used.
- Contributions: CSV with 'Region_Pair' and per-feature columns.

Usage
-----
Adjust paths and run:
    python tests/examples/example_neurotransmitter_pls.py

Outputs
-------
- CSV at output_weights_path with weights and permutation p-values.
- MSE and explained-variance plots if paths provided.
"""
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

from BrainNetAnno.neurotransmitter_pls_cge import run_neurotransmitter_pls_pipeline

fc_csv = str(ROOT / "examples_data" / "fc_weights.csv")
nt_contrib_csv = str(ROOT / "examples_output" / "neurotransmitter_contribution.csv")
weights_csv = str(ROOT / "examples_output" / "neurotransmitter_weights.csv")
mse_plot = str(ROOT / "examples_output" / "neurotransmitter_pls_mse.tif")
var_plot = str(ROOT / "examples_output" / "neurotransmitter_explained_variance.tif")

os.makedirs(os.path.dirname(weights_csv), exist_ok=True)

best_n, df = run_neurotransmitter_pls_pipeline(
    fc_matrix_path=fc_csv,
    nt_contrib_path=nt_contrib_csv,
    top_k=350,
    mse_plot_path=mse_plot,
    explained_variance_plot_path=var_plot,
    output_weights_path=weights_csv,
    max_components=10,
    cv_splits=5,
    random_state=42,
    n_permutations=1000,
)

print("Best components:", best_n)
print("Weights saved to:", weights_csv)
