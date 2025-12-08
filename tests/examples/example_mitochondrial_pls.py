"""
Example: Mitochondrial PLS analysis

Functionality
-------------
Link FC weights/deviations with per-connection mitochondrial contributions using
PLS regression, CV component selection, and permutation tests for weights.

Parameters, formats, usage, outputs
-----------------------------------
- fc_matrix_path: square CSV (no header), upper-triangle non-zero entries
- nt_contrib_path: CSV with 'Region_Pair' and per-feature columns
- top_k: select top-|FC| connections
- plots and weights CSV paths as in neurotransmitter PLS example

Usage:
    python tests/examples/example_mitochondrial_pls.py
"""
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

from BrainNetAnno.mitochondrial_pls_cge import run_mitochondrial_pls_pipeline

fc_csv = str(ROOT / "examples_data" / "fc_weights.csv")
mito_contrib_csv = str(ROOT / "examples_output" / "mitochondrial_contribution.csv")
weights_csv = str(ROOT / "examples_output" / "mitochondrial_weights.csv")
mse_plot = str(ROOT / "examples_output" / "mitochondrial_pls_mse.tif")
var_plot = str(ROOT / "examples_output" / "mitochondrial_explained_variance.tif")

os.makedirs(os.path.dirname(weights_csv), exist_ok=True)

best_n, df = run_mitochondrial_pls_pipeline(
    fc_matrix_path=fc_csv,
    nt_contrib_path=mito_contrib_csv,
    top_k=500,
    mse_plot_path=mse_plot,
    explained_variance_plot_path=var_plot,
    output_weights_path=weights_csv,
    max_components=6,
    cv_splits=5,
    random_state=42,
    n_permutations=1000,
)

print("Best components:", best_n)
print("Weights saved to:", weights_csv)
