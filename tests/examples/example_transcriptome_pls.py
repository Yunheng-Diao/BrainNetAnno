"""
Example: Transcriptome PLS-CGE analysis

Functionality
-------------
Link FC deviations (t-values) with per-connection gene contributions using
Partial Least Squares (PLS) regression, cross-validated component selection,
and permutation-based significance.

Pipeline steps
--------------
1) Read FC deviation matrix (square CSV, no header); use upper-triangle non-zero entries
2) Read per-connection gene contributions (CSV with 'Region_Pair' and gene columns)
3) Align rows by 'Region_Pair' and construct X (genes) and Y (FC values)
4) Select PLS components via CV (MSE) or use provided n_components
5) Fit PLS, compute explained variance (X) and permutation p-values
6) Save significant genes and optional plots (MSE curve, explained variance)

Parameters
----------
- t_values_file_path (str): Path to FC deviations CSV (square, no header).
- gene_expression_file_path (str): Path to per-connection gene contribution CSV.
- output_best_genes_path (str): Path to save significant genes CSV.
- mse_plot_path (Optional[str]): Path to save CV MSE plot; paired *_data.csv is written.
- explained_variance_plot_path (Optional[str]): Path to save explained-variance plot; paired *_data.csv is written.
- max_components (int): Upper bound of components to search (capped internally by data).
- cv_splits (int): KFold splits for CV.
- random_state (int): Random seed for reproducibility.
- n_permutations (int): Permutations for significance testing.
- n_components (Optional[int]): If set, bypass CV and use fixed components.

Input file formats
------------------
- FC deviations: square CSV without header; only upper-triangle non-zero entries are used.
- Gene contributions: CSV with 'Region_Pair' column and per-gene columns.

Usage
-----
Adjust paths and run:
    python tests/examples/example_transcriptome_pls.py

Outputs
-------
- Significant genes CSV at output_best_genes_path
- MSE and explained-variance plots (if paths provided) with paired *_data.csv
"""
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

from BrainNetAnno.transcriptome_pls_cge import run_transcriptome_pls_pipeline

# Example paths
fc_csv = str(ROOT / "examples_data" / "fc_deviation.csv")
contrib_csv = str(ROOT / "examples_output" / "gene_contribution.csv")
best_genes_csv = str(ROOT / "examples_output" / "best_genes.csv")
mse_plot = str(ROOT / "examples_output" / "transcriptome_pls_mse.tif")
var_plot = str(ROOT / "examples_output" / "transcriptome_explained_variance.tif")

os.makedirs(os.path.dirname(best_genes_csv), exist_ok=True)

run_transcriptome_pls_pipeline(
    t_values_file_path=fc_csv,
    gene_expression_file_path=contrib_csv,
    output_best_genes_path=best_genes_csv,
    mse_plot_path=mse_plot,
    explained_variance_plot_path=var_plot,
    max_components=10,
    cv_splits=5,
    random_state=42,
    n_permutations=1000,
)

print("Saved:", best_genes_csv)
