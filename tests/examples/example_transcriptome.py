"""
Example: Transcriptome contribution pipeline

Functionality
-------------
Compute per-connection gene contribution scores using region coordinates and a
transcriptome expression matrix. Internally, the pipeline:
1) Loads coordinates (columns: MNI_X, MNI_Y, MNI_Z)
2) Loads gene expression (rows=regions, columns=genes; optional 'label' dropped)
3) Z-scores expression per region
4) Computes CGE (correlated gene expression)
5) Fits an exponential decay model between distance and CGE
6) Builds expected matrix and derives per-connection contributions for each gene
7) Saves a CSV with 'Region_Pair' and per-gene columns

Parameters
----------
- coordinates_path (str): Path to CSV/Excel coordinates file with required columns.
- gene_expression_path (str): Path to CSV/Excel gene expression matrix.
- output_contribution_path (str): Path to save output CSV (directories auto-created).
- initial_params (tuple[float, float, float]): (A, n, B) initial guess for decay.
- save_plot (bool): Whether to save the fit figure.
- plot_path (Optional[str]): Output path for figure; if None and save_plot=True, no figure is saved.

Input file formats
------------------
- Coordinates: CSV/Excel with headers MNI_X,MNI_Y,MNI_Z.
- Gene expression: CSV/Excel; each row corresponds to a region; each column to a gene.
  A 'label' column is allowed and will be removed.

Usage
-----
Adjust the example paths below to match your environment, then run:
    python tests/examples/example_transcriptome.py

Outputs
-------
- CSV at output_contribution_path containing per-connection contributions.
- Optional figure at plot_path.
"""
import sys
import os
from pathlib import Path

# Ensure local src can be imported when not installed
ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

from BrainNetAnno.transcriptome import run_transcriptome_pipeline

# Example paths (replace with real data)
coords_csv = str(ROOT / "examples_data" / "coords.csv")
# CSV/Excel; rows=regions, columns=genes; optional 'label' column will be dropped
genes_csv = str(ROOT / "examples_data" / "gene_expression.csv")
# Output CSV path (automatically creates parent directories)
output_contrib_csv = str(ROOT / "examples_output" / "gene_contribution.csv")

# Create example dirs if needed
os.makedirs(os.path.dirname(output_contrib_csv), exist_ok=True)

# Run pipeline (adjust initial_params based on your dataset if needed)
# Returns: (A, n, B) fitted decay params and a DataFrame of contributions
params, contrib_df = run_transcriptome_pipeline(
    coordinates_path=coords_csv,
    gene_expression_path=genes_csv,
    output_contribution_path=output_contrib_csv,
    initial_params=(0.64, 90.4, -0.19),
    save_plot=False,
    plot_path=None,
)

print("Fitted parameters:", params)
print("Contribution shape:", contrib_df.shape)
