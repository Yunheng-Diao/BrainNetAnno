"""
Example: Mitochondrial contribution pipeline

Functionality
-------------
Compute per-connection mitochondrial contribution scores analogous to neurotransmitter.

Parameters, formats, usage, and outputs
---------------------------------------
- coordinates_path: CSV/Excel with MNI_X,MNI_Y,MNI_Z
- mitochondrial_expression_path: CSV/Excel; rows=regions, cols=features
- output_contribution_path: Output CSV with 'Region_Pair' and per-feature columns
- initial_params, save_plot, plot_path: same semantics as other contribution pipelines

Usage:
    python tests/examples/example_mitochondrial.py
"""
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

from BrainNetAnno.mitochondrial import run_mitochondrial_pipeline

coords_csv = str(ROOT / "examples_data" / "coords.csv")
mito_xlsx = str(ROOT / "examples_data" / "mitochondrial.xlsx")
mito_contrib_csv = str(ROOT / "examples_output" / "mitochondrial_contribution.csv")

os.makedirs(os.path.dirname(mito_contrib_csv), exist_ok=True)

params, contrib_df = run_mitochondrial_pipeline(
    coordinates_path=coords_csv,
    mitochondrial_expression_path=mito_xlsx,
    output_contribution_path=mito_contrib_csv,
    initial_params=(1.0, 50.0, 0.0),
    save_plot=False,
    plot_path=None,
)

print("Fitted parameters:", params)
print("Contribution shape:", contrib_df.shape)
