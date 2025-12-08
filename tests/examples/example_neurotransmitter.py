"""
Example: Neurotransmitter contribution pipeline

Functionality
-------------
Compute per-connection neurotransmitter contribution scores from region coordinates
and a transmitter/receptor expression matrix.

Pipeline steps
--------------
1) Load coordinates (MNI_X, MNI_Y, MNI_Z)
2) Load neurotransmitter expression (rows=regions, cols=features)
3) Z-score per region, build precision/covariance as needed
4) Fit distance-decay model and compute expected matrix
5) Derive per-connection contributions and save CSV

Parameters
----------
- coordinates_path (str): CSV/Excel coordinates file path.
- neurotransmitter_expression_path (str): CSV/Excel transmitter/receptor matrix path.
- output_contribution_path (str): Output CSV path; directories are created automatically.
- initial_params (tuple[float, float, float]): (A, n, B) initial guess for decay model.
- save_plot (bool): Whether to save the fit figure.
- plot_path (Optional[str]): Path for the figure; if None, no figure is saved.

Input file formats
------------------
- Coordinates: CSV/Excel with headers MNI_X,MNI_Y,MNI_Z.
- Neurotransmitter expression: CSV/Excel; rows=regions, cols=features; optional 'label' dropped.

Usage
-----
Adjust paths and run:
    python tests/examples/example_neurotransmitter.py

Outputs
-------
- CSV at output_contribution_path with 'Region_Pair' and per-feature columns.
- Optional figure if plot_path is provided.
"""
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

from BrainNetAnno.neurotransmitter import run_neurotransmitter_pipeline

coords_csv = str(ROOT / "examples_data" / "coords.csv")
nt_xlsx = str(ROOT / "examples_data" / "neurotransmitter.xlsx")
nt_contrib_csv = str(ROOT / "examples_output" / "neurotransmitter_contribution.csv")

os.makedirs(os.path.dirname(nt_contrib_csv), exist_ok=True)

params, contrib_df = run_neurotransmitter_pipeline(
    coordinates_path=coords_csv,
    neurotransmitter_expression_path=nt_xlsx,
    output_contribution_path=nt_contrib_csv,
    initial_params=(1.0, 50.0, 0.0),
    save_plot=False,
    plot_path=None,
)

print("Fitted parameters:", params)
print("Contribution shape:", contrib_df.shape)
