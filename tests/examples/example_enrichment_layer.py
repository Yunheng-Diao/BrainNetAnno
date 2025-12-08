"""
Example: Cortical layer enrichment analysis

Functionality
-------------
Permutation-based enrichment test to assess whether target genes show elevated
layer-specific t-statistics compared to random gene sets.

Inputs
------
- layer_marker_path: Excel with sheet 'Table S4B' containing 'gene' and t_stat_LayerX columns
- target_genes_path: CSV with column 'Gene Index'

Outputs
-------
- Enrichment CSV per layer with columns: Layer, Z-Score, Raw P-Value, FDR, Mean t (Target Genes), Mean t (Random)

Usage
-----
    python tests/examples/example_enrichment_layer.py
"""
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

from BrainNetAnno.gene_layer import run_pipeline as gene_layer_pipeline

layer_xlsx = str(ROOT / "examples_data" / "layers.xlsx")
target_genes_csv = str(ROOT / "examples_output" / "best_genes.csv")
output_csv = str(ROOT / "examples_output" / "layer_enrichment.csv")

os.makedirs(os.path.dirname(output_csv), exist_ok=True)

res = gene_layer_pipeline(
    layer_marker_path=layer_xlsx,
    target_genes_path=target_genes_csv,
    output_path=output_csv,
)
print("Saved:", output_csv)
print(res.head())
