"""
Example: Cell-type enrichment analysis

Functionality
-------------
Permutation-based enrichment test to assess whether a set of target genes is
overrepresented within cell-type marker sets.

Inputs
------
- celltype_path: CSV with columns 'gene' and 'class'
- target_genes_path: CSV with column 'Gene Index'

Outputs
-------
- Enrichment CSV with columns: Cell Type, Overlap, Z-Score, Raw P-Value, Mean Overlap (Random), FDR

Usage
-----
    python tests/examples/example_enrichment_celltype.py
"""
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

from BrainNetAnno.gene_celltype import run_pipeline as gene_celltype_pipeline

celltype_csv = str(ROOT / "examples_data" / "celltypes.csv")
target_genes_csv = str(ROOT / "examples_output" / "best_genes.csv")
output_csv = str(ROOT / "examples_output" / "celltype_enrichment.csv")

os.makedirs(os.path.dirname(output_csv), exist_ok=True)

res = gene_celltype_pipeline(
    celltype_path=celltype_csv,
    target_genes_path=target_genes_csv,
    output_path=output_csv,
)
print("Saved:", output_csv)
print(res.head())
