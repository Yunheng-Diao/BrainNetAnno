# network-molecular

Molecular network analysis tools for transcriptome, mitochondrial, and neurotransmitter datasets.

## Features
- CPU-only NumPy/SciPy/Scikit-learn implementations
- Shared utilities in `network_molecular.utils`
- PLS pipelines: `*_pls_cge.run_pipeline`
- Enrichment modules: `gene_celltype.run_pipeline`, `gene_layer.run_pipeline`
- Command-line interface (CLI) for common workflows

## Requirements
- Python >= 3.9
- Recommended packages (specified in pyproject.toml):
  - numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib

## Installation
```
pip install network-molecular
```
Or from source:
```
pip install .
```

## Quickstart (Python API)
- Transcriptome PLS-CGE pipeline:
```python
from network_molecular import run_transcriptome_pls_cge
best_n, df = run_transcriptome_pls_cge(
    t_values_file="path/to/fc_deviation.csv",  # square matrix CSV, no header
    gene_expression_file="path/to/gene_contribution_df.csv",  # contains 'Region_Pair'
    output_best_genes_csv="path/to/best_genes.csv",
    fig_outputfile="path/to/rmse_variance.tif",
)
print(best_n)
print(df.head())
```
- Cell-type enrichment:
```python
from network_molecular import run_celltype_enrichment
res = run_celltype_enrichment(
    celltype_csv="path/to/celltypes_PSP.csv",   # columns: gene, class
    target_genes_csv="path/to/best_genes_SCZ.csv",  # column: Gene Index
    output_csv="path/to/celltype_enrichment_results.csv",
)
print(res.head())
```
- Layer enrichment:
```python
from network_molecular import run_layer_enrichment
res = run_layer_enrichment(
    layer_marker_path="path/to/41593_2020_787_MOESM3_ESM.xlsx",  # sheet: Table S4B
    target_genes_path="path/to/best_genes_all.csv",
    output_csv="path/to/ALL_gene_layer_analysis_results.csv",
)
print(res.head())
```

## CLI Usage
After installation, the following commands are available:
- Cell-type enrichment
```
network-molecular-celltype path/to/celltypes.csv path/to/target_genes.csv --output path/to/out.csv --n-perm 5000 --seed 42
```
- Layer enrichment
```
network-molecular-layer path/to/layer_markers.xlsx path/to/target_genes.csv --output path/to/out.csv --n-perm 5000 --seed 42
```
- Transcriptome PLS-CGE
```
network-molecular-pls path/to/fc_deviation.csv path/to/gene_contrib.csv --output-best-genes path/to/out.csv --fig path/to/out.tif
```

## Data Format
- FC deviation CSV: square numeric matrix without header; only upper-triangle non-zero entries are used.
- Gene contribution CSV: must include a `Region_Pair` column formatted as `i-j`, where i < j.
- Cell-type markers CSV: columns `gene` and `class` (gene symbols are uppercased internally).
- Layer markers Excel: sheet contains columns named `t_stat_LayerX` (e.g., Layer1..Layer6), and a `gene` column.

## Development & Testing
- Build package:
```
python -m pip install -U build
python -m build
```
- Install wheel locally:
```
pip install dist/network_molecular-<version>-py3-none-any.whl
```
- Run tests (after install or from source):
```
pip install -U pytest
pytest -q
```

## Versioning
- Semantic versioning (MAJOR.MINOR.PATCH)

## License
MIT
