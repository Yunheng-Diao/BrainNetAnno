# BrainNetAnno

A pure-Python toolkit for molecular annotation of brain networks, integrating transcriptome, neurotransmitters, and mitochondria. It provides PLS-CGE (Partial Least Squares – Co-Gene Expression) pipelines and enrichment analyses to explain functional connectivity (FC) changes at the molecular level.

## Table of Contents
- Motivation & Background
- Features
- Environment & Installation
- Quickstart (Python API)
- Data Format & Preparation
- Modules & API
- Advanced Examples
- FAQ & Troubleshooting
- Development & Testing
- License

## Motivation & Background
In neuroimaging, FC deviations or weights are often used to describe connectivity changes under disease or task conditions. BrainNetAnno integrates transcriptomic, neurotransmitter, and mitochondrial molecular features, using PLS-CGE to map these features to FC changes, enabling you to:
- Identify key genes, neurotransmitter systems, or mitochondrial modules associated with FC changes.
- Provide molecular biological explanations for network-level alterations.
- Support downstream enrichment analyses of cell types and cortical layers.

## Features
- CPU-only implementation relying on NumPy, SciPy, scikit-learn, statsmodels, and matplotlib.
- Shared utility functions in `BrainNetAnno.utils`.
- PLS-CGE pipelines:
  - `transcriptome_pls_cge.run_transcriptome_pls_pipeline`
  - `neurotransmitter_pls_cge.run_neurotransmitter_pls_pipeline`
  - `mitochondrial_pls_cge.run_mitochondrial_pls_pipeline`
- Enrichment analyses:
  - Cell-type enrichment: `gene_celltype.run_pipeline`
  - Cortical layer enrichment: `gene_layer.run_pipeline`
- Export results: weights, best genes, metric plots, etc.

## Environment & Installation
- Python >= 3.8
- Main dependencies in `pyproject.toml` (numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib)

Install from source (recommended in a virtual environment):
```
pip install BrainNetAnno
```
Or build and install a wheel:
```
python -m pip install -U build
python -m build
pip install dist/BrainNetAnno-<version>-py3-none-any.whl
```

## Quickstart (Python API)
Minimal examples for typical tasks. See the next section for data format requirements.

- Transcriptome PLS-CGE pipeline:
```python
from BrainNetAnno.transcriptome_pls_cge import run_transcriptome_pls_pipeline

run_transcriptome_pls_pipeline(
    t_values_file="path/to/fc_deviation.csv",          # Square CSV, no header; only upper-triangle non-zeros are used
    gene_expression_file="path/to/gene_contrib.csv",   # Must include column 'Region_Pair' like 'i-j' with i<j
    output_best_genes_csv="path/to/best_genes.csv",    # Output best genes
    fig_outputfile="path/to/rmse_variance.tif",        # Output RMSE/variance curves figure
)
```

- Neurotransmitter PLS pipeline:
```python
from BrainNetAnno.neurotransmitter_pls_cge import run_neurotransmitter_pls_pipeline

best_n, df = run_neurotransmitter_pls_pipeline(
    fc_matrix_path="path/to/fc_deviation.csv",      # Square CSV, no header
    nt_contrib_csv="path/to/nt_contrib.csv",        # Must include 'Region_Pair'
    output_weights_csv="path/to/nt_weights.csv",    # Output neurotransmitter weights
)
print(best_n)
print(df.head())
```

- Mitochondrial PLS pipeline:
```python
from BrainNetAnno.mitochondrial_pls_cge import run_mitochondrial_pls_pipeline

best_n, df = run_mitochondrial_pls_pipeline(
    fc_matrix_path="path/to/fc_deviation.csv",      # Square CSV, no header
    nt_contrib_csv="path/to/mito_contrib.csv",      # Must include 'Region_Pair'
    output_weights_csv="path/to/mito_weights.csv",  # Output mitochondrial weights
)
print(best_n)
print(df.head())
```

- Cell-type enrichment:
```python
from BrainNetAnno.gene_celltype import run_pipeline

res = run_pipeline(
    celltype_csv="path/to/celltypes_PSP.csv",      # Columns: gene, class; gene symbols are uppercased internally
    target_genes_csv="path/to/best_genes_SCZ.csv", # Columns: Gene Index or Gene; supports index or symbol
    output_csv="path/to/celltype_enrichment_results.csv",
)
print(res.head())
```

- Cortical layer enrichment:
```python
from BrainNetAnno.gene_layer import run_pipeline

res = run_pipeline(
    layer_marker_path="path/to/41593_2020_787_MOESM3_ESM.xlsx",  # sheet: Table S4B; requires t_stat_LayerX and gene columns
    target_genes_path="path/to/best_genes_all.csv",              # Columns: Gene Index or Gene
    output_csv="path/to/ALL_gene_layer_analysis_results.csv",
)
print(res.head())
```

## Data Format & Preparation
- FC deviation/weights CSV:
  - N×N square numeric matrix, no header; diagonal may be 0; only upper-triangle non-zero entries are used.
  - Example: `fc_deviation.csv`, entry (i, j) indicates FC change between regions i and j (i<j).
- Gene/NT/Mito contribution CSV:
  - Must include `Region_Pair` formatted as `i-j` with i<j (aligned to the FC matrix indexing).
  - Other columns contain feature values (e.g., gene contribution scores or NT class intensities).
- Cell-type markers CSV:
  - Must include `gene` and `class` columns; gene symbols will be uppercased internally.
- Cortical layer markers Excel:
  - Specify the sheet (e.g., `Table S4B`), containing `t_stat_Layer1..Layer6` and `gene` columns.

Data consistency tips:
- Ensure `Region_Pair` matches the same region index system used by the FC matrix.
- Handle missing values beforehand (clean or impute appropriately).
- Keep units/scales consistent to avoid instability in modeling.

## Modules & API
- `BrainNetAnno.transcriptome_pls_cge`
  - `run_transcriptome_pls_pipeline(t_values_file, gene_expression_file, output_best_genes_csv, fig_outputfile)`
  - Purpose: PLS modeling between FC deviations and gene co-expression; outputs best genes and RMSE/variance curves.
- `BrainNetAnno.neurotransmitter_pls_cge`
  - `run_neurotransmitter_pls_pipeline(fc_matrix_path, nt_contrib_csv, output_weights_csv)`
  - Purpose: PLS modeling between FC deviations and neurotransmitter contributions; returns best component count and weights table.
- `BrainNetAnno.mitochondrial_pls_cge`
  - `run_mitochondrial_pls_pipeline(fc_matrix_path, nt_contrib_csv, output_weights_csv)`
  - Purpose: PLS modeling between FC deviations and mitochondrial features; returns best component count and weights table.
- `BrainNetAnno.gene_celltype`
  - `run_pipeline(celltype_csv, target_genes_csv, output_csv)`
  - Purpose: Cell-type enrichment for target genes; returns enrichment statistics and significance.
- `BrainNetAnno.gene_layer`
  - `run_pipeline(layer_marker_path, target_genes_path, output_csv)`
  - Purpose: Cortical layer (Layer1..Layer6) enrichment for target genes.
- `BrainNetAnno.utils`
  - Common utilities for I/O, index mapping, upper-triangle extraction, normalization, and statistics.

## Advanced Examples
- Cross-validation and parameter selection:
```python
# Example: read and visualize neurotransmitter PLS weights
from BrainNetAnno.neurotransmitter_pls_cge import run_neurotransmitter_pls_pipeline
import pandas as pd

best_n, weights = run_neurotransmitter_pls_pipeline(
    fc_matrix_path="data/fc_deviation.csv",
    nt_contrib_csv="data/nt_contrib.csv",
    output_weights_csv="outputs/nt_weights.csv",
)
print(f"Best latent components: {best_n}")

# Load exported weights and sort
weights_df = pd.read_csv("outputs/nt_weights.csv")
print(weights_df.sort_values(by="weight", ascending=False).head(10))
```

- Use best genes for downstream enrichment:
```python
from BrainNetAnno.transcriptome_pls_cge import run_transcriptome_pls_pipeline
from BrainNetAnno.gene_celltype import run_pipeline as run_celltype
from BrainNetAnno.gene_layer import run_pipeline as run_layer

best_genes_csv = "outputs/best_genes.csv"
run_transcriptome_pls_pipeline(
    t_values_file="data/fc_deviation.csv",
    gene_expression_file="data/gene_contrib.csv",
    output_best_genes_csv=best_genes_csv,
    fig_outputfile="outputs/rmse_variance.tif",
)

# Perform cell-type and layer enrichment
celltype_res = run_celltype(
    celltype_csv="markers/celltypes_PSP.csv",
    target_genes_csv=best_genes_csv,
    output_csv="outputs/celltype_enrichment.csv",
)
layer_res = run_layer(
    layer_marker_path="markers/41593_2020_787_MOESM3_ESM.xlsx",
    target_genes_path=best_genes_csv,
    output_csv="outputs/layer_enrichment.csv",
)
print(celltype_res.head())
print(layer_res.head())
```

## FAQ & Troubleshooting
- Error "Region_Pair missing or malformed":
  - Ensure `Region_Pair` exists and follows `i-j` with i<j, aligned to the FC matrix indices.
- Error "Matrix dimension mismatch":
  - Confirm the FC matrix is N×N and the feature file covers matching upper-triangle pairs.
  - If some pairs are missing, the code intersects available pairs, but complete coverage is recommended.
- All results are NaN or insignificant:
  - Check if data are all zeros or missing; consider normalization.
  - Adjust the number of PLS components or CV folds (if exposed by the function).

## Development & Testing
- Build and install locally:
```
python -m pip install -U build
python -m build
pip install dist/BrainNetAnno-<version>-py3-none-any.whl
```
- Examples & tests:
  - See `tests/examples/` including:
    - `example_transcriptome.py` / `example_transcriptome_pls.py`
    - `example_neurotransmitter.py` / `example_neurotransmitter_pls.py`
    - `example_mitochondrial.py` / `example_mitochondrial_pls.py`
    - `example_enrichment_celltype.py` / `example_enrichment_layer.py`

## License
MIT
