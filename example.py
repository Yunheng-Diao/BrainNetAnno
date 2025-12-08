# When running as a standalone script, add the src directory to sys.path to avoid relative import issues
# Notes:
# - This example script demonstrates typical usage of BrainNetAnno modules and required data formats.
# - To run without installing into site-packages, if a local src layout exists, we prepend it to sys.path.
# - Paths below are examples; replace them with your local data locations.
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if SRC.exists():
    # Put source directory at the front of import search path to prefer local code over installed package
    sys.path.insert(0, str(SRC))

# Absolute imports (package-level public APIs)
# - run_transcriptome_pipeline: compute transcriptome contributions (distance-correlation fit, expected matrix, per-connection gene contributions)
# - run_mitochondrial_pipeline: compute mitochondrial feature contributions (precision matrix, distance decay, LOFO)
# - run_neurotransmitter_pipeline: compute neurotransmitter feature contributions (similar workflow)
# - run_*_pls_pipeline: PLS analysis (component selection, fit, permutation tests, plotting and saving)
# - gene_celltype_pipeline / gene_layer_pipeline: enrichment analyses for target gene sets by cell types/layers
from BrainNetAnno.transcriptome import run_transcriptome_pipeline
from BrainNetAnno.mitochondrial import run_mitochondrial_pipeline
from BrainNetAnno.neurotransmitter import run_neurotransmitter_pipeline
from BrainNetAnno.transcriptome_pls_cge import run_transcriptome_pls_pipeline
from BrainNetAnno.mitochondrial_pls_cge import run_mitochondrial_pls_pipeline
from BrainNetAnno.neurotransmitter_pls_cge import run_neurotransmitter_pls_pipeline
from BrainNetAnno.gene_celltype import run_pipeline as gene_celltype_pipeline
from BrainNetAnno.gene_layer import run_pipeline as gene_layer_pipeline

print("Running example pipeline...")
# Step 1: Transcriptome contribution calculation
# Inputs:
# - coords_csv: coordinate table (CSV/Excel) with columns MNI_X/MNI_Y/MNI_Z
# - genes_csv: expression matrix (CSV/Excel), rows=regions, cols=genes. A 'label' column is allowed and will be dropped.
# Output:
# - transcriptome_contribution_path: per-connection gene contribution table with 'Region_Pair' and gene contribution columns
coords_csv = r"E:\aliyun_backup\muilt_disorders\01_altas\BNA246_Cerebellum_Stem_Asym.csv"
genes_csv = r"E:\aliyun_backup\muilt_disorders\13_Allen\AHBA_expression_BNA246.csv"
transcriptome_contribution_path = r"C:\Users\diaoy\Desktop\test\gene_contribution.csv"


run_transcriptome_pipeline(
    coordinates_path=coords_csv,
    gene_expression_path=genes_csv,
    output_contribution_path=transcriptome_contribution_path,
    initial_params=(0.64, 90.4, -0.19),  # initial parameters (A, n, B) for distance-correlation exponential decay
    save_plot=True,                       # whether to save the fit visualization
    plot_path=None,                       # if None, no figure is saved (only show or internal default)
)

print("Step 2: Transcriptome PLS analysis")
# Step 2: Transcriptome PLS analysis
# Inputs:
# - t_values_file_path: square CSV (no header), FC deviations/weights; only upper-triangle non-zero entries are used
# - gene_expression_file_path: per-connection gene contributions from Step 1 (must contain 'Region_Pair')
# Output:
# - output_best_genes_path: significant genes in the best component (with z-scores and permutation p-values)
# - mse_plot_path / explained_variance_plot_path: CV MSE curve and explained variance plots
output_best_genes_csv = r"C:\Users\diaoy\Desktop\test\best_genes_SCZ.csv"
transcriptome_pls_mse_plot_path = r"C:\Users\diaoy\Desktop\test\transcriptome_pls_selection.tif"
transcriptome_pls_explained_variance_plot_path = r"C:\Users\diaoy\Desktop\test\transcriptome_explained_variance.tif"

run_transcriptome_pls_pipeline(
    t_values_file_path=r"E:\aliyun_backup\muilt_disorders\11_pls\result\pls_weight\偏离连接_SCZ_top10.csv",
    gene_expression_file_path=transcriptome_contribution_path,
    output_best_genes_path=output_best_genes_csv,
    mse_plot_path=transcriptome_pls_mse_plot_path,
    explained_variance_plot_path=transcriptome_pls_explained_variance_plot_path,
    max_components=10,      # upper bound of components to search; internally further capped by data dimensions
    cv_splits=5,            # KFold splits
    n_permutations=1000,    # number of permutations for significance evaluation
)


print("Step 3: Gene cell-type enrichment analysis")
# Step 3: Cell-type enrichment analysis
# Inputs:
# - celltype_path: cell-type marker table (CSV) with at least 'gene' and 'class' columns
# - target_genes_path: significant genes from Step 2 (CSV) with 'Gene Index' column
# Output:
# - celltype_output_csv: overlap enrichment statistics per cell type (Z-score, raw P, FDR, etc.)
celltype_output_csv = r"C:\Users\diaoy\Desktop\test\celltype_enrichment_results.csv"

gene_celltype_pipeline(
    celltype_path=r"E:\aliyun_backup\muilt_disorders\16_celltype\celltypes_PSP.csv",
    target_genes_path=output_best_genes_csv,
    output_path=celltype_output_csv,
)

print("Step 4: Gene layer analysis")
# Step 4: Cortical layer enrichment analysis
# Inputs:
# - layer_marker_path: layer marker statistics (Excel), default sheet 'Table S4B'; includes 'gene' and per-layer t-stat columns
# - target_genes_path: same significant gene list
# Output:
# - layer_output_csv: enrichment statistics for each layer (Z-score, raw P, FDR, etc.)
layer_output_csv = r"C:\Users\diaoy\Desktop\test\gene_layer_analysis_results.csv"

gene_layer_pipeline(
    layer_marker_path=r"E:\aliyun_backup\muilt_disorders\16_celltype\41593_2020_787_MOESM3_ESM.xlsx",
    target_genes_path=output_best_genes_csv,
    output_path=layer_output_csv,
)

print("Step 5: Neurotransmitter contribution calculation")
# Step 5: Neurotransmitter contribution calculation
# Inputs:
# - coordinates_path: coordinates table (same as Step 1)
# - neurotransmitter_expression_path: expression matrix (CSV/Excel), rows=regions, cols=transmitters or receptors
# Output:
# - neurotransmitter_contribution_path: per-connection transmitter contribution table with 'Region_Pair' and feature columns
neurotransmitter_xlsx = r"E:\aliyun_backup\muilt_disorders\14_neurotransmitter\neurotransmitter_BNA246.xlsx"
neurotransmitter_contribution_path = r"C:\Users\diaoy\Desktop\test\Neurotransmitter_contribution.csv"

run_neurotransmitter_pipeline(
    coordinates_path=coords_csv,
    neurotransmitter_expression_path=neurotransmitter_xlsx,
    output_contribution_path=neurotransmitter_contribution_path,
    initial_params=(1.0, 50.0, 0.0),  # initial parameters (A, n, B) for distance-connection fitting
    save_plot=True,
    plot_path=None,
)

# Step 6: Neurotransmitter PLS analysis
# Inputs:
# - fc_matrix_path: square FC weights/deviations (CSV, no header), upper-triangle non-zero entries
# - nt_contrib_path: per-connection transmitter contributions (must include 'Region_Pair')
# Output:
# - Optional: MSE curve, explained variance plot, CSV of weights and permutation p-values
neurotransmitter_pls_mse_plot_path = r"C:\Users\diaoy\Desktop\test\neurotransmitter_pls_selection.tif"
neurotransmitter_pls_explained_variance_plot_path = r"C:\Users\diaoy\Desktop\test\neurotransmitter_explained_variance.tif"
neurotransmitter_pls_output_weights_path = r"C:\Users\diaoy\Desktop\test\neurotransmitter_pls_weights_with_pvals.csv"

print("Step 6: Neurotransmitter PLS analysis")
run_neurotransmitter_pls_pipeline(
    fc_matrix_path=r"E:\aliyun_backup\muilt_disorders\11_pls\result\pls_weight\偏离连接_SCZ_01.csv",
    nt_contrib_path=neurotransmitter_contribution_path,
    top_k=350,                          # select top-k by |FC| to reduce noise and computation
    mse_plot_path=neurotransmitter_pls_mse_plot_path,
    explained_variance_plot_path=neurotransmitter_pls_explained_variance_plot_path,
    output_weights_path=neurotransmitter_pls_output_weights_path,
    max_components=10,
    cv_splits=5,
    random_state=42,
    n_permutations=1000,
)

print("Step 7: Mitochondrial contribution calculation")
# Step 7: Mitochondrial contribution calculation (similar to neurotransmitter workflow)
mitochondrial_xlsx = r"E:\aliyun_backup\muilt_disorders\15_Mitochondrial\Mitochondrial.xlsx"
mitochondrial_contribution_path = r"C:\Users\diaoy\Desktop\test\mitochondrial_contribution.csv"

run_mitochondrial_pipeline(
    coordinates_path=coords_csv,
    mitochondrial_expression_path=mitochondrial_xlsx,
    output_contribution_path=mitochondrial_contribution_path,
    initial_params=(1.0, 50.0, 0.0),
    save_plot=True,
    plot_path=None,
)

# Step 8: Mitochondrial PLS analysis
# Similar to Step 6; outputs MSE/explained variance plots and permutation significance for weights
mitochondrial_pls_mse_plot_path = r"C:\Users\diaoy\Desktop\test\mitochondrial_pls_selection.tif"
mitochondrial_pls_explained_variance_plot_path = r"C:\Users\diaoy\Desktop\test\mitochondrial_explained_variance.tif"
mitochondrial_pls_output_weights_path = r"C:\Users\diaoy\Desktop\test\mitochondrial_pls_weights_with_pvals.csv"

print("Step 8: Mitochondrial PLS analysis")
run_mitochondrial_pls_pipeline(
    fc_matrix_path=r"E:\aliyun_backup\muilt_disorders\11_pls\result\pls_weight\偏离连接_SCZ_top10.csv",
    nt_contrib_path=mitochondrial_contribution_path,
    top_k=500,
    mse_plot_path=mitochondrial_pls_mse_plot_path,
    explained_variance_plot_path=mitochondrial_pls_explained_variance_plot_path,
    output_weights_path=mitochondrial_pls_output_weights_path,
    max_components=6,
    cv_splits=5,
    random_state=42,
    n_permutations=1000,
)
