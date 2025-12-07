"""
PLS analysis for transcriptome gene contributions against FC deviations.

Provides reusable functions to select PLS components via CV, run PLS, compute
permutation p-values for weights, plotting utilities, and high-level pipeline.
"""
from typing import Tuple, Optional
from scipy.stats import zscore
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import os

from .utils import load_table_generic

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

try:
    from BrainNetAnno.utils import (
        optimal_pls_components_rmse,
        permutation_explained_variance,
        plot_rmse_and_variance,
    )
except ModuleNotFoundError:
    import sys, os
    src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if src_dir not in sys.path:
        sys.path.append(src_dir)
    try:
        from utils import (
            optimal_pls_components_rmse,
            permutation_explained_variance,
            plot_rmse_and_variance,
        )
    except ModuleNotFoundError:
        project_root = os.path.abspath(os.path.join(src_dir, '..'))
        if project_root not in sys.path:
            sys.path.append(project_root)
        from BrainNetAnno.utils import (
            optimal_pls_components_rmse,
            permutation_explained_variance,
            plot_rmse_and_variance,
        )

# -----------------------------
# Data reading
# -----------------------------

def read_data_from_files(t_values_file: str, gene_expression_file: str) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    """Read FC deviation (t_values) and gene contribution matrices, align by Region_Pair.

    Parameters
    ----------
    t_values_file : str
        Path to CSV containing a square matrix of FC deviations (no header).
    gene_expression_file : str
        Path to CSV containing gene contributions with a 'Region_Pair' column.

    Returns
    -------
    Tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]
        (df_pairs_values, t_values_vector, gene_matrix_aligned, gene_names_array)
    """
    t_values = load_table_generic(t_values_file, header=None)
    gene_df = load_table_generic(gene_expression_file)

    # Upper triangle non-zero entries from t_values
    coords = [(c, r) for r in range(t_values.shape[0]) for c in range(r + 1, t_values.shape[1]) if t_values.iat[r, c] != 0]
    formatted = [f"{r}-{c}" for c, r in coords]

    # Align gene_df by Region_Pair
    filtered = gene_df[gene_df["Region_Pair"].isin(formatted)].reset_index(drop=True)
    gene_names = filtered.columns.values
    gene_names = gene_names[gene_names != "Region_Pair"]
    gene_matrix = filtered.drop(columns=["Region_Pair"]).values

    df_pairs_values = pd.DataFrame([(f"{r}-{c}", t_values.at[r, c]) for c, r in coords], columns=["Region_Pair", "Value"])
    t_vector = df_pairs_values[["Value"]].values

    return df_pairs_values, t_vector, gene_matrix, gene_names

# -----------------------------
# Select component and genes
# -----------------------------

def select_best_component_and_genes(explained_variance_ratio: np.ndarray,
                                    p_values: np.ndarray,
                                    weights: np.ndarray,
                                    threshold: float = 0.05):
    """Select best component (prefer significant) and significant genes by simple z-score permutation.

    Parameters
    ----------
    explained_variance_ratio : np.ndarray
        Ratio of variance explained per component.
    p_values : np.ndarray
        P-values per component from permutation.
    weights : np.ndarray
        X weights matrix from PLS, shape (n_features, n_components).
    threshold : float, optional
        Significance threshold, by default 0.05.

    Returns
    -------
    Tuple[int, np.ndarray, np.ndarray, np.ndarray]
        best_component_index, best_genes_indices, z_scores_per_gene, p_values_per_gene
    """
    significant_components = np.where(p_values < threshold)[0]
    if len(significant_components) > 0:
        best_component = significant_components[np.argmax(explained_variance_ratio[significant_components])]
    else:
        best_component = int(np.argmax(explained_variance_ratio))

    best_weights = weights[:, best_component]
    z_scores = zscore(best_weights)
    p_values_genes = np.array([np.mean(np.random.permutation(z_scores) >= abs(z)) for z in z_scores])
    best_genes = np.where(p_values_genes < 0.05)[0]
    return best_component, best_genes, z_scores, p_values_genes

# -----------------------------
# Save and plot
# -----------------------------

def save_best_genes_to_csv(gene_names: np.ndarray,
                           best_genes: np.ndarray,
                           z_scores: np.ndarray,
                           p_values_genes: np.ndarray,
                           output_file: str) -> None:
    """Save significant genes info to CSV file."""
    results = {
        'Gene Index': gene_names[best_genes],
        'Original Z-Score': z_scores[best_genes],
        'P-Value': p_values_genes[best_genes]
    }
    results_df = pd.DataFrame(results)
    results_df.sort_values(by='P-Value', inplace=True)
    results_df.to_csv(output_file, index=False)

# -----------------------------
# High-level pipeline
# -----------------------------

def run_transcriptome_pls_pipeline(
    t_values_file: str,
    gene_expression_file: str,
    output_best_genes_csv: str,
    fig_outputfile: Optional[str] = None,
    max_components: int = 15,
    n_splits: int = 5,
    n_permutations: int = 1000,
    n_components: Optional[int] = None,
) -> None:
    """Run PLS-based transcriptome analysis linking FC deviations to gene contributions.

    This pipeline reads FC deviation (t-values) and per-connection gene
    contribution matrices, aligns them by ``Region_Pair``, standardizes inputs,
    selects the optimal number of PLS components via cross-validated RMSE,
    estimates permutation-based significance for explained variance and feature
    weights, saves the significant genes to CSV, and plots RMSE and explained
    variance curves.

    Parameters
    ----------
    t_values_file : str
        Path to the CSV containing a square matrix of FC deviations (no header).
    gene_expression_file : str
        Path to the CSV containing per-connection gene contributions with a
        ``'Region_Pair'`` column.
    output_best_genes_csv : str
        Path to save the significant genes table. Parent directories are
        created automatically.
    fig_outputfile : Optional[str], optional
        Path to save plots (RMSE curve and explained variance). If ``None``,
        plots are saved next to ``output_best_genes_csv``. Default ``None``.
    max_components : int, optional
        Maximum number of components to consider in CV. Default ``15``.
    n_splits : int, optional
        Number of CV splits. Default ``5``.
    n_permutations : int, optional
        Number of permutations for significance testing. Default ``1000``.
    n_components : Optional[int], optional
        If provided, forces the number of components; otherwise chosen by CV.

    Returns
    -------
    None
        Saves outputs to disk and logs progress.

    Notes
    -----
    - Inputs are standardized before fitting.
    - Significant component selection prefers components with permutation
      p-values below the threshold; otherwise the component with highest
      explained variance is chosen.
    - Output directories are created automatically for CSV and figures.
    """
    logger.info(f"Loading t-values from: {t_values_file}")
    logger.info(f"Loading gene contributions from: {gene_expression_file}")
    df_pairs_values, t_values, gene_matrix, gene_names = read_data_from_files(t_values_file, gene_expression_file)
    df_pairs_values.to_csv(output_best_genes_csv.replace('.csv', '_pairs_values.csv'), index=False)

    # Standardize
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    X = scaler.fit_transform(gene_matrix)
    Y = scaler.fit_transform(t_values)

    # CV errors and best components
    logger.info(f"Selecting optimal components (max={max_components}, splits={n_splits})")
    cv_errors = optimal_pls_components_rmse(X, Y, max_components=max_components, n_splits=n_splits)
    best_n_components = int(np.argmin(cv_errors) + 1) if n_components is None else int(n_components)
    logger.info(f"Best number of components selected: {best_n_components}")

    # 置换检验
    explained_variance_ratio, p_values, weights = permutation_explained_variance(X, Y, best_n_components, n_permutations=n_permutations)
    best_component, best_genes, z_scores, p_values_genes = select_best_component_and_genes(explained_variance_ratio, p_values, weights)

    # Save significant genes
    save_best_genes_to_csv(gene_names, best_genes, z_scores, p_values_genes, output_best_genes_csv)

    # Plot curves
    plot_rmse_and_variance(cv_errors, explained_variance_ratio, fig_outputfile)

    import os
    out_dir = os.path.dirname(output_best_genes_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    if fig_outputfile:
        fig_dir = os.path.dirname(fig_outputfile)
        if fig_dir:
            os.makedirs(fig_dir, exist_ok=True)

    df_best = pd.DataFrame({
        'Gene Index': gene_names[best_genes],
        'Original Z-Score': z_scores[best_genes],
        'P-Value': p_values_genes[best_genes]
    })
    df_best.to_csv(output_best_genes_csv, index=False)
    logger.info(f"Saved best genes CSV to: {output_best_genes_csv}")
    if fig_outputfile:
        logger.info(f"Saved plots to: {fig_outputfile}")
    logger.info("Transcriptome PLS pipeline completed successfully.")