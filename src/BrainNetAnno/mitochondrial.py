"""
Mitochondrial network analysis utilities using NumPy/SciPy/Scikit-learn on CPU.

This module provides functions to build precision matrices, fit distance-decay
models, compute residual networks, and perform leave-one-feature-out (LOFO)
contribution analysis for mitochondrial gene sets.
"""
from typing import Tuple, Sequence, Iterable, Optional
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf
import matplotlib.pyplot as plt
import os

# Import shared utilities with robust fallback for script execution
try:
    from .utils import (
        load_coordinates,
        load_gene_expression_excel,
        zscore_rows,
        compute_distance_matrix,
        exponential_decay,
        fit_exponential_decay,
        compute_expected_matrix,
    )
except ImportError:
    import os, sys
    # Add the parent directory of this file (i.e., the package root's parent) to sys.path
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from BrainNetAnno.utils import (
        load_coordinates,
        load_gene_expression_excel,
        zscore_rows,
        compute_distance_matrix,
        exponential_decay,
        fit_exponential_decay,
        compute_expected_matrix,
    )

# -----------------------------
# Core computations (module-specific)
# -----------------------------

def build_precision_matrix(z_mitochondrial_expression: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Estimate covariance and precision for mitochondrial expression.

    Parameters
    ----------
    z_mitochondrial_expression : np.ndarray
        Z-scored gene expression matrix of shape (n_regions, n_genes). Rows are regions.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (covariance, precision) matrices of shape (n_genes, n_genes).
    """
    model = LedoitWolf()
    model.fit(z_mitochondrial_expression.T)
    cov = model.covariance_
    precision = np.linalg.inv(cov)
    return cov, precision


def compute_residual_matrix(observed_matrix: np.ndarray, expected_matrix: np.ndarray) -> np.ndarray:
    """Compute residual network as observed minus expected.

    Parameters
    ----------
    observed_matrix : np.ndarray
        Observed connection/precision matrix.
    expected_matrix : np.ndarray
        Expected connection matrix computed from distance decay.

    Returns
    -------
    np.ndarray
        Residual matrix of the same shape as inputs.
    """
    return observed_matrix - expected_matrix


def lofo_contribution(z_mitochondrial_expression: np.ndarray,
                      expected_matrix: np.ndarray,
                      region_pairs: Iterable[Tuple[int, int]],
                      base_precision: np.ndarray,
                      gene_names: Sequence[str]) -> pd.DataFrame:
    """Perform Leave-One-Feature-Out (LOFO) analysis to estimate gene contributions.

    For each gene g, remove it from the expression matrix, rebuild the precision
    matrix, and compute the change in residuals per connection.

    Parameters
    ----------
    z_mitochondrial_expression : np.ndarray
        Z-scored gene expression matrix of shape (n_regions, n_genes).
    expected_matrix : np.ndarray
        Expected connection matrix of shape (n_regions, n_regions).
    region_pairs : Iterable[Tuple[int, int]]
        Iterable of (i, j) indices representing upper-triangle connections.
    base_precision : np.ndarray
        Precision matrix computed from the full set of genes (shape (n_genes, n_genes)).
    gene_names : Sequence[str]
        Names of genes corresponding to columns in z_mitochondrial_expression.

    Returns
    -------
    pd.DataFrame
        DataFrame containing per-connection contribution scores per gene with column 'Region_Pair'.
    """
    region_pairs = list(region_pairs)
    num_connections = len(region_pairs)
    num_genes = z_mitochondrial_expression.shape[1]

    # Residual with full gene set
    residual_full = compute_residual_matrix(base_precision, expected_matrix)

    out = pd.DataFrame(np.zeros((num_connections, num_genes)), columns=gene_names)

    for g in range(num_genes):
        # Remove gene g
        reduced_expression = np.delete(z_mitochondrial_expression, g, axis=1)
        try:
            _, precision_lofo = build_precision_matrix(reduced_expression)
        except np.linalg.LinAlgError:
            # Skip if covariance is singular
            continue

        residual_lofo = compute_residual_matrix(precision_lofo, expected_matrix)

        # Delta residual per connection
        for idx, (i, j) in enumerate(region_pairs):
            delta = residual_full[i, j] - residual_lofo[i, j]
            out.iloc[idx, g] = delta

    out.insert(0, 'Region_Pair', [f'{i}-{j}' for i, j in region_pairs])
    return out

# -----------------------------
# Visualization
# -----------------------------

def plot_distance_decay(distances: np.ndarray, conn_values: np.ndarray, params: Tuple[float, float, float]) -> None:
    """Plot empirical connection strength vs. distance and fitted exponential decay curve.

    Parameters
    ----------
    distances : np.ndarray
        1D array of pairwise distances (upper triangle entries).
    conn_values : np.ndarray
        1D array of observed connection strengths matching distances.
    params : Tuple[float, float, float]
        Fitted parameters (A, n, B).
    """
    A, n, B = params
    plt.figure(figsize=(8, 6))
    plt.scatter(distances, conn_values, s=4, color='gray', alpha=0.6, label='Data')
    order = np.argsort(distances)
    plt.plot(distances[order], exponential_decay(distances[order], A, n, B),
             color='red', linewidth=2, label='Fitted Curve')
    plt.xlabel("Distance between regions (mm)")
    plt.ylabel("Connection strength (precision)")
    plt.title("Connection vs. Distance with Fitted Exponential Decay")
    plt.legend()

# -----------------------------
# High-level pipeline
# -----------------------------

def run_mitochondrial_pipeline(coordinates_csv: str,
                 gene_expression_xlsx: str,
                 output_contribution_csv: str,
                 initial_params: Sequence[float] = (1.0, 50.0, 0.0),
                 save_plot: bool = True,
                 plot_path: Optional[str] = None) -> Tuple[Tuple[float, float, float], pd.DataFrame]:
    """Run mitochondrial CGE/precision analysis and optionally save LOFO gene contributions.

    This high-level pipeline loads region coordinates and mitochondrial gene
    expression from Excel, computes distances and a precision-based connection
    measure, fits an exponential distance-decay model, builds an expected
    connection matrix, and performs leave-one-feature-out (LOFO) analysis to
    quantify per-gene contributions per region pair. Results can be saved to
    CSV and an optional distance-decay plot can be saved.

    Parameters
    ----------
    coordinates_csv : str
        Path to the CSV containing region coordinates (columns: ``MNI_X``,
        ``MNI_Y``, ``MNI_Z``).
    gene_expression_xlsx : str
        Path to the Excel file containing mitochondrial gene expression
        (rows=regions, columns=genes). Label columns may be dropped upstream.
    output_contribution_csv : str
        Path to save the per-connection LOFO gene contribution scores as CSV.
        Parent directories are created automatically.
    initial_params : Sequence[float], optional
        Initial guess for decay parameters ``(A, n, B)``, by default
        ``(1.0, 50.0, 0.0)``.
    save_plot : bool, optional
        Whether to save the fitted distance-decay plot. Default ``True``.
    plot_path : Optional[str], optional
        Path to save the plot when ``save_plot=True``. If ``None``, a file next
        to ``output_contribution_csv`` is used.

    Returns
    -------
    Tuple[Tuple[float, float, float], pd.DataFrame]
        The fitted decay parameters ``(A, n, B)`` and the LOFO contribution
        DataFrame with a ``'Region_Pair'`` column followed by per-gene scores.

    Notes
    -----
    - Rows of the expression matrix are z-scored prior to precision modeling.
    - Precision matrix is estimated via Ledoit–Wolf shrinkage.
    - Expected matrix is derived from the fitted decay parameters.
    - Parent directories for ``output_contribution_csv`` are created
      automatically when saving results.
    """
    # Load data
    coords = load_coordinates(coordinates_csv)
    matrix, gene_names = load_gene_expression_excel(gene_expression_xlsx)

    # Preprocess
    dist_mat = compute_distance_matrix(coords)
    z_mat = zscore_rows(matrix)

    # Precision from full genes
    _, precision_full = build_precision_matrix(z_mat)

    # Upper triangle indices for region pairs
    upper = np.triu_indices(dist_mat.shape[0], k=1)
    distances = dist_mat[upper]

    # Connection values per pair from precision matrix: use entries aligned with region pairs
    conn_values = np.array([precision_full[i, j] for i, j in zip(upper[0], upper[1])])

    # Fit exponential decay
    params = fit_exponential_decay(distances, conn_values, initial_params=initial_params)
    A, n, B = params

    # Expected matrix
    expected = compute_expected_matrix(dist_mat, A, n, B)

    # LOFO contributions
    contrib_df = lofo_contribution(z_mat, expected, zip(upper[0], upper[1]), precision_full, gene_names)

    # Save contributions
    d = os.path.dirname(output_contribution_csv)
    if d:
        os.makedirs(d, exist_ok=True)
    contrib_df.to_csv(output_contribution_csv, index=False)

    # Optional plot
    if save_plot:
        if plot_path is None:
            plot_path = os.path.splitext(output_contribution_csv)[0] + "_decay_plot.png"
        pd = os.path.dirname(plot_path)
        if pd:
            os.makedirs(pd, exist_ok=True)
        plt.figure(figsize=(8, 6))
        order = np.argsort(distances)
        plt.scatter(distances, conn_values, s=3, color='gray', alpha=0.5)
        plt.plot(distances[order], exponential_decay(distances[order], A, n, B), color='red', linewidth=2)
        plt.xlabel("Distance between regions (mm)")
        plt.ylabel("Connection strength (precision)")
        plt.title("Connection vs. Distance with Fitted Exponential Decay")
        plt.savefig(plot_path, dpi=300)
        plt.close()

    return params, contrib_df