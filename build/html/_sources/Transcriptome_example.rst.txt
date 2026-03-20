Transcriptome Example
================

A complete pipeline for transcriptome analysis is demonstrated in this example.

1. Transcriptome-based network annotation — downloading Allen Human Brain Atlas data
-----------------------------------------------------------------------------------

Install the `abagen` package first:

.. code-block:: bash

    pip install abagen

After installation, you can download Allen Human Brain Atlas microarray data and map it to a standard brain atlas:

.. code-block:: python

    import abagen
    # Download Allen microarray data
    abagen.datasets.fetch_microarray(donors='all')  # fetch data for all donors
    import pandas as pd
    import inspect

    # 兼容 DataFrame.set_axis(..., inplace=...)
    if 'inplace' not in inspect.signature(pd.DataFrame.set_axis).parameters:
        _orig_set_axis = pd.DataFrame.set_axis
        def _set_axis_compat(self, labels, axis=0, inplace=None):
            return _orig_set_axis(self, labels, axis=axis)
        pd.DataFrame.set_axis = _set_axis_compat

    # 兼容 DataFrame.append(...)
    if not hasattr(pd.DataFrame, 'append'):
        def _df_append(self, other, ignore_index=False, verify_integrity=False, sort=False):
            # pandas<2 的 append 等价于 concat 按行拼接
            res = pd.concat([self, other], axis=0, ignore_index=ignore_index, sort=sort)
            return res
        pd.DataFrame.append = _df_append

    # 兼容 Series.append(...)
    if not hasattr(pd.Series, 'append'):
        def _ser_append(self, to_append, ignore_index=False, verify_integrity=False):
            res = pd.concat([self, to_append], axis=0, ignore_index=ignore_index)
            return res
        pd.Series.append = _ser_append
        
    expression = abagen.get_expression_data(
        atlas="/path/to/standard_brain_atlas.nii.gz",  # path to the template atlas (e.g., BNA246, AAL)
        missing='centroids',                            # strategy to fill missing data; see abagen docs
        data_dir="/path/to/downloaded_data",          # where to store downloaded files
        verbose=True                                    # show download and processing logs
    )

Save the mapped expression table to CSV:

.. code-block:: python

    expression.to_csv('/path/to/gene_expression.csv')

.. note::

    Downloading the full Allen microarray dataset (~4 GB) may take a long time depending on your network connection. Refer to the `abagen` documentation for additional options and caching behavior.

2. Computing gene co-expression and distance dependence
-------------------------------------------------------

One key challenge when annotating connections (edges) is aligning regional gene expression with inter-regional connectivity. BrainNetAnno follows an edge-centric approach that quantifies each gene's contribution to connection strength by modeling the relationship between inter-regional spatial distance and gene co-expression (CGE). The main steps are:

1) Data preparation

- Load MNI coordinates for atlas regions and compute the pairwise Euclidean distance matrix.
- Load Allen Brain Atlas (AHBA) expression matrix (regions × genes).
- Z-score normalize expression values row-wise (across genes) to remove scale differences.

2) Compute gene co-expression

Compute the Pearson correlation (or other similarity metric) between region-wise gene expression profiles to obtain a CGE matrix (regions × regions), reflecting global similarity in gene expression between region pairs.

3) Fit distance–CGE relationship

- Extract upper-triangle elements from the distance and CGE matrices to create paired observations.
- Fit an exponential decay model to describe CGE as a function of distance

    .. math::

        CGE(d) = A \times \exp(-d / B) + C

  where d is inter-regional distance and A, B, C are fitted parameters. Use nonlinear least squares to estimate parameters and evaluate model fit (R^2).

.. note::

    Adjacent brain regions tend to exhibit similar gene-expression profiles. This spatial autocorrelation can inflate apparent associations between regional gene expression and neuroimaging features. It is therefore important to account for spatial autocorrelation to avoid biased inference.

4) Compute per-gene contributions to edges

Using the edge-centric framework, compute for each connection (region pair i-j) the contribution of each gene to the connection's expression-based similarity, yielding a contribution matrix where rows correspond to region pairs (upper triangle) and columns correspond to genes.

5) Save results

Save the per-connection gene contribution matrix to CSV for downstream PLS analysis.

Example call (high-level API):

.. code-block:: python

    run_transcriptome_pipeline(
        coordinates_path="/path/to/coords.csv",
        gene_expression_path="/path/to/gene_expression.csv",
        output_contribution_path="/path/to/transcriptome_contribution.csv",
        initial_params=(0.64, 90.4, -0.19),
        save_plot=False,
        plot_path=None,
    )

3. Connectome–transcriptome association via PLS-CGE
---------------------------------------------------

Identify genes most strongly associated with network abnormalities using Partial Least Squares (PLS) regression combined with permutation testing and multiple-testing correction. Key steps:

1) Data preparation

- Load the connectome difference (or t-values / weights) matrix and the per-connection gene contribution matrix.
- Extract target connections (e.g., non-zero upper-triangle entries or top-N by absolute value) and select corresponding gene-contribution rows.
- Standardize features and target values using `StandardScaler`.

2) Select optimal PLS components

- Use K-fold cross-validation (e.g., 5-fold) to try different numbers of PLS components.
- Fit PLS on the training folds and evaluate prediction MSE on validation folds; choose the component number minimizing average MSE.

3) Permutation testing for significance

- Fit the optimal PLS model on the full data and compute explained variance per component.
- Perform permutations (e.g., 1000): shuffle the target vector, refit the model, and record explained variance to build null distributions.
- Compute empirical P-values and identify significant components (e.g., P < 0.05).

4) Identify significant genes

- Extract weights for the optimal component and compute z-scores for gene weights using the permutation null distribution.
- Select genes with empirical P < 0.05 as significantly associated.

5) Outputs & visualization

- Save significant genes with z-scores and P-values to CSV.
- Plot RMSE and explained-variance curves to inspect model selection and fit.

.. note::

    You may upload significant-gene lists to web tools such as Metascape for functional enrichment and pathway analysis: https://metascape.org/

.. code-block:: python

    run_transcriptome_pls_pipeline(
        t_values_file_path="/path/to/fc_deviation.csv",
        gene_expression_file_path="/path/to/transcriptome_contribution.csv",
        output_best_genes_path="/path/to/transcriptome_best_genes.csv",
        mse_plot_path="/path/to/transcriptome_pls_mse.tif",
        explained_variance_plot_path="/path/to/transcriptome_explained_variance.tif",
        max_components=10,
        cv_splits=5,
        n_permutations=1000,
    )

4. Cell-type enrichment of significant genes
-------------------------------------------

Test whether the set of significant genes is enriched in marker genes for specific cell types using permutation testing.

1) Data preparation

- Load the cell-type marker table containing marker genes for each cell class.
- Load the list of significant genes obtained from PLS-CGE.

.. note::

    Example cell-type marker tables used in the field are available from public repositories (e.g. Hansen et al.).

2) Observed overlap

Compute the size of the intersection between the target gene set and each cell-type marker set.

3) Permutation testing

- Sample random gene sets from the background marker gene pool and compute their overlaps with each cell-type marker set.
- Aggregate permutation overlaps to obtain null distributions per cell type.

4) Statistical evaluation

- Compute z-scores and two-sided empirical P-values from the permutation distributions.
- Apply multiple-testing correction (e.g. FDR-BH) across tested cell types.

.. code-block:: python

    gene_celltype_pipeline(
        celltype_path="/path/to/celltype_markers.csv",
        target_genes_path="/path/to/transcriptome_best_genes.csv",
        output_path="/path/to/celltype_enrichment.csv",
    )

5. Cortical layer enrichment
---------------------------

Assess whether the target gene set shows enrichment or depletion across cortical layers (e.g., layers 1–6) using permutation testing and z-score statistics.

1) Data preparation

- Read layer-marker statistics (e.g., t-statistics per gene for each layer) from the layer marker table.
- Intersect target genes with layer-marker genes to form the analysis set.

.. note::

    Example layer-marker resources come from spatial transcriptomics studies (see Supplementary Table references in the literature).

2) Observed statistics

Compute the mean layer-specific statistic (e.g., mean t) for the target gene set for each cortical layer.

3) Permutation null distributions

- Sample gene sets from the background pool (size matched to the target set) many times (e.g., 5000) and compute their layer means to build null distributions.

4) Z-scores, p-values and FDR

- For each layer compute the z-score and two-sided empirical p-value.
- Apply Benjamini–Hochberg FDR correction across layers.

5) Save results

- Save per-layer observed means, permutation means, z-scores, raw p-values and FDR-adjusted p-values to CSV for reporting.

.. code-block:: python

    gene_layer_pipeline(
        layer_marker_path="/path/to/layer_markers.xlsx",
        target_genes_path="/path/to/transcriptome_best_genes.csv",
        output_path="/path/to/layer_enrichment.csv",
    )

Notes
-----

- The snippet above is purely illustrative. To run these calls you must import
  the corresponding functions from the `BrainNetAnno` package or ensure the
  project's `src/` is on `PYTHONPATH`.
- Replace `/path/to/...` with real, accessible file paths and choose sensible
  parameters for `max_components`, `cv_splits`, and `n_permutations` based on
  your dataset size and computational resources.