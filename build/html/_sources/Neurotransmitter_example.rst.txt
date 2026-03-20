Neurotransmitter Example
================

A complete pipeline for neurotransmitter analysis is demonstrated in this example.

1. Neurotransmitter feature importance via LOFO
----------------------------------------------

Estimate the importance of each neurotransmitter feature using a Leave-One-Feature-Out (LOFO) strategy. The approach:

1) Map the neurotransmitter profile onto the standard brain atlas

.. code-block:: python

    from BrainNetAnno.molecular_mapping import MolecularBrainMapper
    # 1. Initialize the mapper
    mapper = MolecularBrainMapper(
        atlas_path="/path/to/aal116MNI.nii.gz", # Replace with your brain atlas path
        data_dir="/path/to/Neurotransmitter_Profile",  # Replace with your data path
        standard_template='MNI152',
        verbose=True
    )

    # 2. Load mitochondrial data
    mito_data = mapper.load_mitochondrial_data(file_pattern="*.nii.gz")

    # 3. Map to brain atlas
    result_df = mapper.map_to_atlas(
        mito_data,
        missing='centroids',    # Strategy for handling missing values
        statistic='mean'        # Extract mean value
    )

    # 4. Preview results
    print("\nPreview of mapping results:")
    print(result_df.head())

    # 5. Save results
    mapper.save_results(result_df, "/path/to/neurotransmitter_mapping_results.csv")

.. note::

    Neurotransmitter profile data can be obtained from public repositories (for example, https://github.com/juryxy/JuSpace/tree/JuSpace_v2.0/JuSpace_v2/PETatlas).

2) Data loading & preprocessing

- Load atlas MNI coordinates and the neurotransmitter distribution matrix (regions × features).
- Compute pairwise Euclidean distances between regions.
- Row-wise z-score normalize neurotransmitter features.


3) Fit baseline network model

- Use a Ledoit–Wolf shrinkage estimator to obtain a stable covariance estimate and invert it to obtain a precision matrix.
- Extract the upper-triangle elements of the precision matrix as the observed network vector.

4) Fit distance–connection expectation

- Fit an exponential decay model (as above) to model expected connection strength as a function of distance.

5) Compute residual network

- Compute the residuals by subtracting the expected connection matrix from the observed precision-based connections.

6) LOFO iterations

For each neurotransmitter feature k:

- Remove feature k from the neurotransmitter matrix and re-fit the covariance/precision using the remaining features.
- If the reduced data causes singular covariance, record and skip that feature.
- Compute the new residual network and quantify the change in residual magnitude per connection:

    .. math::

        Contribution_{i,j}^{(k)} = |Residual_{i,j}^{(full)}| - |Residual_{i,j}^{(LOFO_k)}|

- Store contributions in a contributions matrix (connections × features).

7) Save results

- Save the neurotransmitter contribution matrix (connections × features) to CSV for downstream PLS and visualization.

.. code-block:: python

    run_neurotransmitter_pipeline(
        coordinates_path="/path/to/coords.csv",
        neurotransmitter_expression_path="/path/to/neurotransmitter_mapping_results.csv",
        output_contribution_path="/path/to/neurotransmitter_contribution.csv",
        initial_params=(1.0, 50.0, 0.0),
        save_plot=False,
        plot_path=None,
    )

2. Neurotransmitter PLS-CGE association
----------------------------------------

Use PLS regression to relate FC deviations to neurotransmitter-derived connection features. Steps mirror the transcriptome PLS workflow:

1) Prepare features and target vector (select top-N connections if desired).
2) Use K-fold CV to select the optimal number of PLS components and plot MSE curve.
3) Train final PLS on full dataset and compute feature weights.
4) Use permutation testing to derive empirical p-values for feature weights.
5) Save weights, explained variance and p-values to CSV.

.. code-block:: python

    run_neurotransmitter_pls_pipeline(
        fc_matrix_path="/path/to/fc_deviation.csv",
        nt_contrib_path="/path/to/neurotransmitter_contribution.csv",
        output_weights_path="/path/to/neurotransmitter_weights.csv",
        max_components=10,
        cv_splits=5,
        n_permutations=1000,
    )

Notes
-----

- The snippet above is purely illustrative. To run these calls you must import
  the corresponding functions from the `BrainNetAnno` package or ensure the
  project's `src/` is on `PYTHONPATH`.
- Replace `/path/to/...` with real, accessible file paths and choose sensible
  parameters for `max_components`, `cv_splits`, and `n_permutations` based on
  your dataset size and computational resources.