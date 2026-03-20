Mitochondrial Example
================

A complete pipeline for mitochondrial analysis is demonstrated in this example.


1. Mitochondrial LOFO and contribution analysis
----------------------------------------------

Apply the LOFO strategy to mitochondrial-related phenotypes to quantify each phenotype's contribution to network residuals. Steps are analogous to the neurotransmitter LOFO workflow:

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

    mitochondrial profile data can be obtained from public repositories (for example, http://humanmitobrainmap.bcblab.com).

2) Data & spatial preparation: load MNI coordinates and mitochondrial expression matrix; compute distances and z-score normalize.
3) Fit baseline network with Ledoit–Wolf and compute residuals relative to distance expectation.
4) For each mitochondrial phenotype, remove it, re-fit, compute new residuals and record contribution changes.
5) Aggregate contributions into a matrix (connections × mitochondrial phenotypes) and save to CSV.

.. code-block:: python

    run_mitochondrial_pipeline(
        coordinates_path="/path/to/coords.csv",
        mitochondrial_expression_path="/path/to/mitochondrial_expression.xlsx",
        output_contribution_path="/path/to/mitochondrial_contribution.csv",
        initial_params=(1.0, 50.0, 0.0),
        save_plot=False,
        plot_path=None,
    )

2. Mitochondrial PLS-CGE association
-----------------------------------

Relate FC deviations to mitochondrial phenotype contributions using PLS-CGE. Main steps:

1) Assemble feature matrix (connections × mitochondrial phenotypes) and target vector.
2) Optionally select top-N connections by absolute effect size.
3) Use K-fold CV to choose the optimal number of PLS components and visualize MSE curves.
4) Fit final PLS, extract weights and compute explained variance.
5) Use permutations to derive empirical P-values for weights and save results.

.. code-block:: python

    run_mitochondrial_pls_pipeline(
        fc_matrix_path="/path/to/fc_deviation.csv",
        nt_contrib_path="/path/to/mitochondrial_contribution.csv",
        output_weights_path="/path/to/mitochondrial_weights.csv",
        max_components=6,
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