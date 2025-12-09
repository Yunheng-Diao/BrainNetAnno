Mitochondrial PLS example
=========================

This example performs PLS analysis between FC and mitochondrial contribution
features, including CV-based component selection and permutation testing.

Run the example script::

   python tests/examples/example_mitochondrial_pls.py

Core snippet
------------

.. code-block:: python

   from BrainNetAnno.mitochondrial_pls_cge import run_mitochondrial_pls_pipeline

   best_n, df = run_mitochondrial_pls_pipeline(
       fc_matrix_path=fc_csv,
       nt_contrib_path=mito_contrib_csv,
       top_k=500,
       mse_plot_path=mse_plot,
       explained_variance_plot_path=var_plot,
       output_weights_path=weights_csv,
       max_components=6,
       cv_splits=5,
       random_state=42,
       n_permutations=1000,
   )
