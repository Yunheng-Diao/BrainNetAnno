Transcriptome PLS-CGE example
=============================

This example shows a complete PLS-CGE analysis that links FC deviations with per-connection
gene contributions. It demonstrates cross-validated component selection, permutation
significance testing, and saving of results and plots.

Usage
-----

Run the example script shipped in the repository::

   python tests/examples/example_transcriptome_pls.py

Key parameters used in the example
---------------------------------
- `t_values_file_path`: FC deviations CSV (square, no header)
- `gene_expression_file_path`: per-connection gene contribution CSV (must contain `Region_Pair`)
- `output_best_genes_path`: where significant genes are saved
- `mse_plot_path`, `explained_variance_plot_path`: optional plots

Core snippet
------------

.. code-block:: python

   from BrainNetAnno.transcriptome_pls_cge import run_transcriptome_pls_pipeline

   run_transcriptome_pls_pipeline(
       t_values_file_path=fc_csv,
       gene_expression_file_path=contrib_csv,
       output_best_genes_path=best_genes_csv,
       mse_plot_path=mse_plot,
       explained_variance_plot_path=var_plot,
       max_components=10,
       cv_splits=5,
       random_state=42,
       n_permutations=1000,
   )

Notes
-----
- Replace the example input paths with your dataset paths. The example assumes
  supporting example data in `examples_data/` and writes outputs to `examples_output/`.
- For quick testing reduce `n_permutations` and `max_components`.
