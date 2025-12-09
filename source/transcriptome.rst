Transcriptome contribution example
=================================

This example computes per-connection gene contribution scores given region
coordinates and a regional gene-expression matrix.

Run the script::

   python tests/examples/example_transcriptome.py

Core snippet
------------

.. code-block:: python

   from BrainNetAnno.transcriptome import run_transcriptome_pipeline

   params, contrib_df = run_transcriptome_pipeline(
       coordinates_path=coords_csv,
       gene_expression_path=genes_csv,
       output_contribution_path=output_contrib_csv,
       initial_params=(0.64, 90.4, -0.19),
       save_plot=False,
       plot_path=None,
   )

Outputs
-------
- Contribution CSV with `Region_Pair` and per-gene columns.
- Optional fit figure if `save_plot=True`.
