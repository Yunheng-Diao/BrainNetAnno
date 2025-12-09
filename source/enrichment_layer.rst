Cortical layer enrichment example
================================

This example shows how to test whether a target gene set shows elevated layer-
specific signals using layer marker statistics.

Run the example script::

   python tests/examples/example_enrichment_layer.py

Core snippet
------------

.. code-block:: python

   from BrainNetAnno.gene_layer import run_pipeline as gene_layer_pipeline

   res = gene_layer_pipeline(
       layer_marker_path=layer_xlsx,
       target_genes_path=target_genes_csv,
       output_path=output_csv,
   )

Outputs
-------
- CSV of per-layer enrichment statistics (Z, p-values, FDR) and mean t-statistics.
