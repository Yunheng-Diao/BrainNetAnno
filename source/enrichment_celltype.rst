Cell-type enrichment example
===========================

This example demonstrates permutation-based cell-type enrichment testing for a
set of target genes against marker gene sets.

Run the example script::

   python tests/examples/example_enrichment_celltype.py

Core snippet
------------

.. code-block:: python

   from BrainNetAnno.gene_celltype import run_pipeline as gene_celltype_pipeline

   res = gene_celltype_pipeline(
       celltype_path=celltype_csv,
       target_genes_path=target_genes_csv,
       output_path=output_csv,
   )

Outputs
-------
- CSV of enrichment results with columns: Cell Type, Overlap, Z-Score, Raw P-Value,
  Mean Overlap (Random), FDR
