Usage
=====

Quickstart (Python API)
-----------------------

Transcriptome PLS-CGE pipeline::

   from BrainNetAnno.transcriptome_pls_cge import run_transcriptome_pls_pipeline

   run_transcriptome_pls_pipeline(
       t_values_file="path/to/fc_deviation.csv",
       gene_expression_file="path/to/gene_contrib.csv",
       output_best_genes_csv="path/to/best_genes.csv",
       fig_outputfile="path/to/rmse_variance.tif",
   )

Cell-type enrichment::

   from BrainNetAnno.gene_celltype import run_pipeline

   res = run_pipeline(
       celltype_path="path/to/celltypes_PSP.csv",
       target_genes_path="path/to/best_genes_SCZ.csv",
       output_path="path/to/celltype_enrichment_results.csv",
   )
   print(res.head())

Cortical layer enrichment::

   from BrainNetAnno.gene_layer import run_pipeline

   res = run_pipeline(
       layer_marker_path="path/to/41593_2020_787_MOESM3_ESM.xlsx",
       target_genes_path="path/to/best_genes_all.csv",
       output_csv="path/to/ALL_gene_layer_analysis_results.csv",
   )
   print(res.head())
