Mitochondrial contribution example
==================================

This example computes mitochondrial per-connection contribution scores using
regional coordinates and mitochondrial expression matrices.

Run the script::

   python tests/examples/example_mitochondrial.py

Core snippet
------------

.. code-block:: python

   from BrainNetAnno.mitochondrial import run_mitochondrial_pipeline

   params, contrib_df = run_mitochondrial_pipeline(
       coordinates_path=coords_csv,
       mitochondrial_expression_path=mito_xlsx,
       output_contribution_path=mito_contrib_csv,
       initial_params=(1.0, 50.0, 0.0),
       save_plot=False,
       plot_path=None,
   )

Outputs
-------
- Contribution CSV with `Region_Pair` and per-feature columns.
