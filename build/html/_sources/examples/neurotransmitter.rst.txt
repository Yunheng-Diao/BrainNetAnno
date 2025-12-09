Neurotransmitter contribution example
====================================

This example computes per-connection neurotransmitter contribution scores from
region coordinates and a transmitter/receptor expression matrix.

Run the script::

   python tests/examples/example_neurotransmitter.py

Core snippet
------------

.. code-block:: python

   from BrainNetAnno.neurotransmitter import run_neurotransmitter_pipeline

   params, contrib_df = run_neurotransmitter_pipeline(
       coordinates_path=coords_csv,
       neurotransmitter_expression_path=nt_xlsx,
       output_contribution_path=nt_contrib_csv,
       initial_params=(1.0, 50.0, 0.0),
       save_plot=False,
       plot_path=None,
   )

Outputs
-------
- Contribution CSV with `Region_Pair` and per-feature columns.
