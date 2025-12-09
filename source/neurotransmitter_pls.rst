Neurotransmitter PLS example
============================

This example links FC weights with per-connection neurotransmitter contributions
using PLS with CV-based component selection and permutation testing of weights.

Run the example script::

   python tests/examples/example_neurotransmitter_pls.py

Core snippet
------------

.. code-block:: python

   from BrainNetAnno.neurotransmitter_pls_cge import run_neurotransmitter_pls_pipeline

   best_n, df = run_neurotransmitter_pls_pipeline(
       fc_matrix_path=fc_csv,
       nt_contrib_path=nt_contrib_csv,
       top_k=350,
       mse_plot_path=mse_plot,
       explained_variance_plot_path=var_plot,
       output_weights_path=weights_csv,
       max_components=10,
       cv_splits=5,
       random_state=42,
       n_permutations=1000,
   )

Notes
-----
- Tune `top_k` to limit the connections used for PLS (e.g., top-K by |FC|).
- Reduce `n_permutations` for faster runs during development.
