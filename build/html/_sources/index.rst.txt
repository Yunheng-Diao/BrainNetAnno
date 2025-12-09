.. BrainNetAnno documentation master file, created by
   sphinx-quickstart on Mon Dec  8 20:39:25 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

BrainNetAnno: A toolkit for molecular annotation of brain networks.
==========================

.. image:: _static/logo.png
   :alt: BrainNetAnno
   :align: center
   :width: 220px

With the rapid advancement of multi-omics technologies and large-scale brain atlas initiatives, molecular-level information about the human brain has become increasingly abundant, including gene expression, neurotransmitter receptor densities, cell-type distributions, mitochondrial functional signatures, synaptic mechanisms, and metabolic features. These resources offer important opportunities to elucidate the biological foundations of brain development and neuropsychiatric disorders. However, brain function is not carried out by isolated regions in isolation; rather, it relies on a highly interconnected network architecture. Inter-regional connections work together to support cognition, emotion, and behavior, and disease-related abnormalities often manifest as systematic alterations in network connectivity. Therefore, understanding the molecular basis of network-level abnormalities is critical for uncovering the pathophysiology of neuropsychiatric disorders.

Although substantial molecular annotation resources now exist at the level of brain regions, enabling characterization of regional molecular properties, methods that focus on the molecular annotation of the connections themselves remain highly limited. In other words, while we know the molecular attributes of individual regions, we still lack answers to fundamental questions such as why certain connections are more vulnerable, what molecular mechanisms drive connection-level abnormalities, or what underlies the susceptibility of specific connections in disease. This gap constrains the integration of multi-omics data with network neuroscience and limits our ability to understand complex brain disorders from a systems-level perspective.

The present toolkit was developed to fill this critical gap. It provides systematic molecular annotation of inter-regional connections, enabling researchers to directly link network connectivity alterations to underlying gene-expression profiles, neurotransmitter systems, mitochondrial features, cell-type compositions, synaptic mechanisms, and metabolic characteristics. By integrating multi-omics information at the level of connections, the toolkit offers a new pathway for understanding the pathophysiology of neuropsychiatric disorders and provides a biologically informed, mechanistically interpretable framework for identifying targets for network-based interventions such as brain stimulation.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   install
   usage
   Example
   api
   changelog
   modules
