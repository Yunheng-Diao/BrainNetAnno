import os
import sys
sys.path.insert(0, os.path.abspath('..'))
# Ensure local source tree (src/) is importable for autodoc
sys.path.insert(0, os.path.abspath('../src'))

project = 'BrainNetAnno'
copyright = '2025, BrainNetAnno'
author = 'BrainNetAnno contributors'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # for Google/NumPy docstrings
    'sphinx.ext.viewcode',
]

# If you want to generate API docs automatically, you can enable sphinx.ext.autosummary
# extensions.append('sphinx.ext.autosummary')
# autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = []

# Use Read the Docs theme if available
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Autodoc settings
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}
