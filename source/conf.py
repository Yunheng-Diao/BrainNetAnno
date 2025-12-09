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
    'sphinx.ext.autosummary',
]

# If you want to generate API docs automatically, you can enable sphinx.ext.autosummary
autosummary_generate = True
# Optional: include imported members in summaries (can be noisy)
# autosummary_imported_members = True

templates_path = ['_templates']
exclude_patterns = []

# Use Read the Docs theme if available
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
# Show logo on all pages
html_logo = '_static/logo.png'
html_theme_options = {
    'logo_only': True,
    'display_version': False,
}
# # Load custom CSS for styling (e.g., logo size)
# html_css_files = ['custom.css']

# Autodoc settings
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

# Mock heavy scientific deps on RTD to prevent build failures
autodoc_mock_imports = [
    'numpy', 'pandas', 'scipy', 'sklearn', 'scikit-learn', 'statsmodels', 'matplotlib', 'torch'
]

# Napoleon (Google/NumPy style docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = True
