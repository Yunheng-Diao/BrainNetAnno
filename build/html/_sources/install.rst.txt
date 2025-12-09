Installation and setup
======================

.. note::

   Using a virtual environment (``venv`` or ``conda``) is highly recommended to
   avoid conflicts with system packages. The examples below assume you are
   working in an isolated environment.


Basic installation
------------------

This package requires Python 3.8 or later. If you already have a supported
Python version installed, the simplest option is to install the latest
released wheel from PyPI:

.. code-block:: bash

   pip install BrainNetAnno

This will install the package and its Python-only dependencies listed in the
project metadata. The wheel distribution avoids source compilation and is the
recommended route for most users.


Install from source (GitHub)
----------------------------

If you prefer to install the most up-to-date version from the GitHub
repository (for example to test a fix or new feature), clone the repository
and install locally:

.. code-block:: bash

   git clone https://github.com/Yunheng-Diao/BrainNetAnno.git
   cd BrainNetAnno
   python -m pip install -U build
   python -m build
   pip install dist/BrainNetAnno-<version>-py3-none-any.whl

Or install editable/develop mode (useful for development):

.. code-block:: bash

   git clone https://github.com/Yunheng-Diao/BrainNetAnno.git
   cd BrainNetAnno
   pip install -e .

.. important::

   When installing directly from the repository, make sure you are following
   the most up-to-date documentation. Installing from source may expose you
   to untested changes.


Optional I/O extras
-------------------

Some features (e.g., faster parquet-backed I/O or optional dependency-driven
converters) require additional packages. These extras are optional and can be
installed with ``pip`` using extras syntax:

.. code-block:: bash

   pip install "BrainNetAnno[io]"

Install the extras from GitHub as follows:

.. code-block:: bash

   git clone https://github.com/Yunheng-Diao/BrainNetAnno.git
   cd BrainNetAnno
   pip install .[io]


Platform notes (Linux)
----------------------

On some Linux distributions you may need to install system-level libraries
before the Python wheels for optional packages can be built. For example,
to enable Snappy/parquet support you may need the following packages:

- Debian/Ubuntu:

  .. code-block:: bash

     sudo apt-get update
     sudo apt-get install -y libsnappy-dev

- Fedora/CentOS:

  .. code-block:: bash

     sudo dnf install -y snappy-devel

If you encounter compilation issues when installing optional dependencies,
consider using prebuilt binary wheels (many are available on PyPI) or install
via Conda where binary packages are commonly provided:

.. code-block:: bash

   conda create -n brainnetanno python=3.10
   conda activate brainnetanno
   pip install BrainNetAnno


Troubleshooting
---------------

- If ``pip install`` fails while building a package wheel, inspect the
  console output for missing system packages (often headers like ``libxml2``
  or ``libsodium``). Install the required OS packages and retry.
- If SciPy or NumPy compile from source and take a long time, prefer installing
  from wheels: ``pip install numpy scipy`` before installing BrainNetAnno.
- Use a fresh virtual environment to rule out dependency conflicts.


Developer / editable install
----------------------------

For contributors who want to develop the package locally, install in editable
mode and run the tests locally:

.. code-block:: bash

   git clone https://github.com/Yunheng-Diao/BrainNetAnno.git
   cd BrainNetAnno
   pip install -e .[dev]

Run tests (requires test dependencies):

.. code-block:: bash

   pytest -q


Local docs build (optional)
---------------------------

To preview the Sphinx documentation locally (useful after editing the
install page), install the docs requirements and build HTML:

.. code-block:: bash

   pip install -r source/requirements.txt
   sphinx-build -b html source source/_build/html

Open ``source/_build/html/index.html`` in your browser to preview the docs.


If you still need help, open an issue on the GitHub repository with a
description of your platform and the error messages you encountered.

