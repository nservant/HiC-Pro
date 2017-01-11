================================================================================
Installation
================================================================================

This package uses distutils, which is the default way of installing
python modules.

The dependencies are:

- python (>= 2.6)
- setuptools
- numpy (>= 1.3)
- scipy (>= 0.7)
- argparse if python <2.7

In addition, pandas is recommanded for fast data loading.


All of these dependencies can be installed at once using `Anaconda
<http://docs.continuum.io/anaconda/install.html>`_

To install in your home directory, use::

    python setup.py install --user

or using pip::

    pip install --user iced

