## Nicolas Servant
##
## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details


## bx
try:
    import bx.intervals
except ImportError:
    raise ImportError('Error - bx-python cannot be imported')

## numpy
try:
    import numpy
except ImportError:
    raise ImportError('Error - numpy cannot be imported')


## scipy
try:
    import scipy
except ImportError:
    raise ImportError('Error - scipy cannot be imported')


## argparse
try:
    import argparse
except ImportError:
    raise ImportError('Error - argparse cannot be imported')

## pysam
try:
    import pysam
except ImportError:
    raise ImportError('Error - pysam cannot be imported')

