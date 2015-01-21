## HiC-Pro
## Copyleft 2015 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007. 

##
## Check Python lib import
##

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


## 
try:
    import scipy
except ImportError:
    raise ImportError('Error - scipy cannot be imported')


## numpy
try:
    import argparse
except ImportError:
    raise ImportError('Error - argparse cannot be imported')


