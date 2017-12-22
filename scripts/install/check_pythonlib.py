## Nicolas Servant
##
## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

import re

## Define cmp() for installation with python >= 3.0.0
def cmp(a, b):
    return (a > b) - (a < b)

def vcmp(version1, version2):
    def normalize(v):
        ## 0.18.0.dev0+609facc
        v = re.sub(r'((\.)?dev).*$','', v)
        ## 0.8.4pre0
        v = re.sub(r'((\.)?pre).*$','', v)
        return [int(x) for x in re.sub(r'(\.0+).*$','', v).split(".")]
    return cmp(normalize(version1), normalize(version2))


## bx
try:
    import bx.intervals
    if vcmp(bx.__version__, '0.5.0') == -1:
        raise ValueError('bx-python '+bx.__version__+' detected. Version >= 0.5.0 required')
except ImportError:
    raise ImportError('Error - bx-python cannot be imported')  


## numpy
try:
    import numpy
    if vcmp(numpy.__version__, '1.8.2') == -1:
        raise ValueError('numpy '+numpy.__version__+' detected. Version >= 1.8.2 required')
except ImportError:
    raise ImportError('Error - numpy cannot be imported')


## scipy
try:
    import scipy
    if vcmp(scipy.version.version, '0.15.1') == -1:
        raise ValueError('scipy '+scipy.version.version+' detected. Version >= 0.15.1 required')
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
    if vcmp(pysam.__version__, '0.8.3') == -1:
        raise ValueError('pysam '+pysam.__version__+' detected. Version >= 0.8.3 required')
except ImportError:
    raise ImportError('Error - pysam cannot be imported')

