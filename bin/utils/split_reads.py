#! /usr/bin/env python

import argparse
import os
import subprocess
import sys
from glob import glob
import shutil


parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--results_folder", "-r")
parser.add_argument("--nreads", "-n", default=10000000)
args = parser.parse_args()

filename = args.filename
nlines = int(args.nreads)*4
out = args.results_folder
try:
    os.makedirs(out)
except OSError:
    pass

cmd = "split -l %d -d %s %s" % (
    nlines, filename,
    os.path.join(out, os.path.basename(filename)))

retcode = subprocess.call(cmd.split())
if retcode != 0:
    print "split file failed with return code %d", retcode
    sys.exit(1)

files = glob(os.path.join(out, os.path.basename(filename) + "*"))
files.sort()
for ifile in files:
    shutil.move(ifile,
                os.path.join(os.path.dirname(ifile),
                             os.path.basename(ifile)[-2:] + "_" +
                             os.path.basename(ifile)[:-2]))
