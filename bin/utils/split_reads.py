#! /usr/bin/env python

import argparse
import os
import subprocess
import sys
from glob import glob
import shutil
import re

parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--results_folder", "-r", default="./")
parser.add_argument("--nreads", "-n", default=20000000)
args = parser.parse_args()

filename = args.filename
nlines = int(args.nreads)*4
out = args.results_folder
try:
    os.makedirs(out)
except OSError:
    pass


if filename.endswith('.gz'):
    prefix = re.sub('((.fastq)|(.fq)).gz','_part', os.path.join(out, os.path.basename(filename)))
    cmd = "zcat -fc {} | split -l {} -d - {}".format(
        filename, nlines, prefix)
else:
    prefix = re.sub('(.fastq)|(.fq)','_part', os.path.join(out, os.path.basename(filename)))
    cmd = "split -l {} -d {} {}".format(
        nlines, filename, prefix)

retcode = subprocess.call(cmd, shell=True)
if retcode != 0:
    print("split file failed with return code {}".format(retcode))
    sys.exit(1)

files = glob(os.path.join(prefix + "*"))
files.sort()

for ifile in files:
    shutil.move(ifile,
                os.path.join(os.path.dirname(ifile),
                             os.path.basename(ifile)[-2:] + "_" +
                             os.path.basename(ifile)[:-7] + ".fastq"))
