#!/usr/bin/env python
'''Reads in a deeptools computeMatrix file and returns the matrix counting all positions with a value (ie all not nan or 0) as 1 all others as 0

Usage: filter_matrix.py 'matrix'.gz --filterCount --filterRegion --filterSamples

writes filtered matrix to file with name added _filtered before .gz, ie matrix_filtered.gz'''


__author__ = 'schmidm'

import sys
import gzip
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('matrix_file')
parser.add_argument('-h', action="store_true", default=False)

parser.add_argument('--filterMinCount', default=0)
parser.add_argument('--filterRegion', default='all')
parser.add_argument('--filterSamples', default='all')

args = parser.parse_args()


if args.h:
    print __doc__
    exit()

fname = args.matrix_file


with gzip.open(fname, 'r') as f:
    for line in f:
        if line[0] == '@':
            content = line
            #print line.rstrip()
        else:
            le = line.rstrip().split('\t')

            line = [ v if (i < 6) else '1' if (v != 'nan' and float(v) > 0) else '0' for i, v in enumerate(le) ]
            #print '\t'.join(line)
            content += '\t'.join(line)
            content += '\n'

outfile = fname.replace('.gz', '_events.gz')

with gzip.open(outfile, 'wb') as f:
    f.write(content)
