#!/usr/bin/env python
'''
Reads in a deeptools computeMatrix file and returns the minimum value in the matrix

Usage: min_value_in_matrix.py 'matrix'.gz greaterX=0

print the minimum value in the matrix, if argument greaterX is added, only min value greater X will be printed

'''


__author__ = 'schmidm'

import sys
import gzip
from math import log

if len(sys.argv) <= 2 or sys.argv[1] == '-h':
    print __doc__
    exit()


fname = sys.argv[1]

if (len(sys.argv) >= 3):
    floor = float( sys.argv[2].split('=')[1] )
else:
    floor = -float('inf')

min_val = float('inf')

with gzip.open(fname, 'r') as f:
    for line in f:
        if line[0] == '@':
            continue
        le = line.rstrip().split('\t')
        for v in le[6:]:
            v = float(v)
            if v > floor and v < min_val:
                min_val = v


print min_val
