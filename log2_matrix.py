#!/usr/bin/env python
'''

Reads in a deeptools computeMatrix file and returns the matrix with log2 values for all positions with a value > 0, all other will be set to nan

Usage: log2_matrix.py matrix.gz outfilename pseudocount

writes log2 matrix to file with name added _log2 before .gz, ie matrix_log2.gz
pseudocount is added to all values of matrix and replaces 0 and nan values.
use min_value_in_matrix.py to get the min value for example.

'''

__author__ = 'schmidm'

import sys
import gzip
import math

if len(sys.argv) < 2 or sys.argv[1] == '-h':
    print __doc__
    exit()

fname = sys.argv[1]

if len(sys.argv) >= 3:
    outfile = sys.argv[2]
else:
    outfile = fname.replace('.gz', '_log2.gz')

if len(sys.argv) >= 4:
    pseudocount = float(sys.argv[3])
    str_log2_pseudocount = str(math.log(pseudocount, 2))
    print 'adding pseudocount ', pseudocount
else:
    pseudocount = 0
    print 'computing log2 without adding pseudocount'


with gzip.open(fname, 'r') as f:
    for line in f:
        if line[0] == '@':
            content = line
            #print line.rstrip()
        else:
            le = line.rstrip().split('\t')
            line = [ v if (i < 6) else str(math.log((float(v)+pseudocount), 2)) if (v != 'nan' and v > 0) else str_log2_pseudocount for i,v in enumerate(le) ]
            #line = [ v if (i < 6) else str(math.log(float(v), 2)) if (v != 'nan' and float(v) > 0) else 'nan' for i, v in enumerate(le) ]

            content += '\t'.join(line)
            content += '\n'



with gzip.open(outfile, 'wb') as f:
    f.write(content)
