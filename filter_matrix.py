#!/usr/bin/env python
'''Reads in a deeptools computeMatrix file and returns the matrix counting all positions with a value (ie all not nan or 0) as 1 all others as 0

Usage: filter_matrix.py 'matrix'.gz --minSum --outFileName




writes filtered matrix to file with name added _filtered before .gz, ie matrix_filtered.gz'''


__author__ = 'schmidm'

import sys
import gzip
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('matrix_file')
parser.add_argument('-h', action="store_true", default=False)

parser.add_argument('--minSum', default=0)
parser.add_argument('--outFileName')
#parser.add_argument('--filterRegion', default='all')
#parser.add_argument('--filterSamples', default='all')

args = parser.parse_args()


if args.h:
    print __doc__
    exit()

fname = args.matrix_file

if args.outFileName:
    outfile = args.outFileName
else:
    outfile = fname.replace('.gz', '_filtered.gz')

def bounds_to_tuple(l):
    '''converts a list [0,10,20] to list of tuple [(0,10),(10,20)]'''
    return [(int(l[i])+1, int(l[i+1])) for i in range(len(l)-1)]


with gzip.open(fname, 'r') as f:
    meta = json.loads( f.readline().strip('@') )
    grp_bnds = bounds_to_tuple(meta['group_boundaries'])

    lnr = 1 #line number
    filt_cnt = 0
    for grp in range(len(grp_bnds)):
        filt_grp_bnds.append((filt_cnt, filt_cnt))
        max_sense_line = grp_bnds[grp][1]
        if lnr != grp_bnds[grp][0]:
            print 'something wrong with indexing sense'
            exit()
        while lnr <= max_sense_line:
            line = f.readline()
            le = line.rstrip().split('\t')
            line_sum = sum([ float(v) for i, v in enumerate(le[6:]) if v != 'nan' ])
            if line_sum >= args.minSum:
                content += line
                filt_cnt += 1
            lnr += 1
        filt_grp_bnds[grp][1] = filt_cnt


meta['group_boundaries'] = filt_grp_bnds


with gzip.open(outfile, 'wb') as f:
    h=json.dumps(meta,separators=(',',':'))
    hx='@'+h+'\n'
    f.write(hx)
    f.write(content)
