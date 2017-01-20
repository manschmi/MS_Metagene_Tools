#!/usr/bin/env python
'''

Reads in a deeptools computeMatrix file and returns the matrix with trimmed range

Usage: python trim_matrix_range.py matrix.gz --upstream --downstream --outFileName

writes matrix to file with name added _trimmed before .gz, ie matrix_trimmed.gz

UPS: only works for reference-point matrices for now
'''

__author__ = 'schmidm'

import sys
import gzip
import json
import argparse


parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument('matrix_file')
parser.add_argument('--upstream', default=0, type=int)
parser.add_argument('--downstream', default=0, type=int)
parser.add_argument('--outFileName')

args = parser.parse_args()


fname = args.matrix_file

if args.outFileName:
    outfile = args.outFileName
else:
    outfile = fname.replace('.gz', '_trimmed.gz')

def bounds_to_tuple(l):
    '''converts a list [0,10,20] to list of tuple [(0,10),(10,20)]'''
    return [(int(l[i]), int(l[i+1])) for i in range(len(l)-1)]

with gzip.open(fname, 'r') as f:
    meta = json.loads(f.readline().strip('@'))
    sample_bnds = bounds_to_tuple(meta['sample_boundaries'])

    current_upstream_bins = meta['upstream'] / meta['bin size']
    if args.upstream < meta['upstream']:
        new_upstream_bins = args.upstream / meta['bin size']
        upstream_bin_sel = range( (current_upstream_bins-new_upstream_bins), current_upstream_bins )
        meta['upstream'] = args.upstream
    else:
        upstream_bin_sel = range(current_upstream_bins)

    #print 'upstream bins: ', upstream_bin_sel

    current_downstream_bins = meta['downstream'] / meta['bin size']
    first_downstream_bin = upstream_bin_sel[-1] + 1
    if args.downstream < meta['downstream']:
        new_downstream_bins = args.downstream / meta['bin size']
        downstream_bin_sel = range( first_downstream_bin, first_downstream_bin+new_downstream_bins )
        meta['downstream'] = args.downstream
    else:
        downstream_bin_sel = range( first_downstream_bin, first_downstream_bin+current_downstream_bins )


    #print 'downstream bins: ', downstream_bin_sel
    bin_sel = upstream_bin_sel + downstream_bin_sel
    #print 'all bins: ', bin_sel
    bin_cnt = len(bin_sel)
    new_sample_bounds = [ bin_cnt * i for i in range(len(meta['sample_labels'])+1)]

    content =''
    for line in f:
        le = line.rstrip().split('\t')
        le_samples = [ le[6+r[0]:6+r[1]] for r in sample_bnds]
        trimmed_samples = [ s[r] for s in le_samples for r in bin_sel ]
        content += '\t'.join(le[0:6])
        content += '\t'
        content += '\t'.join(trimmed_samples)
        content += '\n'


meta['sample_boundaries'] = new_sample_bounds

with gzip.open(outfile, 'wb') as f:
    # write output
    h = json.dumps(meta, separators=(',', ':'))
    hx = '@' + h + '\n'
    f.write(hx)
    f.write(content)
