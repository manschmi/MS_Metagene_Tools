#!/usr/bin/env python
'''Reads in a deeptools computeMatrix file and returns the matrix filter for all rows with sum below filterValue

Usage: filter_matrix.py 'matrix'.gz --filterType [exact,below,above,below_or_equal,above_or_equal] --filterValue --filterBed --filterString --regexMatch --filterColumn --outFileName


default mode: filters matrix by sum in each row using filterValue above, below or exact, ...
ie above: only return values greater to filterValues returned

alternative: filter entries in the bed ie --filterBed --filterString "SUT" --filterColumn 5 --regexMatch:
 removes all rows containing regex SUT in column 5 of the bed file (using python re.search() function )
--filterBed: filter using the bed interval information (no arguments)
    --filterString: a string to filter the bed intervals for (ie "SUT")
    --filterColumn: the column of the bed information to filter (ie typically, 4 or 5)
    --regexMatch: use regex matching with filterString instead of exact match

writes filtered matrix to file with name added _filtered before .gz, ie matrix_filtered.gz'''


__author__ = 'schmidm'

import sys
import gzip
import argparse
import json
import re


parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument('matrix_file')
parser.add_argument('--filterValue', default=0, type=float)
parser.add_argument('--filterBed', action="store_true")
parser.add_argument('--filterString', default='', type=str)
parser.add_argument('--filterColumn', default=0, type=int)
parser.add_argument('--regexMatch', default=False, action="store_true")
parser.add_argument('--outFileName')
parser.add_argument('--filterType', required=True, default='below', choices=['below', 'exact', 'above'])

args = parser.parse_args()


fname = args.matrix_file

if args.outFileName:
    outfile = args.outFileName
else:
    outfile = fname.replace('.gz', '_filtered.gz')

print 'filtering ', args.filterType, ' for value: ', args.filterValue
print 'writing filtered to ', outfile

def bounds_to_tuple(l):
    '''converts a list [0,10,20] to list of tuple [(0,10),(10,20)]'''
    return [(int(l[i])+1, int(l[i+1])) for i in range(len(l)-1)]


with gzip.open(fname, 'r') as f:
    meta = json.loads( f.readline().strip('@') )
    grp_bnds = bounds_to_tuple(meta['group_boundaries'])

    lnr = 1 #line number
    filt_cnt = 0
    filt_grp_bnds = [0]
    content = ''
    for grp in range(len(grp_bnds)):
        print 'filtering grp', meta['group_labels'][grp]
        max_sense_line = grp_bnds[grp][1]
        if lnr != grp_bnds[grp][0]:
            print 'something wrong with indexing sense'
            exit()
        while lnr <= max_sense_line:
            line = f.readline()
            le = line.rstrip().split('\t')
            if args.filterBed:
                if args.regexMatch:
                    if not re.search(args.filterString, le[args.filterColumn]):
                        content += line
                        filt_cnt += 1
                elif le[args.filterColumn] != args.filterString:
                    content += line
                    filt_cnt += 1
            else:
                line_sum = sum([ float(v) for v in le[7:] if v != 'nan' ])
                #print line_sum
                if args.filterType == 'below':
                    if line_sum >= args.filterValue:
                        content += line
                        filt_cnt += 1
                elif args.filterType == 'below_or_equal':
                    if line_sum > args.filterValue:
                        content += line
                        filt_cnt += 1
                elif args.filterType == 'exact':
                    if line_sum != args.filterValue:
                        content += line
                        filt_cnt += 1
                elif args.filterType == 'above':
                    if line_sum <= args.filterValue:
                        content += line
                        filt_cnt += 1
                elif args.filterType == 'above_or_equal':
                    if line_sum < args.filterValue:
                        content += line
                        filt_cnt += 1
            lnr += 1
        filt_grp_bnds.append(filt_cnt)
        print 'unfiltered: ', meta['group_boundaries'][grp+1] - meta['group_boundaries'][grp], 'after filter: ', filt_grp_bnds[grp+1] - filt_grp_bnds[grp]



meta['group_boundaries'] = filt_grp_bnds


with gzip.open(outfile, 'wb') as f:
    h=json.dumps(meta,separators=(',',':'))
    hx='@'+h+'\n'
    f.write(hx)
    f.write(content)
