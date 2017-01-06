#!/usr/bin/env python
'''

Reads in a deeptools computeMatrix file and sorts it according to the order in a provided bed file.


Only makes sense if the matrix and bed files contain the same regions, or if bed is a subset of regions in matrix.
Any regions present in bed but not in matrix will be dropped with a warning, specifying the region.
Regions present in matrix but not bed will be counted and dropped with a warning but not specified.

Groupings in the matrix are lost in this operation!!

Groupings can be derived from the bedfile using a specific column from the bedfile
ie bed line: chr1  100  10000  ENSG001231232   protein-coding  +
can be split using --groupByColumn 4
default is not to group the output

Usage: sort_matrix_using_bed.py matrix.gz bed.bed matrix_bedsorted.gz [--groupByColumn]

Example: python /Users/schmidm/Documents/MS_Metagene_Tools/sort_matrix_using_bed.py matrix.gz order.bed -o matrix_bedsorted.gz

writes to matrix_bedsorted.gz if not provided it concatenated the names and uses suffix _bedsorted.gz


!!!UPS: groupByColumn currently requires the bedfile to be sorted by that column!
'''


__author__ = 'schmidm'

import sys
import gzip
import json
import argparse

parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument('matrix_file')
parser.add_argument('bed_file')
parser.add_argument('--groupByColumn', default=None, type=int)
parser.add_argument('--outFileName')


args = parser.parse_args()


matrix = args.matrix_file
bed = args.bed_file

if args.outFileName:
    outfile = args.outFileName
else:
    outfile = matrix.replace('.gz', '_bedsorted.gz')

print 'writing resorted to ', outfile



if args.groupByColumn:
    group_column = args.groupByColumn - 1
    print 'grouping by column ', args.groupByColumn

def bounds_to_tuple(l):
    '''converts a list [0,10,20] to list of tuple [(0,10),(10,20)]'''
    return [(int(l[i])+1, int(l[i+1])) for i in range(len(l)-1)]


def bed_lines(fname):
    '''iterator through lines of a bed file'''
    with open(fname, 'r') as bed:
        for line in bed:
            if line[0] == "#":
                continue
            yield line.rstrip().split('\t')


#open sense and antisense file and iterate through an merge the contents
with gzip.open(matrix, 'r') as f:

    #parse the headers
    meta = json.loads( f.readline().strip('@') )

    #sample boundaries sense and antisense
    sample_bnds = bounds_to_tuple(meta['sample_boundaries'])

    ##loop through the bed files and find matching row of matrix for each

    matrix = [line.split('\t') for line in f]
    matrix_region_cnt = len(matrix)


sorted_content = []
bed_region_cnt=0
found_region_cnt=0
current_group_label=''
current_group_start=0
group_labels=[]
group_boundaries=[0]
#bed_regions_found=0

for region in bed_lines(bed):
    bed_region_cnt += 1

    if args.groupByColumn and current_group_label != region[group_column]:
        if current_group_label != '' and group_boundaries[-1] != bed_region_cnt :
            print 'added group: ', current_group_label, ' at:', current_group_start, ' to ', bed_region_cnt
            group_labels.append(current_group_label)
            group_boundaries.append(bed_region_cnt-1)
        current_group_label = region[group_column]

    found = False
    for line in matrix:
        if line[0] == region[0] and line[1] == region[1] and line[2] == region[2] and line[5] == region[5]:
            found_region_cnt += 1
            found = True
            sorted_content.append('\t'.join(line))
            break

    if not found:
        print 'WARNING region from bed file not found in matrix: ', region


if args.groupByColumn:
    if current_group_label != '' and group_boundaries[-1] != found_region_cnt:
        print 'added group: ', current_group_label, ' at:', current_group_start, ' to ', bed_region_cnt
        group_labels.append(current_group_label)
        group_boundaries.append(bed_region_cnt)
else:
    group_boundaries = [0,found_region_cnt]
    group_labels = ['genes']

meta['group_boundaries'] = group_boundaries
meta['group_labels'] = group_labels


print 'bed file contained: ', bed_region_cnt, 'regions'
print 'matrix file contained: ', matrix_region_cnt, 'regions'
print 'found ', found_region_cnt, ' of the bed regions in the matrix'

#write output
with gzip.open(outfile, 'wb') as f:
    h=json.dumps(meta,separators=(',',':'))
    hx='@'+h+'\n'
    f.write(hx)
    f.write(''.join(sorted_content))
