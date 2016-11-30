#!/usr/bin/env python
'''

Reads in a deeptools computeMatrix file and sorts it according to the order in a provided bed file.


Only makes sense if the matrix and bed files contain the same regions, or if bed is a subset of regions in matrix.
Any regions present in bed but not in matrix will be dropped with a warning, specifying the region.
Regions present in matrix but not bed will be counted and dropped with a warning but not specified.

Groupings are lost in this operation!!

Usage: sort_matrix_using_bed.py matrix.gz bed.bed matrix_bedsorted.gz

Example: python /Users/schmidm/Documents/MS_Metagene_Tools/sort_matrix_using_bed.py matrix.gz order.bed matrix_bedsorted.gz

writes to matrix_bedsorted.gz if not provided it concatenated the names and uses suffix _bedsorted.gz

'''


__author__ = 'schmidm'

import sys
import gzip
import json


##number to add for correct indexing of columns given that the first 6 columns are the region
BED_INFO_OFFSET = 5


if len(sys.argv) <= 2 or sys.argv[1] == '-h':
    print '!! NEED MORE ARGUMENTS !!'
    print __doc__
    exit()

matrix = sys.argv[1]
bed = sys.argv[2]

if len(sys.argv) >= 4:
    outfilename = sys.argv[3]
else:
    outfilename = matrix.replace('.gz', '_bedsorted.gz')

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

print 'writing resorted to ', outfilename

#open sense and antisense file and iterate through an merge the contents
with gzip.open(matrix, 'r') as f:

    #parse the headers
    meta = json.loads( f.readline().strip('@') )

    #sample boundaries sense and antisense
    sample_bnds = bounds_to_tuple(meta['sample_boundaries'])

    ##loop through the bed files and find matching row of matrix for each

    matrix = [line.split('\t') for line in f]



sorted_content = []
bed_region_cnt=0
matrix_region_cnt=0
#bed_regions_found=0
for region in bed_lines(bed):
    bed_region_cnt += 1
    # print 'REGION'
    #if bed_region_cnt % 10: print '.'
    #if bed_region_cnt % 100: print bed_region_cnt / 100
    found = False
    for line in matrix:
        if line[0] == region[0] and line[1] == region[1] and line[2] == region[2] and line[5] == region[5]:
            #print 'FOUND!'
            found = True
            sorted_content.append('\t'.join(line))
            matrix_region_cnt += 1
            break

    if not found:
        print 'WARNING region from bed file not found in matrix: ', region
        matrix_region_cnt += 1
        bed_region_cnt -= 1



if matrix_region_cnt != bed_region_cnt:
    print 'WARNING matrix contained ', matrix_region_cnt-bed_region_cnt, ' regions not present in bed file'

print 'found', bed_region_cnt, ' regions in the matrix out of ', meta['group_boundaries'][1], 'regions in the bed file'

meta['group_boundaries'] = [0,bed_region_cnt]

#write output
with gzip.open(outfilename, 'wb') as f:
    h=json.dumps(meta,separators=(',',':'))
    hx='@'+h+'\n'
    f.write(hx)
    f.write(''.join(sorted_content))
