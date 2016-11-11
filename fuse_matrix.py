#!/usr/bin/env python
'''

Reads in 2 deeptools computeMatrix files and fuses them (binds columns as opposed to join_matrix which binds rows).
Simply avoids having to call computeMatrix again for very large matrices


Only makes sense if the matrix files were produced using same bed files ...


Usage: fuse_matrix.py [-invert2] matrix1.gz matrix2.gz matrix_fused.gz

Example: python /Users/schmidm/Documents/MS_Metagene_Tools/fuse_matrix.py matrix_plus.gz matrix_minus.gz matrix_fused.gz

writes to matrix_fused.gz if not provided it concatenated the names and uses suffix _fused.gz

'''


__author__ = 'schmidm'

import sys
import gzip
import json
import os


##number to add for correct indexing of columns given that the first 6 columns
BED_INFO_OFFSET = 6


if len(sys.argv) <= 2 or sys.argv[1] == '-h':
    print '!! NEED MORE ARGUMENTS !!'
    print __doc__
    exit()

invert2 = False
if sys.argv[2] == 'invert 2':
    invert2 = True

matrix1 = sys.argv[1]
matrix2 = sys.argv[2]

if len(sys.argv) >= 4:
    outfilename = sys.argv[3]
else:
    outfilename = matrix1.replace('.gz', '_') + os.path.split(matrix2)[1].replace('.gz', '_fused.gz')

def bounds_to_tuple(l):
    '''converts a list [0,10,20] to list of tuple [(0,10),(10,20)]'''
    return [(int(l[i])+1, int(l[i+1])) for i in range(len(l)-1)]

print 'writing joined to ', outfilename

#open sense and antisense file and iterate through an merge the contents
with gzip.open(matrix1, 'r') as f, gzip.open(matrix2, 'r') as f2:

    #parse the headers
    meta = json.loads( f.readline().strip('@') )
    meta2 = json.loads( f2.readline().strip('@') )

    #sample boundaries sense and antisense
    sample_bnds = bounds_to_tuple(meta['sample_boundaries'])
    sample_bnds2 = bounds_to_tuple(meta2['sample_boundaries'])

    grp_bnds = bounds_to_tuple(meta['group_boundaries'])

    ##loop through the 2 files alternating, line by line

    joined_content = ''
    for line in f:
        line = line.rstrip().split('\t')
        region = [line[0], line[1], line[2], line[5]]
        #print region

        line2 = f2.readline()
        line2 = line2.rstrip().split('\t')
        region2 = [line2[0], line2[1], line2[2], line2[5]]
        #print region2

        if region != region2:
            raise Exception, ('bed info for matrix1 does not fit bed info for matrix 1 region: ' + ' '.join(str(x) for x in region) + ' vs matrix2 region: ' +  ' '.join(str(x) for x in region2))

        f2_content = line2[BED_INFO_OFFSET:]

        line.extend(f2_content)
        joined_content += '\t'.join(line)
        joined_content += '\n'


#adjust sample_boundaries (ie rows for each 'group')
#print meta['sample_boundaries']

sample_bnds1 = meta['sample_boundaries']
max_sample_bnd = sample_bnds1[-1]

sample_bnds2_shifted = [ max_sample_bnd + b for b in meta2['sample_boundaries'][1:]]

joined_sample_bnds = sample_bnds1 + sample_bnds2_shifted

joined_meta = meta
joined_meta['sample_boundaries'] = joined_sample_bnds
joined_meta['sample_labels'] = meta['sample_labels'] + meta2['sample_labels']

print 'sample bounds after join: ', joined_meta['sample_boundaries']
print 'samples after join: ', joined_meta['sample_labels']
#write output
with gzip.open(outfilename, 'wb') as f:
    h=json.dumps(joined_meta,separators=(',',':'))
    hx='@'+h+'\n'
    f.write(hx)
    f.write(joined_content)
