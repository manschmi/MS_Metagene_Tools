#!/usr/bin/env python
'''

Reads in 2 deeptools computeMatrix files and joins them.
Only makes sense if the matrix files were produced using bed and bw files from plus and minus strand...

Note deeptools does not take strandedness of a bw file into account...
Useful for ie stranded RNAseq data.
Note deeptools treats the strand information from the bed file correct (ie TSS of a minus strand gene will be the end of the annotation ...)

If option -invert2 -> inverts the rows from the 2nd file.
This is typically not necessary if the bed file used to compute the matrix has strand information,
ie deeptools computeMatrix then inverts the ranges by default.


Usage: join_matrix.py [-invert2] matrix1.gz matrix2.gz matrix_joined.gz

Example: python /Users/schmidm/Documents/MS_Metagene_Tools/join_matrix.py matrix_plus.gz matrix_minus.gz matrix_joined.gz

writes to matrix_joined.gz if not provided it concatenated the names and uses suffix _joined.gz

'''


__author__ = 'schmidm'

import sys
import gzip
import json


##number to add for correct indexing of columns given that the first 6 columns
BED_INFO_OFFSET = 5


if len(sys.argv) <= 2 or sys.argv[1] == '-h':
    print '!! NEED MORE ARGUMENTS !!'
    print __doc__
    exit()

invert2 = False
if sys.argv[2] == 'invert 2':
    invert2 = True

sense_name = sys.argv[1]
antisense_name = sys.argv[2]

if len(sys.argv) >= 4:
    outfilename = sys.argv[3]
else:
    outfilename = sense_name.replace('_plus.gz', '_') + antisense_name.replace('_minus.gz', '_joined.gz')

def bounds_to_tuple(l):
    '''converts a list [0,10,20] to list of tuple [(0,10),(10,20)]'''
    return [(int(l[i])+1, int(l[i+1])) for i in range(len(l)-1)]

print 'writing joined to ', outfilename

#open sense and antisense file and iterate through an merge the contents
with gzip.open(sense_name, 'r') as f, gzip.open(antisense_name, 'r') as f_as:

    #parse the headers
    meta = json.loads( f.readline().strip('@') )
    as_meta = json.loads( f_as.readline().strip('@') )

    if meta['group_labels'] != as_meta['group_labels']:
        print 'cannot join files with different group labels'

    #group boundaries sense and antisense
    grp_bnds = bounds_to_tuple(meta['group_boundaries'])
    as_grp_bnds = bounds_to_tuple(as_meta['group_boundaries'])


    #sample boundaries sense and antisense
    sample_bnds = bounds_to_tuple(meta['sample_boundaries'])
    as_sample_bnds = bounds_to_tuple(as_meta['sample_boundaries'])

    if sample_bnds != as_sample_bnds:
        print 'cannot join files with different samples...'


    ##loop through the 2 files alternating, group by group
    lnr = 1 #line number
    as_lnr = 1

    joined_content = ''
    for grp in range(len(grp_bnds)):
        max_sense_line = grp_bnds[grp][1]
        if lnr != grp_bnds[grp][0]:
            print 'something wrong with indexing sense'
            exit()
        while lnr <= max_sense_line:
            joined_content += f.readline()
            lnr += 1

        max_as_line = as_grp_bnds[grp][1]
        if as_lnr != as_grp_bnds[grp][0]:
            print 'something wrong with indexing antisense'
            exit()
        while as_lnr <= max_as_line:
            if not invert2:
                joined_content += f_as.readline()
            else:
                line = f_as.readline()
                line = line.rstrip().split('\t')
                for i in range( len(sample_bnds) ):
                    left = BED_INFO_OFFSET + int(sample_bnds[i][0])
                    right = BED_INFO_OFFSET + int(sample_bnds[i][1])
                    line[left:right] = reversed(line[left:right])
                joined_content += '\t'.join(line)
                joined_content += '\n'
            as_lnr += 1


#adjust group_boundaries (ie lines for each 'group')
bs1 = meta['group_boundaries']
bs2 = as_meta['group_boundaries']

joined_group_bnds = [(b + bs2[i]) for i, b in enumerate(bs1)]

joined_meta = meta
joined_meta['group_boundaries'] = joined_group_bnds
joined_meta['group_labels'] = [str.replace('_plus', '') for str in joined_meta['group_labels']]
joined_meta['sample_labels'] = [str.replace('_plus', '') for str in joined_meta['sample_labels']]


#write output
with gzip.open(outfilename, 'wb') as f:
    h=json.dumps(joined_meta,separators=(',',':'))
    hx='@'+h+'\n'
    f.write(hx)
    f.write(joined_content)
