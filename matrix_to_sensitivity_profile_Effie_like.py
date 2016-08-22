#!/usr/bin/env python
'''

Reads in a deeptools computeMatrix file and computes average Profile as done for Meola et al:
1. add pseudocount to all values
2. make average tpm for replicates
3. make ratio kd over control
4. set ratio 0 for all cells in matrix with pseudocount in both kd and control
5. average over all genes

Usage: matrix_to_profile_Effie_like.py matrix.gz control_label [sample_label1,sample_label2,...] pseudocount outname

writes to matrix"_Effie_profiles.txt"

'''


__author__ = 'schmidm'

import sys
import gzip
import json
import numpy as np


if len(sys.argv) <= 3 or sys.argv[1] == '-h':
    print __doc__
    exit()

matrix_name = sys.argv[1]
ctrl_name = sys.argv[2]
sample_names = sys.argv[3].split(',')
pseudocount = float(sys.argv[4])

if len(sys.argv) >= 6:
    outfilename = sys.argv[5]
else:
    outfilename = matrix_name.replace('.gz', '_sensitivity.gz')

def bounds_to_tuple(l):
    '''converts a list [0,10,20] to list of tuple [(0,10),(10,20)]'''
    return [(int(l[i])+1, int(l[i+1])) for i in range(len(l)-1)]

#sense file ... nothing to be done except keep the contents
with gzip.open(matrix_name, 'r') as f:

    #parse the headers
    meta = json.loads( f.readline().strip('@') )

    print meta

    #group boundaries sense and antisense
    grp_bnds = bounds_to_tuple(meta['group_boundaries'])

    #sample boundaries sense and antisense
    sample_bnds = bounds_to_tuple(meta['sample_boundaries'])

    controls = [ (i, sample) for i, sample in enumerate(meta['sample_labels']) if ctrl_name in sample]
    ctrl_idx_starts = [sample_bnds[i[0]][0]+5 for i in controls]

    grp_kds = [ [(i, sample) for i, sample in enumerate(meta['sample_labels']) if s in sample] for s in sample_names ]
    kd_idx_starts = [ [sample_bnds[i[0]][0]+5 for i in kd] for kd in grp_kds ]

    print sample_names
    print 'controls found: ', controls
    print 'grp_kds found: ', grp_kds
    print 'using start idx for control: ', ctrl_idx_starts #[sample_bnds[i[0]] for i in controls]
    print 'using start idx for kds: ', kd_idx_starts #[[sample_bnds[i[0]] for i in kd] for kd in grp_kds]
    print 'sample_bounds: ', sample_bnds

    bin = int(meta['bin size'])
    up = -int(meta['upstream'])
    dn = int(meta['downstream'])
    print 'bin, up and down: ', bin, up, dn


    #adjust group_boundaries (ie lines for each 'group')
    sens_meta = meta

    sens_meta['sample_labels'] = sample_names
    sens_meta['sample_boundaries'] = meta['sample_boundaries'][0:len(sample_names)+1]

    print 'sens_meta_sample_bounds', sens_meta['sample_boundaries']


    with gzip.open(outfilename, 'wb') as outfile:
        h=json.dumps(sens_meta,separators=(',',':'))
        hx='@'+h+'\n'
        outfile.write(hx)

        lnr = 1
        for line in f:
            line = line.strip('\n').split()

            #ctrl_means = [0 for i in range(up, dn, bin)]
            #print len(ctrl_means)
            sensitivities = [[0 for i in range(up, dn-1, bin)] for k in grp_kds]
            for i, rel_pos in enumerate(range(up, dn-1, bin)):

                idx_ctrls = [idx + i for idx in ctrl_idx_starts]
                #print 'idx_ctrls', idx_ctrls
                ctrl = np.mean([ float(line[idx]) + pseudocount if line[idx] != 'nan' else pseudocount for idx in idx_ctrls ])


                #for each kd:
                sensitivity = [0 for k in grp_kds]
                for j, kd_idx in enumerate(kd_idx_starts):
                    idx_kd = [idx + i for idx in kd_idx]
                    kd_mean = np.mean([ float(line[idx]) + pseudocount if line[idx] != 'nan' else pseudocount for idx in idx_kd ])

                    sensitivity[j] = kd_mean / ctrl
                    sensitivities[j][i] = kd_mean / ctrl

            new_line = '\t'.join(line[0:6])
            new_line += '\t'
            new_line += '\t'.join(str(v) for sens in sensitivities for v in sens)
            new_line += '\n'
            outfile.write(new_line)

            lnr += 1