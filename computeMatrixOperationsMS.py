#!/usr/bin/env python
import deeptools.heatmapper as heatmapper
import re
import numpy as np
import argparse
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
This tool performs a variety of custom operations not supported by the official deeptools computeMatrixOperations
on files produced by computeMatrix.

detailed help:
    AVAILABLE TOOLS:
        - scaleBy
        - negateMinusStrandValues
        - nanToValue
        - addPseudoCount
        - binarize
        - log2
        - averageSamples
        - ratioToSample
        - sortUsingBed
        - subset
        - trimLabels


DETAILS:

  computeMatrixOperationsMS scaleBy <value>
        multiplies all values by a specific value, requires a value, ie -1 for negation

  computeMatrixOperationsMS negateMinusStrandValues
        multiplies all values from minus strand regions with -1 for negation

  computeMatrixOperationsMS nanToValue <value>
        sets all nan bins to a specific value, requires a value, typically 0

  computeMatrixOperationsMS addPseudoCount <value>
        if value is skipped it will use the minimum positive value in the matrix

  computeMatrixOperationsMS binarize
        all positive values set to 1 and all negative and nans to 0

  computeMatrixOperationsMS log2
        all positive values transformed log2 all others to nan

  computeMatrixOperationsMS averageSamples=regex1,regex2,...
        searches for samples matching each regex and creates average for each bin in matrix
        samples not matched by any of the regexes are removed
        output sample names are the regex used

  computeMatrixOperationsMS diffToSample=regex
        display matrix as difference relative to the sample
        UPS: only regex matching a single sample are allowed!

  computeMatrixOperationsMS ratioToSample=regex
        display matrix as ratio relative to the sample
        UPS: only regex matching a single sample are allowed!

  computeMatrixOperationsMS sortUsingBed=file1.bed,file2.bed;groupByColumn:6
        resort using the specified bed files.
        optional groupByColumn, if the bed files have column with group identifiers these will
        be used for grouping (ups: group column index is 0 based)

  computeMatrixOperationsMS subset=samples:sample_regex1,sample_regex2;groups:group_regex1,group_regex2
        subsets the matrix for samples and groups, order or argument is obeyed, ie can be used for resorting
        note: regex will NOT be used as sample and group names in the output
        (these are preserved as multiple hits are possible)
        UPS: not that commas are now allowed in the regexes for now

  computeMatrixOperationsMS trimLabels=samples:sample_regex1,sample_regex2;groups:group_regex1,group_regex2
        trims the specified regexes from sample and group labels

OPERATION CHAINING:
-> it is possible to chain multiple operations (ie to avoid saving intermediate matrices!!)
  computeMatrixOperations addPseudoCount averageSamples=eGFP,KD1,KD2 ratioToSample=eGFP log2

""",
        epilog='example usages:\n'
               'python computeMatrixOperationsMS -m mat.gz -o mat_mod.gz addPseudoCount averageSamples=eGFP,KD1,KD2 ratioToSample=eGFP log2 \n\n'
               ' \n\n')

    parser.add_argument('--matrixFile', '-m',
                          help='Matrix files from the computeMatrix tool.',
                          required=True)

    parser.add_argument('--outFileName', '-o',
                          help='Output file name',
                          required=True)

    parser.add_argument('operations', type=str, nargs='+',
                        help='the operations to perform')

    return parser


def perform_operations(args, matrix):
    operations_dict = {'scaleBy': scaleValue,
                       'negateMinusStrandValues': negateMinusStrandValues,
                       'nanToValue': nanToValue,
                       'addPseudoCount': add_pseudocount,
                       'log2': log2,
                       'binarize': binarize,
                       'averageSamples': averageSamples,
                       'diffToSample': diffToSample,
                       'ratioToSample': ratioToSample,
                       'sortUsingBed': sortUsingBed,
                       'subset': subset,
                       'trimLabels': trimLabels}
    for arg in args.operations:
        arg_list = arg.split('=')
        op = arg_list[0]

        if not op in operations_dict:
            raise Exception('found unsupported operation: '+op)

        print('performing operation: ' + op)
        if len(arg_list) > 1:
            op_arg_str = arg_list[1]
        else:
            op_arg_str = None
        operations_dict[op](matrix, op_arg_str)


def scaleValue(matrix, scaleFactor):
    """
    converts all nan bins to a specific value
    """
    if not scaleFactor:
        raise Exception('cannot perform scaling without a value')

    print('  scaling by value: ' + str(scaleFactor))

    matrix.matrix *= float(scaleFactor)

    return


def negateMinusStrandValues(matrix, pad=None):
    """
    multiplies all values for minus strand regions *-1
    """

    neg_strand_indeces = [i for i,r in enumerate(matrix.regions) if r[4] == '-']
    print('  negating values from ' + str(len(neg_strand_indeces)) + ' minus strand regions')
    matrix.matrix[neg_strand_indeces,] = matrix.matrix[neg_strand_indeces,] * (-1)

    return


def nanToValue(matrix, nanValue):
    """
    converts all nan bins to a specific value
    """
    if not nanValue:
        raise Exception('cannot perform nanToValue without a value')

    matrix.matrix[np.isnan(matrix.matrix)] = float(nanValue)

    return


def add_pseudocount(matrix, pseudoCountStr=None):
    """
    Add pseudoCount to Matrix, if pseudCount is None finds minimum positive value in matrix and adds that.
    Note: nans are kept as nans!!
    """

    if not pseudoCountStr:
        pseudoCount = min_pos_val(matrix.matrix)
    else:
        pseudoCount = float(pseudoCountStr)

    print('  adding pseudocount: ' + str(pseudoCount))
    matrix.matrix[matrix.matrix > 0] += pseudoCount
    matrix.matrix[matrix.matrix <= 0] = pseudoCount

    return


def min_pos_val(mat):
    """
    returns minimum positive value in a numpy matrix
    """
    return np.nanmin(mat[mat > 0])


def log2(matrix, pad=None):
    """
    Converts all bins to log2 values, 0 and negatives are nan
    Note: nans are kept as nans!!
    """
    matrix.matrix = np.log2(matrix.matrix)

    return


def binarize(matrix, pad=None):
    """
    binarize matrix, ie all positive bins are set to 1 and all 0, negative and nan bins are set to 0
    """
    matrix.matrix[matrix.matrix > 0] = 1
    matrix.matrix[matrix.matrix <= 0] = 0
    matrix.matrix[np.isnan(matrix.matrix)] = 0

    return


def averageSamples(matrix, regex_str):
    """
    Given the sample labels (substrings), average over all samples fitting each of them
    :param regex_str comma separated regexes ie 'KD1.*in,KD2.*in,KD3.*in'
    """
    bounds = matrix.sample_boundaries
    sample_labels = matrix.sample_labels
    ncol_per_sample = bounds[1]
    nrow = matrix.matrix.shape[0]
    mat = None

    labels_list = regex_str.split(',')
    for label in labels_list:
        samples_match = [(i, sample) for i, sample in enumerate(sample_labels) if re.search(label, sample)]
        print('  averaging samples:\n    ' + '\n    '.join([s[1] for s in samples_match]))
        start_cols = np.array([bounds[s[0]] for s in samples_match])
        mats = [np.nanmean(matrix.matrix[:,start_cols+col], axis=1, keepdims=True) for col in range(ncol_per_sample)]
        mats = np.hstack(mats)
        mats.shape = (nrow, ncol_per_sample)
        if mat is not None:
            mat = np.hstack([mat, mats])
        else:
            mat = mats
            mat.shape = (nrow, ncol_per_sample)

    matrix.matrix=mat
    matrix.sample_boundaries = [ncol_per_sample*i for i in range(len(labels_list)+1)]
    matrix.sample_labels = labels_list

    return


def diffToSample(matrix, label):
    """
    Given a sample label (substring interpreted as regex), substracts values from sample from all others and drops it
    """

    bounds = matrix.sample_boundaries
    sample_labels = matrix.sample_labels
    ref_idx = [i for i, sample in enumerate(sample_labels) if re.search(label, sample)]

    if len(ref_idx) == 0:
        raise Exception('found no sample matching diffToSample argument: ', label)
    elif len(ref_idx) > 1:
        raise Exception('found more than one sample to subtract: ', ','.join(str(i) for i in ref_idx))
    else:
        ref_idx = ref_idx[0]

    ref_cols = np.array(range(bounds[ref_idx],bounds[ref_idx+1]))
    ref_mat = matrix.matrix[:, ref_cols]
    matrix.matrix = np.delete(matrix.matrix, ref_cols, axis=1)

    col_list = [np.array(range(bounds[i], bounds[i+1])) for i in range(len(sample_labels)-1)]

    mats = [matrix.matrix[:, cols] - ref_mat for cols in col_list]
    matrix.matrix = np.concatenate(mats, axis=1)

    matrix.sample_boundaries.pop()
    matrix.sample_labels.pop(ref_idx)
    matrix.sample_labels = [l + '-' + label for l in matrix.sample_labels]

    return


def ratioToSample(matrix, label):
    """
    Given a sample label (substring interpreted as regex), normalize all others to this sample and drop it
    """

    bounds = matrix.sample_boundaries
    sample_labels = matrix.sample_labels
    ref_idx = [i for i, sample in enumerate(sample_labels) if re.search(label, sample)]

    if len(ref_idx) == 0:
        raise Exception('found no sample matching ratioToSample argument: ', label)
    elif len(ref_idx) > 1:
        raise Exception('found more than one sample to normalize to: ', ','.join(str(i) for i in ref_idx))
    else:
        ref_idx = ref_idx[0]

    ref_cols = np.array(range(bounds[ref_idx],bounds[ref_idx+1]))
    ref_mat = matrix.matrix[:, ref_cols]
    matrix.matrix = np.delete(matrix.matrix, ref_cols, axis=1)

    col_list = [np.array(range(bounds[i], bounds[i+1])) for i in range(len(sample_labels)-1)]

    mats = [matrix.matrix[:, cols]/ref_mat for cols in col_list]
    matrix.matrix = np.concatenate(mats, axis=1)

    matrix.sample_boundaries.pop()
    matrix.sample_labels.pop(ref_idx)
    matrix.sample_labels = [l + '/' + label for l in matrix.sample_labels]

    return


def sortUsingBed2(matrix, arg_str):
    """
    DEPRECATED VERSION!!
    Given a bedfile name resort the matrix using the bed file intervals, assumes unique intervals or matching names
    """

    args = arg_str.split(';')
    bed_fname = args[0]
    print('  resort matrix using BED: ' + bed_fname)

    regions = loadBED(bed_fname)
    mat_keys = [deeptools_region_str(r) for r in matrix.regions]
    try:
        rstrs = [tabbed_BED_region_str(r) for r in regions]
    except:
        raise Exception('could not interpret BED regions used for sorting')

    try:
        mat_order = [mat_keys.index(rstr) for rstr in rstrs if rstr in mat_keys ]
        if len(rstrs) > len(mat_order):
            failed = [rstr for rstr in rstrs if rstr not in mat_keys]
            print('NOTE: regions from bed file not found in matrix: ' + '\n'.join(failed))
        elif len(mat_order) < len(mat_keys):
            failed = [mat_key for mat_key in mat_keys if mat_key not in rstrs]
            print('NOTE: regions from matrix not found in bed files: ' + '\n'.join(failed))
        matrix.matrix = matrix.matrix[np.ma.array(mat_order),]
        matrix.regions = [matrix.regions[i] for i in mat_order]
    except:
        raise Exception('reordering of regions failed')

    #appending group information
    if len(args) > 1:
        arg1 = args[1].split(':')
        if arg1[0] == 'groupByColumn':
            groupColumn = int(arg1[1])
            print('  grouping using column (ups: 0-based!!): ' + str(groupColumn))
            #try:
            region_to_group = {tabbed_BED_region_str(r): r[groupColumn] for r in regions}
            #groups=set([r[4] for r  in regions])
            #group_boundaries = [0]
            #group_labels=[]
            #for group in groups:
            #   grouped_regions[group] = [(i,tabbed_BED_region_str(r)) for i, r in enumerate(regions) if r[4] == group]
            #   group_boundaries.append(len(grouped_regions[group]+group_boundaries[-1])
            #   group_labels.append(group)
            #mat_order = [x[0] for x in grouped_regions[group] for group in groups]
            #matrix.matrix = matrix.matrix[mat_order,]
            #matrix.regions = [matrix.regions[i] for i in mat_order]
            #
            i = 0
            matrix.group_boundaries = []
            grp_name = ''
            matrix.group_labels = []
            for region in matrix.regions:
                rstr = deeptools_region_str(region)
                try:
                    region_grp = region_to_group[rstr]
                    if region_grp != grp_name:
                        matrix.group_labels.append(region_grp)
                        matrix.group_boundaries.append(i)
                        grp_name = region_grp
                    i+=1
                except:
                    raise Exception('did not find region group for region: ' + rstr + ' !! This should never happen!')
            matrix.group_boundaries.append(i)
            #except:
            #    raise Exception('regrouping failed')
    return


def sortUsingBed(matrix, arg_str):
    """
    Given a bedfile name resort the matrix using the bed file intervals, assumes unique intervals or matching names
    """

    arg_list = arg_str.split(';')
    bed_fname = arg_list[0]
    print('  resort matrix using BED: ' + bed_fname)

    regions = loadBED(bed_fname)
    mat_keys = [deeptools_region_str(r) for r in matrix.regions]
    try:
        rstrs = [tabbed_BED_region_str(r) for r in regions]
    except:
        raise Exception('could not interpret BED regions used for sorting')

    if len(arg_list) == 1:
        try:
            mat_order = [mat_keys.index(rstr) for rstr in rstrs if rstr in mat_keys ]
            if len(rstrs) > len(mat_order):
                failed = [rstr for rstr in rstrs if rstr not in mat_keys]
                print('NOTE: regions from bed file not found in matrix: ' + '\n'.join(failed))
            elif len(mat_order) < len(mat_keys):
                failed = [mat_key for mat_key in mat_keys if mat_key not in rstrs]
                print('NOTE: regions from matrix not found in bed files: ' + '\n'.join(failed))
            matrix.matrix = matrix.matrix[np.ma.array(mat_order),]
            matrix.regions = [matrix.regions[i] for i in mat_order]
        except:
            raise Exception('reordering of regions failed')

    #appending group information
    else:
        arg1 = arg_list[1].split(':')
        if arg1[0] == 'groupByColumn':
            groupColumn = int(arg1[1])
            print('  grouping using column (ups: 0-based!!): ' + str(groupColumn))
            try:
                groups=set(r[4] for r  in regions)
                print('  found groups: ' + ' '.join(groups))
                region_group = {tabbed_BED_region_str(r): r[4] for r in regions}
                print(region_group)
                mat_row_regions = [region_group[deeptools_region_str(region)] for region in matrix.regions]
                print(mat_row_regions)
                group_boundaries=[0]
                group_labels=list(groups)
                grouped_row_order = []
                for group in groups:
                   rows_for_group = [i for i, grp in enumerate(mat_row_regions) if grp == group]
                   print(len(rows_for_group))
                   grouped_row_order.extend(rows_for_group)
                   group_boundaries.append(len(grouped_row_order))
                matrix.matrix = matrix.matrix[grouped_row_order,]
                matrix.regions = [matrix.regions[i] for i in grouped_row_order]
                matrix.group_boundaries = group_boundaries
                matrix.group_labels = group_labels
            except:
                raise Exception('regrouping failed')
    return


def subset(matrix, arg_str):
    args = arg_str.split(';')
    sample_regexes = group_regexes = None
    for arg in args:
        arg_name, regexes = arg.split(':')
        if arg_name == 'samples':
            sample_regexes = regexes.strip().split(',')
            print(sample_regexes[0])
            print('  trying to find sample with regexes: ' + ', '.join(sample_regexes))
        elif arg_name == 'groups':
            group_regexes = regexes.split(',')

    subset_matrix(matrix, sample_regexes, group_regexes)


def trimLabels(matrix, arg_str):
    '''
    :param matrix: the matrix
    :param arg_str: of type ""
    :return:
    '''
    args = arg_str.split(';')
    sample_regexes = group_regexes = None
    for arg in args:
        arg_name, regexes = arg.split(':')
        if arg_name == 'samples':
            sample_regexes = regexes.strip().split(',')
            for regex in sample_regexes:
                matrix.sample_labels = [re.sub(regex, '', sample) for sample in matrix.sample_labels]
        elif arg_name == 'groups':
            group_regexes = regexes.split(',')
            for regex in group_regexes:
                matrix.group_labels = [re.sub(regex, '', group) for group in matrix.group_labels]


def subset_matrix(matrix, sample_regex_list=None, group_regex_list=None):
    """
    Given the sample_regex_str and group_regex_str return matrix containing only samples and groups fitting each of them
    if None, nothing is done, ie all samples or groups are selected
    :param sample_regex_str, group_regex_str comma separated regexes ie 'KD1.*in,KD2.*in,KD3.*in'
    """

    if sample_regex_list:
        col_range = []
        sample_sel = []
        sample_bounds = matrix.sample_boundaries
        sample_labels = matrix.sample_labels
        ncol_per_sample = sample_bounds[1]
        for label in sample_regex_list:
            for i, sample in enumerate(sample_labels):
                if re.search(label, sample):
                    col_range.extend(range(sample_bounds[i], sample_bounds[i+1]))
                    sample_sel.append(sample)

        matrix.matrix = matrix.matrix[:,np.array(col_range)]
        matrix.sample_labels = sample_sel
        matrix.sample_boundaries = [ncol_per_sample * i for i in range(len(sample_sel) + 1)]

        print('  subset to samples: ' + ', '.join(sample_sel))

    if group_regex_list:
        row_range = []
        group_sel = []
        group_sel_bounds = [0]
        group_bounds = matrix.group_boundaries
        group_labels = matrix.group_labels
        for label in group_regex_list:
            for i, group in enumerate(group_labels):
                if re.search(label, group):
                    row_range.extend(range(group_bounds[i], group_bounds[i+1]))
                    group_sel.append(group)
                    group_sel_bounds.append(group_sel_bounds[-1] + (group_bounds[i+1] - group_bounds[i]) )

        row_range = np.array(row_range)
        matrix.matrix = matrix.matrix[row_range,]
        matrix.regions = [matrix.regions[i] for i in row_range]
        matrix.group_labels = group_sel
        matrix.group_boundaries = group_sel_bounds

        print('  subset to groups: ' + ', '.join(group_sel))

    return



def loadBED(fname):
    with open(fname, 'r') as f:
        regions = [line.rstrip().split('\t') for line in f if line[0] != '#']
    return regions


def tabbed_BED_region_str(region):
    '''
    tab-separated region to string
    '''
    return region[0]+':'+region[1]+'-'+region[2]+'('+region[5]+')_'+region[3]

def tabbed_BED_region_to_deeptools(region):
    '''
    tab-separated region to string
    '''
    return [region[0], region[1], region[2], region[3]]

def deeptools_region_str(region):
    '''
    tab-separated region to string
    '''
    if isinstance(region, dict):
        return region['chrom']+':'+str(region['start'])+'-'+str(region['end'])+'('+region['strand']+')_'+region['name']
    elif isinstance(region, list):
        return region[0] + ':' + str(region[1][0][0]) + '-' + str(region[1][-1][1]) + '(' + region[4] + ')_' + \
               region[2]



if __name__ == '__main__':

    args = parse_arguments().parse_args(sys.argv[1:])
    hm = heatmapper.heatmapper()
    hm.read_matrix_file(args.matrixFile)

    perform_operations(args, hm.matrix)

    hm.parameters['group_labels'] = hm.matrix.group_labels
    hm.parameters['group_boundaries'] = hm.matrix.group_boundaries

    hm.save_matrix(args.outFileName)
