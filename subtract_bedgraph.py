#!/usr/bin/env python
'''

subtracts values from file2 from file1

Usage: subtract_bedgraph.py file1 file2 --chunk_size --value_col -o/--outFileName

example: python subtract_bedgraph.py ChIP.bedgraph input.bedgraph --chunk_size 1000000 --value_col 5 -o ChIP_minus_input.bedgraph


'''


__author__ = 'schmidm'


import sys
import argparse
import numpy as np


parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument('file1')
parser.add_argument('file2')
parser.add_argument('--chunk_size', default=1000000, type=int, help="size of the internal array, allows for minor optimization, probably most efficient to use highest end position in bed file (optional, default = 1000000)")
parser.add_argument('--value_col', default=4, type=int, help="UPS: IMPORTANT, which column to get values for subtraction (default = 4)")
parser.add_argument('-o', '--outFileName', type=str)

args = parser.parse_args()


def print_subtracted(chr, f1_values, f2_values):
    '''
    :param chr: name of chr for output
    :param f1_values: np.array raw values
    :param f2_values: np.array values to be subtracted
    :return: nothing

    Prints the values in collapsed form to args.outFileName
    '''
    if len(f1_values) > len(f2_values):
        f2_values = np.append(f2_values, np.zeros(len(f1_values)-len(f2_values)))
    elif len(f2_values) > len(f1_values):
        f1_values = np.append(f1_values, np.zeros(len(f2_values)-len(f1_values)))

    out_values = f1_values - f2_values

    with open(args.outFileName, 'a') as outfile:
        start = 0
        value = out_values[0]
        i = 0
        while i < len(out_values):
            while i < len(out_values) and out_values[i] == value:
                i += 1
            if value != 0:
                outfile.write(chr + '\t' + str(start) + '\t' + str(i) + '\t' + str(value) + '\n')
            if i == len(out_values):
                return
            start = i
            value = out_values[i]


def read_chr(file):
    '''
    :param file: bedgraph file open for read
    :return: iterator of tuple chr, array of values
    '''
    ar = np.zeros(args.chunk_size)
    line = file.next().split('\t')
    chr = line[0]
    while True:
        try:
            while line[0] == chr:
                if int(line[2]) >= len(ar):
                    ar = np.append(ar, np.zeros(len(ar)))
                ar[int(line[1]): int(line[2])] += float(line[3])
                line = file.next().split('\t')
            yield chr, ar
            ar = np.zeros(args.chunk_size)
            chr = line[0]
        except StopIteration:
            yield chr, ar
            return


with open(args.outFileName, 'w') as f:
    f.seek(0)
    f.truncate()


with open(args.file1, 'r') as file1, open(args.file2, 'r') as file2:
    f2_iter = read_chr(file2)
    for chr1, f1_values in read_chr(file1):
        try:
            chr2, f2_values = f2_iter.next()
            while chr2 < chr1:
                chr2, f2_values = f2_iter.next()
            if chr2 == chr1:
                print_subtracted(chr1, f1_values, f2_values)
        except:
            f2_values = np.zeros(args.chunk_size)
            print_subtracted(chr1, f1_values, f2_values)
            exit()