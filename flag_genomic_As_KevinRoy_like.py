#!/usr/bin/env python
'''

Scans the genome for genomic A stretches and prints them in bed format.

Usage:

    flag_genomic_As_KevinRoy_like.py -i genome.fa -l [word_len] -A [min As within word] -maxCT [max C+Ts within word] -el [extend output left] -er [extend output right] -o [output_file]

-i ... genome in fasta format
-l ... length of motif, ie word length (default = 6)
-A ... minimum number of A for word to be flagged (default = 4)
-maxCT ... maximum number of C and Ts for word to be flagged (default = 2)
-el ... extend each single nucleotide hit by this number upstream (default = 0)
-er ... extend each single nucleotide hit by this number downstream (default = 0)
-o ... output file
-q ... quiet, ie suppress print of progression to stdout

Intended for combination as input for bedtools subtract of A-tail based 3p end sequencing.

'''


__author__ = 'schmidm'

import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fname', action='store', help="input genome file")
    parser.add_argument("-l", "--word_len", action='store', help="word length", type=int, default=5)
    parser.add_argument("-A", "--A", action='store', help="minimum As within word", type=int, default=0)
    parser.add_argument("-maxCT", "--maxCT", action='store', help="maximum C plus T within word", type=int, default=0)
    parser.add_argument("-el", "--extend_left", action='store', default=0, help="extend output hits by this upstream", type=int)
    parser.add_argument("-er", "--extend_right", action='store', default=0, help="extend output hits by this downstream", type=int)
    parser.add_argument("-o", "--output_file", action='store', default="flagged.bed", help="output bed file")
    parser.add_argument("-q", "--quiet", action='store_true', default=True, help="quiet, suppress print of progress info")

    args=parser.parse_args()
    return args

class Sequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

def read_fasta(fname):
    '''
    :param fname: fasta file
    :return: Generator of Sequence objects for the individual sequences in the file
    '''
    try:
        with open(fname, 'r') as f:
            name = ''
            seq = ''
            for line in f:
                if line[0] == '>':
                    if len(seq) > 0:
                        yield Sequence(name, seq)
                    name = line.lstrip('>').rstrip()
                    seq = ''
                else:
                    seq += line.strip()
            yield Sequence(name, seq)
    except Exception as e:
        print 'ERROR: genome fasta file ', fname, 'not found!'

    return



if __name__ == '__main__':

    args = parse_args()

    flagged_regions = []
    for seq in read_fasta(args.fname):
        if not args.quiet: print 'scanning: ', seq.name
        start = 0
        end = start + args.word_len
        cnts = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}
        cnts['A'] = sum(1 for c in seq.seq[start:end] if c == 'A' )
        cnts['C'] = sum(1 for c in seq.seq[start:end] if c == 'C' )
        cnts['G'] = sum(1 for c in seq.seq[start:end] if c == 'G' )
        cnts['T'] = sum(1 for c in seq.seq[start:end] if c == 'T' )

        while end < len(seq.seq):
            #print 'A', cnts, cnts['C'] + cnts['T']
            if cnts['A'] >= args.A and (cnts['C'] + cnts['T']) <= args.maxCT:
                flagged_regions.append( (seq.name, start-1-args.extend_left, start+args.extend_right, seq.seq[start:end], '.', '+' ))
                #print 'A', cnts, cnts['C'] + cnts['T']

            if cnts['T'] >= args.A and (cnts['G'] + cnts['A']) <= args.maxCT:
                flagged_regions.append( (seq.name, end-args.extend_right, end+1+args.extend_left, seq.seq[start:end], '.', '-' ))
                #print 'T', cnts

            cnts[seq.seq[start]] -= 1
            cnts[seq.seq[end]] += 1

            start += 1
            end += 1

    if not args.quiet: print 'scanning done, found ', len(flagged_regions), ' hits'


    with open(args.output_file, 'w') as f:
        for r in flagged_regions:
            f.write('\t'.join(str(x) for x in r))
            f.write('\n')