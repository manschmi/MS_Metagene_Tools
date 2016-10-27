#!/usr/bin/env python
'''

usage:
        python map_motif.py <motif> <fasta_file> [-bed:bed_file] [-count]

finds motif in fasta_file and reports the positions found for each sequence in fasta file

output is of type
bed_info\tab\comma-delimited positions
when working without bed file
chr\tab\1\tab\chr_length\tab\comma-delimited positions


-count not implemented yet:: >>if option  only sum of all hits to a single 'count'

if a bed file is specified only searches in bed file intervals instead of entire sequence in fasta file.

automatically takes strand information into account if its present in the bed file, otherwise only scans plus strand of each interval.

!! Input bed file is assumed to be and output positions are 0-based indices !!
!! on the minus strand the position refers to the hit of the first based in the motif (ie right end on minus strand and left end on plus strand)
!!sequence
>chr1
GATCATG

bed
chr1 0 6 . . +
chr1 0 6 . . -

motif AT is reported at
chr1 0 6 . . + 1,4  (G-A-TC-A-TG)
chr1 0 6 . . - 2,5  (GA-T-CA-T-G)

'''

__author__ = 'schmidm'

import sys
import argparse



class Sequence:
    '''single sequence with name'''
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.com_base = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    def rev_com(self, s):
        '''reverse complement string of sequence s'''
        s = s[::-1]
        return ''.join(self.com_base[c] for c in s if c!='\n' )

    def get_interval(self, strand, start, stop):
        '''returns sequence interval from chromosome, end included'''
        if start < 0 or stop > len(self.seq)-1:
            return
        s = self.seq[start:stop]
        if strand == '-':
            return self.rev_com(s)
        else:
            return s


class Interval:
    '''a genomic interval with chr start end col4 col5 strand and extra info '''
    def __init__(self, line):
        line = line.split()
        self.name = line[0]
        self.start = int(line[1])
        self.end = int(line[2])
        if len(line) >= 4:
            self.col4 = line[3]
        else: self.col4 = '.'
        if len(line) >= 5:
            self.col5 = line[4]
        else: self.col5 = '.'
        if len(line) >= 6:
            self.strand = line[5]
        else:
            self.strand = '+'
        if len(line) > 6:
            self.extra = line[6:]
        else: self.extra = []

    def __str__(self):
        s=self.name+'\t'+str(self.start)+'\t'+str(self.end)
        if self.col4:
            s += '\t'+self.col4
        if self.col5:
            s += '\t' +self.col5
        if self.strand:
            s+= '\t'+self.strand
        if self.extra:
            s += '\t'+'\t'.join(self.extra)
        return s


def read_fasta(file):
    '''iterator for Sequence objects in a fasta file'''
    seq=[]
    name = None
    for line in file:
        if line[0]=='>':
            if name: #skip first line
                yield Sequence(name, seq)
            name = line.lstrip('>').rstrip()
            seq=''
        else:
            seq += line.strip()

    yield Sequence(name, seq) #last seq
    return


def read_bed(file):
    '''iterator for list of elements of each line in a bed file'''

    for line in file:
        yield Interval(line)

    return



def find_all(a_str, sub):
    '''index of all occurrences of sub in a_str'''
    start = 0
    result = []
    while True:
        start = a_str.find(sub, start)
        if start == -1: return result
        result.append(start)
        start += len(sub)
    return result




########################################################################
########INPUT###########################################################
########################################################################

#python map_motif.py <motif> <fasta_file> [-bed:bed_file]  [-s] [-count]

parser = argparse.ArgumentParser()

parser.add_argument('motif')
parser.add_argument('fasta_file')

parser.add_argument('-bed', default=None)
parser.add_argument('-s', action="store_true", default=True)
parser.add_argument('-count', action="store_true", default=False)

args = parser.parse_args()

if args.bed:
    with open(args.bed, 'r') as b:
        bed_intervals = [i for i in read_bed(b)]


with open(args.fasta_file, 'r') as f:
    for seq in read_fasta(f):
        if args.bed:
            intvls = [i for i in bed_intervals if i.name == seq.name]
        else:
            intvls = Interval(seq.name+' 0 '+str(len(seq.seq)))

        if len(intvls) == 0: continue
        for intvl in intvls:
            s = seq.get_interval(intvl.strand, intvl.start, intvl.end)
            hits = find_all(s, args.motif)
            if intvl.strand == '-':
                offset = intvl.end - intvl.start + 1
                hits = [offset - hit for hit in reversed(hits)]
            print str(intvl) + '\t' + ','.join(str(i) for i in hits)
