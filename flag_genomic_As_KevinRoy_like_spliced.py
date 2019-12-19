#!/usr/bin/env python
'''

Scans the genome for genomic A stretches and prints them in bed format to stdout.

Usage:

    flag_genomic_As_KevinRoy_like.py -fa [genome.fa] -fi [genome.fa.fai] -l [word_len] -A [min As within word] -maxCT [max C+Ts within word] -up [extend upstream splice junction] -dn [extend downstream splice junction]

-fa ... genome in fasta format
-f ... index for fasta file
-l ... length of motif, ie word length (default = 6)
-A ... minimum number of A for word to be flagged (default = 4)
-maxCT ... maximum number of C and Ts for word to be flagged (default = 2)
-up ... extend each single nucleotide hit by this number upstream (default = 0)
-dn ... extend each single nucleotide hit by this number downstream (default = 0)
-q ... quiet, ie suppress print of progression to stdout

Intended for combination as input for bedtools subtract of A-tail based 3p end sequencing.

Example usage:
python flag_genomic_As_KevinRoy_like_spliced.py -fa /Users/au220280/Documents/genomewide_datasets/annotations/hg38/Genome/GRCh38.fa -fi /Users/au220280/Documents/genomewide_datasets/annotations/hg38/Genome/GRCh38.fa.fai -b /Users/au220280/Documents/genomewide_datasets/annotations/hg38/RefSeq_GRCh38/RefSeqNCBIAll_GRCh38.bed -l 6 -A 4 -maxCT 0 -up 20 -dn 20 -q > Refseq_spliced_flagA4in6.bed

python flag_genomic_As_KevinRoy_like_spliced.py -fa /Users/au220280/Documents/genomewide_datasets/annotations/hg38/Genome/GRCh38.fa -fi /Users/au220280/Documents/genomewide_datasets/annotations/hg38/Genome/GRCh38.fa.fai -b /Users/au220280/Documents/genomewide_datasets/annotations/hg38/RefSeq_GRCh38/RefSeqNCBIAll_GRCh38.bed -l 18 -A 12 -maxCT 6 -up 20 -dn 20 -q > Refseq_spliced_flagA12in18.bed

cat Refseq_spliced_flagA4in6.bed Refseq_spliced_flagA12in18.bed | sort -k1,1 -k2,2n > Refseq_spliced_flagA4in6_and_A12in18.bed

'''


__author__ = 'schmidm'

import sys
import argparse


DNA_RC = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-fa", "--fasta", action='store', type=str, help="genome fasta file")
    parser.add_argument("-fi", "--faidx", action='store', type=str, help="index for fasta file")
    parser.add_argument("-b", "--bed", action='store', help="annotation file for splice junctions")
    parser.add_argument("-l", "--word_len", action='store', help="word length", type=int, default=5)
    parser.add_argument("-A", "--A", action='store', help="minimum As within word", type=int, default=0)
    parser.add_argument("-maxCT", "--maxCT", action='store', help="maximum C plus T within word", type=int, default=0)
    parser.add_argument("-up", "--upstreamjunction", action='store', help="scan in sequences up nts upstream exon-exon junctions", type=int, default=15)
    parser.add_argument("-dn", "--downstreamjunction", action='store', help="scan in sequences up nts downstream exon-exon junctions", type=int, default=15)
    parser.add_argument("-q", "--quiet", action='store_true', default=True, help="quiet, suppress print of progress info")

    args=parser.parse_args()
    return args

class Sequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq.upper()
        
        
class Transcript:
    '''
    >>> t=Transcript("\t".join(['chr1', '10', '100', 'x', '0', '+', '10', '100', '0', '3', '10,50,200', '0,20,100']))
    >>> t.exon_sizes
    [10, 50, 200]
    '''
    def __init__(self, line):
        line = line.split()
        self.chr = line[0]
        self.start = int(line[1])
        self.end = int(line[2])
        self.name = line[3]
        self.strand = line[5]
        self.exon_sizes = [int(s) for s in line[10].split(',')]
        self.exon_starts = [int(s) for s in line[11].split(',')]
        

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

def read_faidx(fname):
    with open(fname, 'r') as fi:
        lines = [line.rstrip().split('\t') for line in fi]
        faidx = {l[0]: [int(li) for li in l[1:]] for l in lines}
    return faidx


def read_bed(fname):
    '''
    :param fname: bed12 file
    :return: Dict of Transcript objects with key: chr name and value list of transcripts for that chromosome
    '''
    tr_dict = {}
    try:
        with open(fname, 'r') as f:
            for line in f:
                tr = Transcript(line.rstrip())
                try:
                    tr_dict[tr.chr].append(tr)
                except:
                    tr_dict[tr.chr] = [tr]
            
    except Exception as e:
        print 'ERROR: bedfile ', fname, 'not found or something else happened!'

    return tr_dict


def junctions_in_bed(fname):
    '''
    :param fname: bed12 file
    :return: Dict with chr names as keys with values as list of exon exon junction duples, ie (up_exon_end, dn_exon_start)
    >>> bedfile = "/Users/au220280/Documents/genomewide_datasets/annotations/hg38/RefSeq_GRCh38/RefSeqNCBIAll_GRCh38.bed"
    >>> junctions = junctions_in_bed(bedfile)
    '''
    junction_dict = {'+':{}, '-':{}}
    try:
        with open(fname, 'r') as f:
            for line in f:
                chr, strand, junctions = junctions_in_bedline(line)
                try:
                    junction_dict[strand][chr].extend(junctions)
                except:
                    junction_dict[strand][chr] = junctions
    
    except Exception as e:
        print 'ERROR: bedfile ', fname, 'not found or something else happened!'

    for strand in junction_dict.keys():
        for chr in junction_dict[strand].keys():
            junction_dict[strand][chr] = list(set(junction_dict[strand][chr]))
            
    return junction_dict


def junctions_in_bedline(bedline):
    '''
    :param fname: bed12 file line as string
    :return: Dict with chr names as keys with values as list of exon exon junction duples, ie (up_exon_end, dn_exon_start)
    >>> bedline = "\t".join(['chr1', '10', '100', 'x', '0', '+', '10', '100', '0', '3', '10,50,200', '0,20,100'])
    >>> junctions_in_bedline(bedline)
    ('chr1', '+', [(10, 20), (70, 100)])
    '''
    tr = bedline.rstrip().split()
    chr, start, strand = tr[0], int(tr[1]), tr[5]
    exon_nr = int(tr[9])
    if exon_nr > 1:
        exon_sizes = [int(s) for s in tr[10].split(',') if s != '']
        exon_starts = [int(s) for s in tr[11].split(',') if s != '']
        if strand == '-':
            start -= 1
        junctions = [ (start+exon_size+exon_starts[i], start+exon_starts[i+1]) for i, exon_size in enumerate(exon_sizes[:-1])]
    else:
        junctions = []
        
    return chr, strand, junctions


def junction_seq(chr, strand, junction, upstream, downstream, fasta, faidx):
    '''
    :param chr:
    :param strand:
    :param junction:
    :param upstream: upstream exon nt before junction to include
    :param downstream: downstream exon nt after junction to include
    :return: raw sequence
    '''
    if strand == '-':
        seq_up = get_seq(chr, junction[1], junction[1] + upstream, strand, faidx, fasta)
        seq_dn = get_seq(chr, junction[0] - downstream, junction[0], strand, faidx, fasta)
    else:
        seq_up = get_seq(chr, junction[0]-upstream, junction[0], strand, faidx, fasta)
        seq_dn = get_seq(chr, junction[1], junction[1] + downstream, strand, faidx, fasta)
        
    return seq_up + seq_dn
    
    
def reverse_complement(seq):
    return ''.join(REVCOM[c] for c in reversed(seq))
    

def flag_junctions(tr, seq, upstream, downstream, word_len, A, maxCT):
    '''
    list of flags over exon exon junctions but in genomic coordingates
    :param tr: Transcript
    :param seq: Sequence
    :param upstream: nt upstream to add
    :param downstream: nt downstream to add
    :return: list of sequences to flag
    >>> t=Transcript("\t".join(['chr1', '3', '30', 'x', '0', '+', '10', '100', '0', '3', '3,5,5,6', '0,5,12,20']))
    >>> s=Sequence("chr1", 'tttAAAttCCAAAttGAGGAtttAAAATTttt')
    >>> flag_junctions(t, s, 3, 3, 6, 4, 0)
    [9, 16]
    >>> tm=Transcript("\t".join(['chr1', '3', '30', 'x', '0', '-', '10', '100', '0', '3', '3,5,5,6', '0,5,12,20']))
    >>> sm=Sequence("chr1", 'tttTTTttGGTTTttCTCCTtttTTTTAAttt')
    >>> flag_junctions(tm, sm, 3, 3, 6, 4, 0)
    [18, 26]
    '''

    flagged_positions = []
    for i, ee in enumerate(tr.exon_starts[1:]):
        up_exon_end = tr.start + tr.exon_starts[i] + tr.exon_sizes[i]
        up_exon_start = up_exon_end - upstream if tr.strand == '+' else up_exon_end - downstream
        dn_exon_start = tr.start + ee
        dn_exon_end = dn_exon_start + downstream if tr.strand == '+' else dn_exon_start + upstream
        junction_seq = seq.seq[up_exon_start:up_exon_end] + seq.seq[dn_exon_start:dn_exon_end]

        if tr.strand == '-':
            junction_seq = ''.join(DNA_RC[s] for s in junction_seq[::-1])

        flagged_pos = flagged_in_seq(junction_seq, word_len, A, maxCT)

        if tr.strand == '-':
            flagged_positions.extend(dn_exon_end - f - 1 if f <= upstream else up_exon_end - f - 1for f in flagged_pos)
        else:
            flagged_positions.extend(f+up_exon_start if f <= upstream else f+dn_exon_start for f in flagged_pos )
        
    return flagged_positions


def flagged_in_seq(seq, word_len, A, maxCT):
    '''
    :param seq: seq object with name and seq field
    :param args: parsed arguments
    :return: flagged regions list for this sequence
    >>> flagged_in_seq('CGTAAAAAA', 6, 4, 0)
    [2]
    >>> flagged_in_seq('CGTAAAAAAGG', 6, 4, 0)
    [2, 3, 4]
    '''
    flagged_regions = []
    
    start = 0
    end = start + word_len
    cnts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    cnts['A'] = sum(1 for c in seq[start:end] if c == 'A')
    cnts['C'] = sum(1 for c in seq[start:end] if c == 'C')
    cnts['G'] = sum(1 for c in seq[start:end] if c == 'G')
    cnts['T'] = sum(1 for c in seq[start:end] if c == 'T')
    
    while end < len(seq):
        if cnts['A'] >= A and (cnts['C'] + cnts['T']) <= maxCT:
            flagged_regions.append(start - 1)
        
        cnts[seq[start]] -= 1
        cnts[seq[end]] += 1
        
        start += 1
        end += 1
    
    if cnts['A'] >= A and (cnts['C'] + cnts['T']) <= maxCT:
        flagged_regions.append(start - 1)
    
    return flagged_regions


def get_seq(chr, start, end, strand, faidx, fasta):
    if strand == '-':
        start, end = start + 1, end + 1
    width = end - start
    
    cidx = faidx[chr]
    chr_len, chr_start, seq_per_line_len, line_len = cidx[0], cidx[1], cidx[2], cidx[3]
    
    linebreaks_before_start = int(start / seq_per_line_len)
    linebreaks_before_end = int(end / seq_per_line_len)
    
    start_pos = chr_start + start + linebreaks_before_start
    chunk_size = width + linebreaks_before_end - linebreaks_before_start
    
    fasta.seek(start_pos)
    seq = fasta.read(chunk_size)
    seq = seq.replace('\n', '').upper()
    
    if strand == '-':
        seq = ''.join(DNA_RC[s] for s in seq[::-1])
    
    return seq


if __name__ == '__main__':

    args = parse_args()
    word_len, A, maxCT = args.word_len, args.A, args.maxCT
    upstream, downstream = args.upstreamjunction, args.downstreamjunction
    faidx_file = args.faidx
    fasta_file = args.fasta
    bedfile = args.bed
    
    #faidx_file = "/Users/au220280/Documents/genomewide_datasets/annotations/hg38/Genome/GRCh38.fa.fai"
    #fasta_file = "/Users/au220280/Documents/genomewide_datasets/annotations/hg38/Genome/GRCh38.fa"
    #bedfile = "/Users/au220280/Documents/genomewide_datasets/annotations/hg38/RefSeq_GRCh38/RefSeqNCBIAll_GRCh38.bed"
    #word_len, A, maxCT = 6, 4, 0
    #upstream, downstream = 20,20
    
    faidx = read_faidx(faidx_file)
    junctions = junctions_in_bed(bedfile)
    with open(fasta_file, 'r') as fa:
        for strand in junctions.keys():
            for chr in junctions[strand].keys():
                if chr not in faidx: continue
                for junction_list in junctions[strand][chr]:
                    if type(junction_list) == tuple:
                        junction_list = [junction_list]
                    for junction in junction_list:
                        seq = junction_seq(chr, strand, junction, upstream, downstream, fa, faidx)
                        for flag in flagged_in_seq(seq, word_len, A, maxCT):
                            if strand == '-':
                                if flag > upstream:
                                    pos = junction[0] - (flag - upstream)
                                else:
                                    pos = junction[1] + upstream - flag
                            else:
                                if flag > upstream:
                                    pos = flag + junction[1] - upstream
                                else:
                                    pos = flag + junction[0] - upstream
                            
                            print chr + "\t" + str(pos) + "\t" + str(pos + 1) + "\t" + seq + ":" + str(flag) +"\t0\t" + strand
        
                