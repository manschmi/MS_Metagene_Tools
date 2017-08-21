#!/usr/bin/env python
'''

Counts reads for each exon-exon, exon-intron and intron-exon junctions for all introns in a BED file

Usage: count_junctions.py bamfile introns_bedfile [--min_overlap, --expand_size, --reportOverlapLength]

For this way of counting its easiest to provide a bed file for the introns.
The bam file needs to be sorted and have an index file. The bed file does not need to be sorted

writes simply the 6-column bed file adding columns EE_count, SD_count, SA_count, intronic, ambigous counts as additional columns

where
. SD is splice-donor ie 5'SS (EI on plus strand and IE on minus strand)

. SA is splice-acceptor ie 3'SS (IE on plus strand and EI on minus strand).

if --reportOverlapLength is TRUE:
adds additional columns consisting of comma-separated lists for length of EE_up EE_dn SD_up SD_dn SA_up SA_dn overlaps
'''

__author__ = 'schmidm'

import sys
import argparse
import HTSeq


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('bam_file')
parser.add_argument('bed_file')
parser.add_argument('--min_overlap', default=2, type=int, help="minimum base pair overlap required upstream and downstream of junction")
parser.add_argument('--expand_size', default=100, type=int, help="expand intron by this size on each size in first pass filtering of reads and only consider those reads for classification")
parser.add_argument('--reportOverlapLength', default=False, action="store_true", help="report length of upstream and downstream overlaps for each overlap type")
args = parser.parse_args()


bam_file = sys.argv[1]
bed_file = sys.argv[2]

MIN_OVERLAP = args.min_overlap #minimum amount of bases on each side of a junction that need to be present for a call
WINDOW_SIZE = args.expand_size # only consider reads with start or end within intron +/- WINDOW_SIZE usually read-length


class INTERVAL:

    def __init__(self, line):
        self.chr = line[0]
        self.start = int(line[1])
        self.end = int(line[2])
        self.id = line[3]
        self.score = line[4]
        self.strand = line[5]
        self.EE = [] #exon-exon splice junction reads
        self.SA = [] #splice acceptor ie 5'SS exon-intron junction reads
        self.SD = [] #splice donor ie 5'SS exon-intron junction reads
        self.ambigous = 0
        self.intronic = 0

    def add_read(self, iv):
        if iv.start < self.end:
            # upstream or overlapping
            if (self.start - WINDOW_SIZE) < iv.start < (self.start - MIN_OVERLAP) and iv.end > self.start:
                # alignment starts upstream and ends internal or downstream
                if (self.start + MIN_OVERLAP) < iv.end < self.end:
                    self.add_EI(iv)
                elif (self.end + MIN_OVERLAP) < iv.end < (self.end + WINDOW_SIZE) :
                    self_found = False
                    # make sure there is no internal alignment to intron
                    for cigop in aligned_read.cigar:
                        if cigop.type == "M" and self.overlaps_internal(cigop.ref_iv):
                            self_found = True
                            break
                    if not self_found:
                        self.add_EE(iv)
                    else:
                        self.ambigous += 1
                else:
                    self.ambigous += 1

            elif self.start < iv.start < (self.end - MIN_OVERLAP):
                # alignment starts inside
                if iv.end < self.end:
                    # alignments ends internal --> intronic
                    self.intronic += 1
                elif (self.end + MIN_OVERLAP) < iv.end < (self.end + WINDOW_SIZE) :
                    self.add_IE(iv)
                else:
                    self.ambigous += 1


    def add_EE(self, iv):
        self.EE.append(iv)
        #self.EE += 1
        #self.EEdists.append((self.start - iv.start, iv.end - self.end))

    def add_IE(self, iv):
        if self.strand == '+':
            self.SA.append(iv)
            #self.SA += 1
            #self.SAdists.append((self.end - iv.start, iv.end - self.end))
        else:
            self.SD.append(iv)
            #self.SD += 1
            #self.SDdists.append((self.start - iv.start, iv.end - self.start))

    def add_EI(self, iv):
        if self.strand == '+':
            self.SD.append(iv)
            #self.SA += 1
            #self.SAdists.append((self.end - iv.start, iv.end - self.end))
        else:
            self.SA.append(iv)
            #self.SD += 1
            #self.SDdists.append((self.start - iv.start, iv.end - self.start))

    def overlaps_internal(self, iv):
        ''' checks whether interval has any overlap with the bedinterval,
            ie whether there are intronic parts in the interval '''
        if iv.start < self.start and iv.end > self.start:
            return True
        elif iv.start < self.end and iv.end > self.end:
            return True
        elif iv.start > self.start and iv.end < self.end:
            return True
        return False

    def summary_str(self):
        '''
        :return: summary of object, ie counts per junction as string
        '''
        return '\t'.join([self.chr, str(self.start), str(self.end), self.id, self.score, self.strand, str(len(self.EE)), str(len(self.SD)), str(len(self.SA)), str(self.intronic), str(self.ambigous)])

    def overlaps_as_str(self):
        '''
        :return: tab separate list of comma-separated lists for upstream and downstream overlaps for EE, SD, SA
        '''

        if self.strand == '+':
            ee_up = ','.join(str(self.start - ee.start) for ee in self.EE)
            ee_dn = ','.join(str(ee.end - self.end) for ee in self.EE)
            sd_up = ','.join(str(self.start - ee.start) for ee in self.SD)
            sd_dn = ','.join(str(ee.end - self.start) for ee in self.SD)
            sa_up = ','.join(str(self.end - ee.start) for ee in self.SA)
            sa_dn = ','.join(str(ee.end - self.end) for ee in self.SA)

        else:
            sa_dn = ','.join(str(self.start - ee.start) for ee in self.SA)
            sa_up = ','.join(str(ee.end - self.start) for ee in self.SA)
            sd_dn = ','.join(str(self.end - ee.start)  for ee in self.SD)
            sd_up = ','.join(str(ee.end - self.end) for ee in self.SD)
            ee_dn = ','.join(str(self.start - ee.start)  for ee in self.EE)
            ee_up = ','.join(str(ee.end - self.end) for ee in self.EE)

        return ee_up + '\t' + ee_dn + '\t' + sd_up + '\t' + sd_dn + '\t' + sa_up + '\t' + sa_dn


##parse the bed file
bed_lines = [ INTERVAL(line.rstrip().split( "\t" )) for line in open( bed_file ) ]

##load the bam file
almnt_file = HTSeq.BAM_Reader( bam_file )

unaligned_count = 0
outside_range_count = 0
i = 0

if not args.reportOverlapLength:
    print 'chr', '\t', 'start', '\t', 'end', '\t', 'id', '\t', 'score', '\t', 'strand', '\t', 'EE', '\t', 'SD', '\t', 'SA', '\t', 'intronic', '\t','ambigous'
else:
    print 'chr', '\t', 'start', '\t', 'end', '\t', 'id', '\t', 'score', '\t', 'strand', '\t', 'EE', '\t', 'SD', '\t', 'SA', '\t', 'intronic', '\t', 'ambigous', '\t', 'EE_up', '\t', 'EE_dn',' \t', 'SD_up', '\t', 'SD_dn', '\t', 'SA_up', '\t', 'SA_dn'

for intron in bed_lines:
    if intron.strand == '+':
        read_strand = '-'
    else:
        read_strand = '+'
    intron_window = HTSeq.GenomicInterval( intron.chr, intron.start - WINDOW_SIZE, intron.end + WINDOW_SIZE, read_strand )

    for aligned_read in almnt_file[intron_window]:
        if not aligned_read.aligned:
            unaligned_count += 1
            continue
        iv = aligned_read.iv

        if not aligned_read.optional_field('NH') == 1:
            continue

        if not iv.strand == read_strand:
            #not: almnt_file[intron_window] ignores the strand info even though its there
            continue

        intron.add_read(iv)

    if not args.reportOverlapLength:
        print intron.summary_str()
    else:
        print intron.summary_str() + '\t' + intron.overlaps_as_str()
