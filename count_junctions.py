#!/usr/bin/env python
'''

Counts reads for each exon-exon, exon-intron and intron-exon junctions for all introns in a BED file

Usage: count_junctions.py bamfile introns_bedfile

For this way of counting its easiest to provide a bed file for the introns.
The bam file needs to be sorted and have an index file. The bed file does not need to be sorted

writes simply a 6-column bed file with, EE_count, EI_count, IE_count, intronic, ambigous counts as additional columns

where EI is exon-intron and IE is intron-exon. The strand information is not taken into account here, an EI is a 5SS on + strand and 3'ss on plus strand.

'''

__author__ = 'schmidm'

import sys
import collections
import HTSeq

bam_file = sys.argv[1]
bed_file = sys.argv[2]

MIN_OVERLAP = 2 #minimum amount of bases on each side of a junction that need to be present for a call
WINDOW_SIZE = 100 # only consider reads with start or end within intron +/- WINDOW_SIZE


class INTERVAL:

    def __init__(self, line):
        self.chr = line[0]
        self.start = int(line[1])
        self.end = int(line[2])
        self.id = line[3]
        self.score = line[4]
        self.strand = line[5]
        self.EE = 0
        self.EI = 0
        self.IE = 0
        self.ambigous = 0
        self.intronic = 0

    def __str__(self):
        return '\t'.join([self.chr, str(self.start), str(self.end), self.id, self.score, self.strand, str(self.EE), str(self.IE), str(self.EI), str(self.intronic), str(self.ambigous)])

    def overlaps_internal(self, interval):
        ''' checks whether interval has any overlap with the bedinterval,
            ie whether there are intronic parts in the interval '''
        if interval.start < self.start and interval.end > self.start:
            return True
        elif interval.start < self.end and interval.end > self.end:
            return True
        elif interval.start > self.start and interval.end < self.end:
            return True
        return False


##parse the bed file
bed_lines = [ INTERVAL(line.rstrip().split( "\t" )) for line in open( bed_file ) ]

##load the bam file
almnt_file = HTSeq.BAM_Reader( bam_file )

unaligned_count = 0
outside_range_count = 0
i = 0

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

        if iv.start < intron.end:
            # upstream or overlapping
            if iv.start < (intron.start - MIN_OVERLAP) and iv.end > intron.start:
                # alignment starts upstream and ends internal or downstream

                if (intron.start + MIN_OVERLAP) < iv.end < intron.end:
                    # alignments ends internal --> 1st junction overlap
                    # print 'added EI'
                    intron.EI += 1
                elif iv.end > (intron.end + MIN_OVERLAP):
                    # potential exon-exon
                    intron_found = False
                    # make sure there is no internal alignment to intron
                    for cigop in aligned_read.cigar:
                        if cigop.type == "M" and intron.overlaps_internal(cigop.ref_iv):
                            intron_found = True
                            # print 'intronic in EE candidate'
                            break
                    if not intron_found:
                        # print 'added EE'
                        intron.EE += 1
                        # if intron.id == 'YFL039C':
                        #  print iv
                    else:
                        # print 'added ambigous'
                        intron.ambigous += 1
                else:
                    # ambigous
                    # print 'added ambigous'
                    intron.ambigous += 1

            elif intron.start < iv.start < (intron.end - MIN_OVERLAP):
                # alignment starts inside
                if iv.end < intron.end:
                    # alignments ends internal --> intronic
                    # print 'added intronic'
                    intron.intronic += 1
                elif iv.end > (intron.end + MIN_OVERLAP):
                    # intron-exon
                    # print 'added IE', intron.id
                    intron.IE += 1
                else:
                    # ambigous
                    # print 'added ambigous'
                    intron.ambigous += 1

            else:
                outside_range_count += 1



print 'chr', '\t', 'start', '\t', 'end', '\t', 'id', '\t', 'score', '\t', 'strand', '\t', 'EE', '\t', 'IE', '\t', 'EI', '\t', 'intronic', '\t','ambigous'
for intron in bed_lines:
    print intron
