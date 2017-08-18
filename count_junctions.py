#!/usr/bin/env python
'''

Counts reads for each exon-exon, exon-intron and intron-exon junctions for all introns in a BED file

Usage: count_junctions.py bamfile introns_bedfile

For this way of counting its easiest to provide a bed file for the introns.
The bam file should be sorted.

writes simply a 6-column bed file with, EE_count, EI_count, IE_count, intronic, ambigous counts as additional columns

where EI is exon-intron and IE is intron-exon. The strand information is not taken into account here, an EI is a 5SS on + strand and 3'ss on plus strand.

'''

__author__ = 'schmidm'

import sys
import collections
import HTSeq

bam_file = sys.argv[1]
bed_file = sys.argv[2]

MIN_OVERLAP = 5 #minimum amount of bases on each side of a junction that need to be present for a call



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
features =  HTSeq.GenomicArrayOfSets( "auto", stranded=True )
bed_lines = ( line.rstrip().split( "\t" ) for line in open( bed_file ) )


##load the bam file
almnt_file = HTSeq.BAM_Reader( bam_file )


#count up the junctions reads
counts = collections.Counter( )
aligned_reads = (almnt for almnt in almnt_file if almnt.aligned)

unaligned_count = 0
outside_range_count = 0
i = 0
bed_i = 0


interval = INTERVAL(bed_lines.next())


for aligned_read in aligned_reads:
    if not aligned_read.aligned:
        unaligned_count += 1
        continue
    iv = aligned_read.iv
    #print i, iv

    i += 1
    #if i % 1000 == 0:
    #    print i, iv
    #if i > 260000:
    #    break

    # loop through introns upstream
    while interval and (iv.chrom != interval.chr or iv.start > interval.end ):
         print str(interval)
         try:  # wrap all in a try / except
             interval = INTERVAL(bed_lines.next())
         except StopIteration, e:  # make the process on the last element here
             exit()


    #categorize the interval
    if iv.chrom == interval.chr and iv.start < interval.end:
        #upstream or overlapping

        if iv.strand == interval.strand:
            #wrong strand at least for those kind of RNAseq bam files
            continue

        if iv.start < (interval.start - MIN_OVERLAP) and iv.end > interval.start:
                #alignment starts upstream and ends internal or downstream

                if (interval.start + MIN_OVERLAP) < iv.end < interval.end:
                    #alignments ends internal --> 1st junction overlap
                    interval.EI += 1
                elif iv.end > (interval.end + MIN_OVERLAP):
                    #potential exon-exon
                    intron_found = False
                    #make sure there is no internal alignment to intron
                    for cigop in aligned_read.cigar:
                        if cigop.type == "M" and interval.overlaps_internal(cigop.ref_iv):
                            intron_found = True
                            break
                    if not intron_found:
                        interval.EE += 1
                    else:
                        interval.ambigous += 1
                else:
                    #ambigous
                    interval.ambigous += 1

        elif interval.start < iv.start < (interval.end - MIN_OVERLAP):
                # alignment starts inside
                if iv.end < interval.end:
                    # alignments ends internal --> intronic
                    interval.intronic += 1
                elif iv.end > (interval.end + MIN_OVERLAP):
                    # intron-exon
                    interval.IE += 1
                else:
                    # ambigous
                    interval.ambigous += 1

        else:
            outside_range_count += 1


# print 'bam file: ', bam_file
# print '  interval file: ', bed_file
# print '  total reads in bam file: ', i
# print '  not aligned :', unaligned_count
# print '  overlap bed intervals: ', (i - outside_range_count)