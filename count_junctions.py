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

MIN_OVERLAP = 2 #minimum amount of bases on each side of a junction that need to be present for a call



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
bed_lines = [ INTERVAL(line.rstrip().split( "\t" )) for line in open( bed_file ) ]
bed_dict = {}
for line in open( bed_file ):
    iv = INTERVAL(line.rstrip().split( "\t" ))
    try:
        bed_dict[iv.chr][iv.strand].append(iv)
    except:
        bed_dict[iv.chr] = {'+': [], '-': []}
        bed_dict[iv.chr][iv.strand] = [iv]


##load the bam file
almnt_file = HTSeq.BAM_Reader( bam_file )


#count up the junctions reads
counts = collections.Counter( )
aligned_reads = (almnt for almnt in almnt_file if almnt.aligned)

unaligned_count = 0
outside_range_count = 0
i = 0
bed_i_plus = 0
bed_i_minus = 0
bed_chr = ''
bed_strand = ''

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

    #if iv.chrom == 'chrIII': break

    # get overlapping introns if any
    if iv.chrom != bed_chr:
        try:
            plus_bed_ar = bed_dict[iv.chrom]['+']
            minus_bed_ar = bed_dict[iv.chrom]['-']
            bed_chr = iv.chrom
            plus_bed_i = 0
            minus_bed_i = 0
             #print 'switching to chr ', bed_chr, ' strand', bed_strand
            # #print plus_bed_ar
            #print minus_bed_ar
        except:
             #print 'failed to find chr or strand', iv.chrom, iv.strand
            continue

    #note on those bams we are reverse complement to genome
    if iv.strand == '-':
        bed_ar = plus_bed_ar
        bed_i = plus_bed_i
    elif iv.strand == '+':
        bed_ar = minus_bed_ar
        bed_i = minus_bed_i
    else:
         #print 'bam entry missing with strand not + or - !!'
        continue

    if bed_i >= len(bed_ar):
        bed_i = len(bed_ar) - 1

    # move backward to before first potential overlapping intron
    min_start = (iv.start - MIN_OVERLAP)
    while 0 < bed_i and min_start > bed_ar[bed_i].end:
        bed_i -= 1

    #print 'backward to: ', bed_i #, str(bed_ar[bed_i])
    # move forward to first potential overlapping intron
    while bed_i < len(bed_ar) and iv.end < bed_ar[bed_i].start:
        bed_i += 1

    #print 'forward to: ', bed_i
    if bed_i >= len(bed_ar):
        #print 'out of luck for interval range', iv.start, iv.end
        #break
        continue

    # loop through overlaps
    while bed_i < len(bed_ar) and iv.start < bed_ar[bed_i].end:
        interval = bed_ar[bed_i]
        #print 'considering ', iv.start, '-', iv.end, 'for interval: ', interval.start, '-', interval.end
        #categorize the interval
        if iv.start < interval.end:
            #upstream or overlapping
            if iv.start < (interval.start - MIN_OVERLAP) and iv.end > interval.start:
                    #alignment starts upstream and ends internal or downstream

                    if (interval.start + MIN_OVERLAP) < iv.end < interval.end:
                        #alignments ends internal --> 1st junction overlap
                        #print 'added EI'
                        interval.EI += 1
                    elif iv.end > (interval.end + MIN_OVERLAP):
                        #potential exon-exon
                        intron_found = False
                        #make sure there is no internal alignment to intron
                        for cigop in aligned_read.cigar:
                            if cigop.type == "M" and interval.overlaps_internal(cigop.ref_iv):
                                intron_found = True
                                 #print 'intronic in EE candidate'
                                break
                        if not intron_found:
                             #print 'added EE'
                            interval.EE += 1
                            #if interval.id == 'YFL039C':
                            #  print iv
                        else:
                            #print 'added ambigous'
                            interval.ambigous += 1
                    else:
                        #ambigous
                        #print 'added ambigous'
                        interval.ambigous += 1

            elif interval.start < iv.start < (interval.end - MIN_OVERLAP):
                    # alignment starts inside
                    if iv.end < interval.end:
                        # alignments ends internal --> intronic
                        #print 'added intronic'
                        interval.intronic += 1
                    elif iv.end > (interval.end + MIN_OVERLAP):
                        # intron-exon
                        #print 'added IE', interval.id
                        interval.IE += 1
                    else:
                        # ambigous
                        #print 'added ambigous'
                        interval.ambigous += 1

            else:
                outside_range_count += 1

        bed_i += 1
        if iv.strand == '+':
            plus_bed_i = bed_i
        else:
            minus_bed_i = bed_i

# print 'bam file: ', bam_file
# print '  interval file: ', bed_file
# print '  total reads in bam file: ', i
# print '  not aligned :', unaligned_count
# print '  overlap bed intervals: ', (i - outside_range_count)

print 'chr', '\t', 'start', '\t', 'end', '\t', 'id', '\t', 'score', '\t', 'strand', '\t', 'EE', '\t', 'IE', '\t', 'EI', '\t', 'intronic', '\t','ambigous'

for chr in bed_dict:
    for strand in ['+', '-']:
        for iv in bed_dict[chr][strand]:
            print iv
