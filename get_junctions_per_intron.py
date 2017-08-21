#!/usr/bin/env python
'''

Get reads for each exon-exon, exon-intron and intron-exon junctions for all introns in a BED file

Usage: get_junctions_per_intron.py bamfile introns_bedfile MIN_BP_OVERLAP

MIN_BP_OVERLAP minimum amount of bases on each side of a junction that need to be present for a call (default=5)

For this way of counting its easiest to provide a bed file for the introns.
The bam file should be sorted.

writes a sam file with the junctions for EI, IE and EE

where EI is exon-intron and IE is intron-exon. Hence strand information is not taken into account here.


example:
'''

__author__ = 'schmidm'

import sys
import collections
import HTSeq

bam_file = sys.argv[1]
bed_file = sys.argv[2]

if len(sys.argv) >= 4:
    MIN_OVERLAP = int(sys.argv[3]) #minimum amount of bases on each side of a junction that need to be present for a call
else:
    MIN_OVERLAP = 5

outputfile = bam_file.replace('.bam', '_junctions_downstream_exon_overlap_length.txt')


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



def get_junction_reads(aligned_reads, bed_lines):
    i = 0
    bed_i = 0
    unaligned_count = 0
    outside_range_count = 0

    EE = []
    EI = []
    IE = []

    interval = INTERVAL(bed_lines.next())
    strand = interval.strand
    interval_dict = {}
    interval_dict[interval] = {'EE':[], 'EI':[], 'IE':[]}

    for aligned_read in aligned_reads:
        if not aligned_read.aligned:
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
             #print str(interval)
             try:  # wrap all in a try / except
                 interval = INTERVAL(bed_lines.next())
                 interval_dict[interval] = {'EE': [], 'EI': [], 'IE': []}
             except StopIteration, e:  # make the process on the last element here
                 return(interval_dict)


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
                        interval_dict[interval]['EI'].append(aligned_read)
                    elif iv.end > (interval.end + MIN_OVERLAP):
                        #potential exon-exon
                        intron_found = False
                        #make sure there is no internal alignment to intron
                        for cigop in aligned_read.cigar:
                            if cigop.type == "M" and interval.overlaps_internal(cigop.ref_iv):
                                intron_found = True
                                break
                        if intron_found:
                            if strand == '+':
                                interval_dict[interval]['EI'].append(aligned_read)
                            else:
                                interval_dict[interval]['IE'].append(aligned_read)
                        else:
                            interval_dict[interval]['EE'].append(aligned_read)
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
                        if strand == '+':
                            interval_dict[interval]['IE'].append(aligned_read)
                        else:
                            interval_dict[interval]['EI'].append(aligned_read)
                    else:
                        # ambigous
                        interval.ambigous += 1

            else:
                outside_range_count += 1

    return( interval_dict )

# print 'bam file: ', bam_file
# print '  interval file: ', bed_file
# print '  total reads in bam file: ', i
# print '  not aligned :', unaligned_count
# print '  overlap bed intervals: ', (i - outside_range_count)

##parse the bed file
features =  HTSeq.GenomicArrayOfSets( "auto", stranded=True )
bed_lines = [ line.rstrip().split( "\t" ) for line in open( bed_file ) ]


##load the bam file
almnt_file = HTSeq.BAM_Reader( bam_file )


aligned_reads_plus = (almnt for almnt in almnt_file if almnt.aligned and almnt.iv.strand == '+')
aligned_reads_minus = (almnt for almnt in almnt_file if almnt.aligned and almnt.iv.strand == '-')

bed_plus = (line for line in bed_lines if line[5] == '+')
bed_minus = (line for line in bed_lines if line[5] == '-')

plus_intron_reads = get_junction_reads(aligned_reads_minus, bed_plus)

minus_intron_reads = get_junction_reads(aligned_reads_plus, bed_minus)


f = open(outputfile, 'w')


for k, v in plus_intron_reads.items():
    # print k.chr, k.start, k.end, k.id, k.strand
    #
    # for e in v['EE']:
    #     print '  EE ', e.iv.start, ' ', e.iv.end
    #
    # for e in v['EI']:
    #     print '  EI ', e.iv.start, ' ', e.iv.end
    #
    # for e in v['IE']:
    #     print '  IE ', e.iv.start, ' ', e.iv.end
    for e in v['IE']:
        if k.strand == '+':
            f.write('IE', e.iv.end - k.end)
        else:
            f.write('IE', k.start - e.iv.start)

    for e in v['EE']:
        if k.strand == '+':
            f.write('EE', e.iv.end - k.end)
        else:
            f.write('EE', k.start - e.iv.start)

