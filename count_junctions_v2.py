'''

Counts reads for each exon-exon, exon-intron and intron-exon junctions for all introns in a BED12 file

Usage: count_junctions.py bamfile introns_bedfile [--min_overlap, --expand_size]

The bam file needs to be sorted and have an index file. The bed file does not need to be sorted

writes simply the 6-column bed file adding columns EE_count, SD_count, SA_count, intronic, ambigous counts as additional columns

where
. EE is exon-exon junction without intron coverage

. SD is splice-donor ie 5'SS (exon-intron on plus strand and intron-exon on minus strand)

. SA is splice-acceptor ie 3'SS (intron-exon on plus strand and exon-intron on minus strand).


'''

__author__ = 'schmidm'



import HTSeq
import argparse


MIN_OVERLAP = 0 #minimum amount of bases on each side of a junction that need to be present for a call
WINDOW_SIZE = 100 # only consider reads with start or end within intron +/- WINDOW_SIZE usually read-length


class GENE:
	'''
	test
	>>> from collections import Counter
	>>> bed = "/Users/au220280/Documents/genomewide_datasets/annotations/sacCer3/sacCer3_ORFT_multiexonics.bed"
	>>> bedf = open(bed, 'r')
	>>> lines = [line.split() for line in bedf]
	>>> bedf.close()
	>>> bam_file = "/Users/au220280/Documents/Results/Lexogen_RNAseq_2/bams/42471_GGTATA_C9P6RANXX_5_20160808B_20160808_trimmed_cleanAligned.sortedByCoord.out.bam"
	>>> almnt_file = HTSeq.BAM_Reader( bam_file )
	>>> gene = [GENE(line) for line in lines if 'SNC1' in line[3]][0]
	>>> introns = gene.introns()
	>>> i1_reads = [read for read in almnt_file[introns[0]]]
	>>> Counter([gene.add_read(i1_read, introns[0]) for i1_read in i1_reads])
	Counter({'EE': 3, 'intronic': 3, 'ambigous': 1, None: 1})
	>>> gene = [GENE(line) for line in lines if 'MOB2' in line[3]][0]
	>>> introns = gene.introns()
	>>> i1_reads = [read for read in almnt_file[introns[0]]]
	>>> Counter([gene.add_read(i1_read, introns[0]) for i1_read in i1_reads])
	Counter({'SA': 3, 'intronic': 3, 'SD': 2, None: 1})
	>>> gene = [GENE(line) for line in lines if 'MCM21' in line[3]][0]
	>>> introns = gene.introns()
	>>> i1_reads = [read for read in almnt_file[introns[0]]]
	>>> Counter([gene.add_read(i1_read, introns[0]) for i1_read in i1_reads])
	Counter({'SA': 10, 'intronic': 8, 'SD': 3, 'EE': 3, None: 1})
	'''
	
	def __init__(self, line):
		self.chr = line[0]
		self.start = int(line[1])
		self.end = int(line[2])
		self.id = line[3]
		self.score = line[4]
		self.strand = line[5]
		self.n_exons = int(line[9])
		self.starts = [int(i) for i in line[11].split(',') if i != '']
		self.sizes = [int(i) for i in line[10].split(',') if i != '']
		self.EE = 0  # exon-exon splice junction reads
		self.SA = 0  # splice acceptor ie 5'SS exon-intron junction reads
		self.SD = 0  # splice donor ie 5'SS exon-intron junction reads
		self.ambigous = 0  # usually only those reads with ambigous overlaps
		self.intronic = 0  # fully intronic read
		self.exonic = 0  # fully exonic read
	
	
	def interval(self):
		return (HTSeq.GenomicInterval(self.chr, self.start, self.end))
			
	def introns(self):
		'''
		collects all introns for this BED interval
		'''
		read_strand = '+' if self.strand == '-' else '-'
		introns = [
			HTSeq.GenomicInterval(self.chr, self.start + self.starts[i] + self.sizes[i], self.start + self.starts[i+1],
			                      read_strand) for i in range(self.n_exons - 1)]
		
		return (introns)
	
	
	def count_junctions(self, almnt_file, min_overlap=2):
		'''
		counts all junctions
		:param almnt_file: HTSeq.BAM_Reader object
		:return: nothing, junction counts are added to instance
		'''
		for read in almnt_file[self.interval()]:
			if read.iv.strand == self.strand: #note reads are considered revcom to annotation, the intersection does not take strand into account
				continue
			for intron in self.introns():
				if self.add_read(read, intron, min_overlap): #only count each read once, even if it may overlap more than 1 intron
					break

			
	def add_read(self, read, intron, min_overlap=2):
		'''
		
		:param read: single HTSeq.BAM_Reader
		:param intron: HTSeq.GenomicInterval of an intron
		:return: nothing, junction counts are added to instance
		'''
		if read.iv.end < (intron.start+min_overlap) or read.iv.start >= (intron.end-min_overlap):
			#read fully up or downstream of intron
			return(None)
		elif read.iv.start < intron.end:
			# alignment starts upstream or inside intron
			if read.iv.start <= (intron.start - min_overlap) and read.iv.end >= intron.start:
				# alignment starts upstream and ends internal or downstream
				if (intron.start + min_overlap) <= read.iv.end <= intron.end:
					return(self.add_EI(read.iv))
				elif (intron.end + min_overlap) < read.iv.end:
					intron_found = False
					# make sure there is no internal alignment to intron
					for cigop in read.cigar:
						if cigop.type == "M" and self.overlaps_internal(cigop.ref_iv, intron):
							intron_found = True
							break
					if not intron_found:
						self.EE += 1
						return('EE')
					else:
						self.ambigous += 1
						return('ambigous')
				else:
					self.ambigous += 1
					return('ambigous')
			
			elif intron.start <= read.iv.start < (intron.end - min_overlap):
				# alignment starts inside
				if read.iv.end <= intron.end:
					# alignments ends internal --> intronic
					self.intronic += 1
					return('intronic')
				elif (intron.end + min_overlap) < read.iv.end:
					return(self.add_IE(read.iv))
				else:
					self.ambigous += 1
					return('ambigous')
		return(None)
	
	
	def add_IE(self, iv):
		if self.strand == '+':
			self.SA += 1
			return('SA')
		else:
			self.SD += 1
			return('SD')
	
	def add_EI(self, iv):
		if self.strand == '+':
			self.SD += 1
			return('SD')
		else:
			self.SA += 1
			return('SA')
	
	def overlaps_internal(self, iv1, iv2):
		''' checks whether interval iv1 has any overlap with iv2,
			ie whether there are intronic parts in the interval
			assuming they come from same sequence'''
		if iv1.end <= iv2.start or iv1.start >= iv2.end:
			return(False)
		return(True)
			
			
	def __str__(self):
		return('\t'.join([self.chr, str(self.start), str(self.end), self.id, self.score, self.strand, str(self.EE), str(self.SD), str(self.SA), str(self.intronic), str(self.ambigous)]))
	
	
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(usage=__doc__)
	parser.add_argument('bam_file')
	parser.add_argument('bed_file')
	parser.add_argument('--min_overlap', default=2, type=int,
	                    help="minimum base pair overlap required upstream and downstream of junction")
	parser.add_argument('--expand_size', default=100, type=int,
	                    help="expand intron by this size on each size in first pass filtering of reads and only consider those reads for classification")
	
	args = parser.parse_args()
	
	bam_file = args.bam_file
	bed_file = args.bed_file
	MIN_OVERLAP = args.min_overlap  # minimum amount of bases on each side of a junction that need to be present for a call
	WINDOW_SIZE = args.expand_size  # deprecated ... only consider reads with start or end within intron +/- WINDOW_SIZE usually read-length
	
	#bam = "/Users/schmidm/Documents/Results/Lexogen_RNAseq_2/bams/42471_GGTATA_C9P6RANXX_5_20160808B_20160808_trimmed_cleanAligned.sortedByCoord.out.bam"
	
	#bed = "/Users/au220280/Documents/genomewide_datasets/annotations/sacCer3/sacCer3_ORFT_multiexonic.bed"
	
	##load the bam file
	bam_reader = HTSeq.BAM_Reader(bam_file)
	
	
	print('chr', '\t', 'start', '\t', 'end', '\t', 'id', '\t', 'score', '\t', 'strand', '\t', 'EE', '\t', 'SD', '\t', 'SA', '\t', 'intronic', '\t', 'ambigous')
	
	with open(bed_file, 'r') as bed:
		for line in bed:
			iv = GENE(line.rstrip().split("\t"))
			iv.count_junctions(bam_reader, MIN_OVERLAP)
			print(str(iv))

