#!/usr/bin/env python
'''

Collects transcripts isoforms from GFF3 file and puts each isoform as line in bed12 file.

thick_start and thick_end (column 7,8) are from start_codon and stop_codon of gff3 files.

column4 is seq_id:geneID:transcriptID

output to stdout

'''

__author__ = 'Manfred Schmid ms@mbg.au.dk'

import sys
import argparse



parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('gff3_file')
args = parser.parse_args()

class GFF_Entry:

    def parse_attributes(self, gff_attr):
        '''
        :param gff_attr: attributes from gff line as single string
        :return: dict of named attributes
        '''
        attr = [att.split('=')for att in gff_attr.split(';')]
        return { att[0]:att[1] for att in attr }

    def __init__(self, line):
        '''
        uses a tab-split line from a gff3 file to build the object
        '''
        self.seq_id = line[0]
        self.source = line[1]
        self.type = line[2]
        self.start = int(line[3])
        self.end = int(line[4])
        self.score = line[5]
        self.strand = line[6]
        self.phase = line[7]
        self.attributes = self.parse_attributes(line[8])

    def bedify(self):
        '''
        convert from 1-based closed to 0-based half-open
        '''
        self.start -= 1
        return self


class Gene:
    def __init__(self, gff_entry):
        self.id = gff_entry.attributes['ID']
        self.entry = gff_entry
        self.transcripts = []

    def add_transcript(self, gff_entry):
        self.transcripts.append(gff_entry)

    def bed_line(self):
        return '\t'.join([self.entry.seq_id,
                          str(self.entry.start),
                          str(self.entry.end),
                          'gene:'+self.id.split(':')[1],
                          '0',
                          self.entry.strand])



class Transcript:
    def __init__(self, gff_entry):
        self.seq_id = gff_entry.seq_id
        self.start = gff_entry.start
        self.end = gff_entry.end
        self.strand = gff_entry.strand

        self.id = gff_entry.attributes['ID']
        self.parent = gff_entry.attributes['Parent']

        self.exons = [(gff_entry.start, gff_entry.end)]
        self.thick_start = gff_entry.start
        self.thick_end =  gff_entry.end

    def add_exon(self, gff_entry):
        self.exons.append((gff_entry.start, gff_entry.end))
        if gff_entry.start < self.start:
            self.start = gff_entry.start
        if gff_entry.end > self.end:
            self.end = gff_entry.end

    def set_thick_start(self, pos):
        #if not self.start <= pos < self.end:
        #    raise Exception('trying to set illegal thick start')
        self.thick_start = pos

    def set_thick_end(self, pos):
        #if not self.start < pos <= self.end:
        #    raise Exception('trying to set illegal thick end')
        self.thick_end = pos

    def bed_line(self):
        #Transcript object as string of a bed line
        exon_cnt = len(self.exons)
        exon_starts = [ exon[0]-self.start for exon in self.exons ]
        exon_sizes = [ exon[1]-exon[0] for exon in self.exons ]
        if self.strand == '-':
            exon_starts = reversed(exon_starts)
            exon_sizes = reversed(exon_sizes)
        exon_starts = ','.join(str(e) for e in exon_starts) + ','
        exon_sizes = ','.join(str(e) for e in exon_sizes) + ','
        return '\t'.join([self.seq_id,
                          str(self.start),
                          str(self.end),
                          'transcript:' + self.parent.split(':')[1] + ':' + self.id,
                          '0',
                          self.strand,
                          str(self.thick_start),
                          str(self.thick_end),
                          '0',
                          str(exon_cnt),
                          exon_sizes,
                          exon_starts])


#load the gff file and convert to bed intervals
f = open(args.gff3_file, 'r')
gff = [ GFF_Entry(line.split('\t')) for line in f ]
f.close()
bed = [ entry.bedify() for entry in gff ]

#load the genes into an array and dict
genes = [ Gene(line) for line in bed if line.type == 'gene' ]
genes_dict = { gene.id:gene for gene in genes }

#load the exons
exons = [ line for line in bed if line.type == 'exon' ]
transcripts_dict = {}

#combine exons for each transcript isoform into Transcripts and add them to the genes
for exon in exons:
    transcript_id = exon.attributes['ID']
    transcript_parent = exon.attributes['Parent']
    parent_gene = genes_dict[transcript_parent]
    if not transcript_id in transcripts_dict:
        transcript = Transcript(exon)
        transcripts_dict[transcript_id] = transcript
        parent_gene.add_transcript(transcript)
    else:
        transcripts_dict[transcript_id].add_exon(exon)


#start and stop codons
start_codons = [ line for line in gff if line.type == 'start_codon' ]
stop_codons = [ line for line in gff if line.type == 'stop_codon' ]
for start_codon in start_codons:
    transcript_id = start_codon.attributes['ID']
    if transcript_id in transcripts_dict:
        transcript = transcripts_dict[transcript_id]
    elif 'unknown_transcript' in transcript_id:
        continue
    else:
        raise Exception('transcript for start codon with transcript ID '+ transcript_id+ ' not found!!')

    transcript_parent = start_codon.attributes['Parent']
    parent_gene = genes_dict[transcript_parent]

    if transcript.strand == '+':
        transcript.set_thick_start(start_codon.start)
    elif transcript.strand == '-':
        transcript.set_thick_end(start_codon.end)

for stop_codon in stop_codons:
    transcript_id = stop_codon.attributes['ID']
    if transcript_id in transcripts_dict:
        transcript = transcripts_dict[transcript_id]
    elif 'unknown_transcript' in transcript_id:
        continue
    else:
        raise Exception('transcript for start codon with transcript ID '+ transcript_id+ ' not found!!')

    transcript_parent = stop_codon.attributes['Parent']
    parent_gene = genes_dict[transcript_parent]

    if transcript.strand == '+':
        transcript.set_thick_end(stop_codon.end)
    elif transcript.strand == '-':
        transcript.set_thick_start(stop_codon.start)


#print the genes
for gene in genes:
    print gene.bed_line()


#print the transcripts
for transcript in sorted(transcripts_dict.values(), key=lambda t: t.start):
    print transcript.bed_line()