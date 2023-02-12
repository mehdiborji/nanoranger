# extract reads containing differentially expressed exons
#
# exon list to be tested should be supplied in file in this format:
#
# CTLA4,exon2,chr2,203871378,203871487# 
# IL7R,exon6,chr5,35874449,35874542
# PTPRC,exon4,chr1,198696712,198696909

import argparse
import csv
import pandas as pd
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--bam', type=str, required=True, help='tagged bam file for processing')
parser.add_argument('--sample', type=str, required=True, help='sample name')
parser.add_argument('--exons', type=str, required=True, help='file with exons to be extracted')
parser.add_argument('--output', type=str, required=False, help='output file', default='none')

args = parser.parse_args()

r_seqs=[]

with open(args.exons, newline='') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for line in reader:
		print (line[0], ' ', line[1])
		gene = line[0]
		exon = line[1]
		chromosome = line[2]
		start = line[3]
		end = line[4]


		samfile = pysam.AlignmentFile(args.bam, 'rb')

		i=0
		for read in samfile.fetch(chromosome, int(start), int(end) ):
			if i > 1000000: # modify if required - maximum number of reads to be processed
				break
		
			r_seqs.append([args.sample, gene, exon, read.reference_start,read.reference_end,read.qlen,read.rlen,\
					read.get_tag('CB'),read.get_tag('UB'), read.get_overlap(int(start), int(end)), int(end)-int(start) ])
			i+=1
	
seqs=pd.DataFrame(r_seqs)
seqs.columns=['sample','gene', 'exon','ref_start','ref_end','query_length','read_length','bc','umi', 'overlap', 'length']
seqs[seqs.query_length==seqs.read_length]
if args.output != 'none':
	seqs[['sample', 'gene', 'exon', 'ref_start', 'ref_end', 'query_length', 'bc', 'umi', 'overlap', 'length']].to_csv(args.output,index=None)
else:
	print(seqs.to_string())
