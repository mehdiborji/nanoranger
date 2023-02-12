import pandas as pd
import pysam
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bam', type=str, required=True, help='tagged bam file for processing')
parser.add_argument('--gene', type=str, required=True, help='gene name')
parser.add_argument('--output', type=str, required=True, help='output file')
parser.add_argument('--begin', type=str, required=True, help='start of gene')
parser.add_argument('--end', type=str, required=True, help='end of gene')

args = parser.parse_args()

samfile = pysam.AlignmentFile(args.bam, 'rb')
r_seqs=[]
for read in samfile.fetch(args.gene, int(args.begin), int(args.end)):
	r_seqs.append([read.reference_start,read.reference_end,read.qlen,read.rlen,read.get_tag('CB'),read.get_tag('UB')])
	
seqs=pd.DataFrame(r_seqs)
seqs.columns=['ref_start','ref_end','query_length','read_length','bc','umi']
seqs[seqs.query_length==seqs.read_length]

seqs['gene'] = args.gene

seqs[['gene','ref_start','ref_end','query_length','bc','umi']].to_csv(args.output,index=None)
