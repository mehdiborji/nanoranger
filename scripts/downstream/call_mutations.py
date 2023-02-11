import matplotlib
matplotlib.use('Agg')

import argparse
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import sys
import time
from pyliftover import LiftOver
lo = LiftOver('hg19', 'hg38')
import subprocess
import csv

def sort_cnt(arr):
    npcnt=np.array(np.unique(arr,return_counts=True)).T
    dfcnt=pd.DataFrame(npcnt)
    dfcnt[1]=dfcnt[1].astype('int')
    dfcnt=dfcnt.sort_values(by=1,ascending=False)
    return dfcnt

parser = argparse.ArgumentParser()
parser.add_argument('--outdir', type=str, required=True, help='directory where output will be stored')
parser.add_argument('--mutations', type=str, required=True, help='file containing mutations')
parser.add_argument('--bam', type=str, required=True, help='tagged bam file for processing')
parser.add_argument('--sample', type=str, required=True, help='sample name')
args = parser.parse_args()

outdir = args.outdir
sam_tag = args.bam
sample = args.sample
sam_tag = args.bam

muts=pd.read_csv(args.mutations,index_col=0)

muts=muts.loc[sample].copy()
try:
    print(muts.shape[1])
except:
    print('one mut only')
    muts=pd.DataFrame(muts).T
    
samfile = pysam.AlignmentFile(sam_tag, 'rb',threads=4)
for gene in muts.gene.to_list():
    with open(f'./{outdir}/{sample}_pileup_{gene}.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['bc','umi','base','Q','indel'])
        chrm=muts[muts.gene==gene].chr.to_list()[0]
        start=muts[muts.gene==gene].pos.to_list()[0]-1
        end=start+1
        print(gene,chrm,start,end)
        for pileupcolumn in samfile.pileup(chrm,start,end,min_base_quality=0,max_depth=3000000):
            if pileupcolumn.pos==start:
                print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        
                        line = [pileupread.alignment.get_tag('CB'),pileupread.alignment.get_tag('UB'),
                                pileupread.alignment.query_sequence[pileupread.query_position],
                                pileupread.alignment.query_qualities[pileupread.query_position],pileupread.indel]
                        writer.writerow(line)
                break
    subprocess.call([ 'pigz', '-f', f'./{outdir}/{sample}_pileup_{gene}.csv'])

sys.exit(0)

plt.rcParams['figure.figsize'] = (5, 5)

for gene in muts.gene.to_list():
    piles=pd.read_csv(f'./{sample}/{sample}_pileup_{gene}.csv.gz')
    
    cnt=sort_cnt(piles['bc'])

    cnt.columns=['bc','cnt']
    print(cnt[cnt.cnt>10].shape[0])
    plt.figure(figsize=(5, 5))
    plt.plot(np.log10(np.arange(1,len(cnt.cnt)+1)),np.log10(cnt.cnt));
    plt.ylabel('log10 read counts');
    plt.xlabel('log10 cell rank');
    plt.title(gene);
    plt.savefig(f'./{sample}/{sample}_knee_{gene}.pdf',bbox_inches='tight');
