import subprocess
import argparse
import utils
import os
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('--cores', type=str)
parser.add_argument('--trns_ref', type=str)
parser.add_argument('--genome_ref', type=str)
parser.add_argument('--infile', type=str)
parser.add_argument('--outdir', type=str)
parser.add_argument('--expname', type=str)
parser.add_argument('--barcodes', type=str)
parser.add_argument('--split', default=False, action='store_true')

args = parser.parse_args()
cores = args.cores
trns_ref = args.trns_ref
genome_ref = args.genome_ref
infile = args.infile
outdir = args.outdir
sample = args.expname
barcodes = args.barcodes
split = args.split

pwd=os.path.dirname(os.path.abspath(__file__))

#print('pwd =',pwd)

if not os.path.exists(outdir): os.makedirs(outdir)

if trns_ref is None:
    #trns_ref = f'{pwd}/data/MT_trns.fa'
    trns_ref = f'{pwd}/data/panel_MT_trns.fa'

print('\n\n alignment to transcriptome reference and defusing/deconcatenation \n\n')

if split:
    
    subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
    fqs=[f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')]
    
    def align_trns(i):
        subprocess.call([ f'{pwd}/scripts/align_trns.sh','1', trns_ref,fqs[i],f'{outdir}/split',f'part_00{i+1}'])

    pool = Pool(int(cores))
    results = pool.map(align_trns, range(int(cores)))
    pool.close()
    
    #sams=[f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('sam')]

    args=[]
    for i in range(int(cores)): args.append((f'part_00{i+1}',f'{outdir}/split'))
    
    pool = Pool(int(cores))
    results = pool.starmap(utils.deconcatenate, args)
    pool.close()
    
    subprocess.call(f'cat {outdir}/split/*_decon*.gz > {outdir}/{sample}_deconcat.fastq.gz',shell=True)
    subprocess.call(f'cat {outdir}/split/*_BCUMI*.gz > {outdir}/{sample}_BCUMI.fasta.gz',shell=True)
    subprocess.call(f'rm -r {outdir}/split/',shell=True)

else:
    subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
    utils.deconcatenate(sample,outdir)
    subprocess.call(f'rm {outdir}/{sample}/',shell=True)
    
print('\n\n alignment of BC-UMI candidates to a reference of barcodes with STAR  \n\n')

if barcodes is None:
    barcodes=f'{pwd}/data/737K-august-2016.txt.gz'

utils.write_bc_10X5p(sample,outdir,barcodes)


subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])

subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
       f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])

print('\n\n alignment to genome and generation of BC-UMI-Transcript tagged BAM \n\n')

if genome_ref is None:
    genome_ref = f'{pwd}/data/MT_chr.fa'
    
subprocess.call([ f'{pwd}/scripts/align_genome.sh', cores, genome_ref, f'{outdir}/{sample}_deconcat.fastq.gz', outdir, sample])
"""
"""
utils.process_matching_10X5p(sample,outdir)


subprocess.call(f'samtools index -@{cores} {outdir}/{sample}_genome_tagged.bam',shell=True)
subprocess.call(f'rm -r {outdir}/{sample}_ref',shell=True)
subprocess.call(f'rm -r {outdir}/{sample}_matching_*',shell=True)
subprocess.call(f'rm -r {outdir}/{sample}_bcreads*',shell=True)
subprocess.call(f'rm -r {outdir}/{sample}_genome.*',shell=True)