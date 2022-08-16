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
parser.add_argument('--mode', type=str)
parser.add_argument('--human', default=False, action='store_true')

args = parser.parse_args()
cores = args.cores
trns_ref = args.trns_ref
genome_ref = args.genome_ref
infile = args.infile
outdir = args.outdir
sample = args.expname
barcodes = args.barcodes
split = args.split
mode = args.mode
human = args.human

pwd=os.path.dirname(os.path.abspath(__file__))

#print('pwd =',pwd)

if not os.path.exists(outdir): os.makedirs(outdir)

def align_trns(i):
    subprocess.call([ f'{pwd}/scripts/align_trns.sh',cores, trns_ref,fqs[i],f'{outdir}/split',f'part_00{i+1}'])

if mode == '5p10XGEX':
    if trns_ref is None:
        #trns_ref = f'{pwd}/data/MT_trns.fa'
        trns_ref = f'{pwd}/data/panel_MT_trns.fa'
    """"""
    print('\n\n alignment to transcriptome reference and defusing/deconcatenation \n\n')

    if split:
        
        subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
        fqs=[f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')]
        
        for i in range(int(cores)): align_trns(i)
        #pool = Pool(int(int(cores)/2))
        
        #results = pool.map(align_trns, range(int(int(cores)/2)))
        #pool.close()

        #sams=[f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('sam')]

        args=[]
        for i in range(int(cores)): args.append((f'part_00{i+1}',f'{outdir}/split'))

        pool = Pool(int(cores))
        results = pool.starmap(utils.decon_5p10XGEX, args)
        pool.close()

        subprocess.call(f'cat {outdir}/split/*_decon*.gz > {outdir}/{sample}_deconcat.fastq.gz',shell=True)
        subprocess.call(f'cat {outdir}/split/*_BCUMI*.gz > {outdir}/{sample}_BCUMI.fasta.gz',shell=True)
        #subprocess.call(f'rm -r {outdir}/split/',shell=True)

    else:
        subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
        utils.decon_5p10XGEX(sample,outdir)
        subprocess.call(f'rm {outdir}/*.sam',shell=True)

    print('\n\n alignment of BC-UMI candidates to a reference of barcodes with STAR  \n\n')


    if barcodes is None:
        barcodes=f'{pwd}/data/737K-august-2016.txt.gz'

    utils.write_bc_5p10X(sample,outdir,barcodes)


    subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])

    subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
           f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])
    
    print('\n\n alignment to genome and generation of BC-UMI-Transcript tagged BAM \n\n')

    if genome_ref is None:
        genome_ref = f'{pwd}/data/MT_chr.fa'

    subprocess.call([ f'{pwd}/scripts/align_genome.sh', cores, genome_ref, f'{outdir}/{sample}_deconcat.fastq.gz', outdir, sample])
    """
    """
    utils.process_matching_5p10X(sample,outdir)

    subprocess.call(f'samtools index -@{cores} {outdir}/{sample}_genome_tagged.bam',shell=True)
    subprocess.call(f'rm -r {outdir}/{sample}_ref',shell=True)
    subprocess.call(f'rm -r {outdir}/{sample}_matching_*',shell=True)
    subprocess.call(f'rm -r {outdir}/{sample}_bcreads*',shell=True)
    subprocess.call(f'rm -r {outdir}/{sample}_genome.*',shell=True)
    
    
if mode == '5p10XTCR':
    if trns_ref is None:
        #trns_ref = f'{pwd}/data/MT_trns.fa'
        trns_ref = f'{pwd}/data/TR_V_human.fa'
    """"""
    print('\n\n alignment to transcriptome reference and defusing/deconcatenation \n\n')

    if split:
        subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
        fqs=[f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')]

        pool = Pool(int(cores))
        results = pool.map(align_trns, range(int(cores)))
        pool.close()
        args=[]
        
        for i in range(int(cores)): args.append((f'part_00{i+1}',f'{outdir}/split'))

        pool = Pool(int(cores))
        results = pool.starmap(utils.decon_5p10XTCR, args)
        pool.close()

        subprocess.call(f'cat {outdir}/split/*_decon*.gz > {outdir}/{sample}_deconcat.fastq.gz',shell=True)
        subprocess.call(f'cat {outdir}/split/*_BCUMI*.gz > {outdir}/{sample}_BCUMI.fasta.gz',shell=True)
        subprocess.call(f'cat {outdir}/split/*_eds*.csv > {outdir}/{sample}_eds.csv',shell=True)
        
        subprocess.call(f'rm -r {outdir}/split/',shell=True)

    else:
        subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
        utils.decon_5p10XTCR(sample,outdir)
        #subprocess.call(f'rm {outdir}/*.sam',shell=True)
    
    if not human: species='mmu'; print(species)
    else: species='hsa'; print(species)

    subprocess.call([ f'{pwd}/scripts/mixcr.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_deconcat.fastq.gz',species, cores ])
    
    #subprocess.call([ f'{pwd}/scripts/mixcr_asmbl.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_deconcat.fastq.gz',species, cores ])
    
    if barcodes is None:
        barcodes=f'{pwd}/data/737K-august-2016.txt.gz'

    utils.write_bc_5p10X(sample,outdir,barcodes)


    subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])

    subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
           f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])
    """"""
    print('\n\n generate clone-barcode-UMI table \n\n')
    
    utils.clone_filt_5p10X(sample,outdir)
    
    utils.process_matching_5p10XTCR(sample,outdir)
    """"""
    
    
if mode == 'RTX':
    if trns_ref is None:
        #trns_ref = f'{pwd}/data/MT_trns.fa'
        trns_ref = f'{pwd}/data/TR_V_human.fa'
    """"""
    print('\n\n alignment to transcriptome reference and defusing/deconcatenation \n\n')

    if split:
        
        subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
        fqs=[f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')]

        pool = Pool(int(cores))
        results = pool.map(align_trns, range(int(cores)))
        pool.close()
        
        args=[]
        for i in range(int(cores)): args.append((f'part_00{i+1}',f'{outdir}/split'))

        pool = Pool(int(cores))
        results = pool.starmap(utils.decon_RTX, args)
        pool.close()

        subprocess.call(f'cat {outdir}/split/*_decon*.gz > {outdir}/{sample}_deconcat.fastq.gz',shell=True)
        
        #subprocess.call(f'rm -r {outdir}/split/',shell=True)

    else:
        subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
        utils.decon_RTX(sample,outdir)
        #subprocess.call(f'rm {outdir}/*.sam',shell=True)
    
    if not human: species='mmu'; print(species)
    else: species='hsa'; print(species)

    subprocess.call([ f'{pwd}/scripts/mixcr.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_deconcat.fastq.gz',species, cores ])
    
    
if mode == '3pXCR_slideseq':
    
    if trns_ref is None:
        #trns_ref = f'{pwd}/data/IG_C_human.fa'
        trns_ref = f'{pwd}/data/XR_C_mouse.fa'
    """
    print('\n\n alignment to C gene with minimap2 \n\n')
    if split:
        
        
        subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
        fqs=[f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')]

        pool = Pool(int(cores))
        results = pool.map(align_trns, range(int(cores)))
        pool.close()
        
        args=[]
        for i in range(int(cores)): args.append((f'part_00{i+1}',f'{outdir}/split'))

        pool = Pool(int(cores))
        results = pool.starmap(utils.decon_3pXCR_slideseq, args)
        pool.close()
        
        
        subprocess.call(f'cat {outdir}/split/*_VDJ.fastq.gz > {outdir}/{sample}_VDJ.fastq.gz',shell=True)
        subprocess.call(f'cat {outdir}/split/*_BCUMI.fasta.gz > {outdir}/{sample}_BCUMI.fasta.gz',shell=True)
        subprocess.call(f'cat {outdir}/split/*_polyA.fasta.gz > {outdir}/{sample}_polyA.fasta.gz',shell=True)
        
        #subprocess.call(f'rm -r {outdir}/split/',shell=True)

    else:
        subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
        utils.decon_3pXCR_slideseq(sample,outdir)
        #subprocess.call(f'rm {outdir}/*.sam',shell=True)
        
    
    print('\n\n align VDJ with MiXCR and extract clones \n\n')

    if not human: species='mmu'; print(species)
    else: species='hsa'; print(species)

    subprocess.call([ f'{pwd}/scripts/mixcr.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_VDJ.fastq.gz',species, cores ])

    print('\n\n align BC-UMI to a reference of barcodes with STAR  \n\n')

    utils.write_bc_slideseq(sample,outdir,barcodes)

    subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])

    subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
           f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])
           
    """

    print('\n\n generate clone-barcode-UMI table \n\n')

    clones,cloneID=utils.clone_filt_slideseq(sample,outdir)

    utils.process_matching_slideseq_XCR(sample,outdir,cloneID)
    """ """
    
if mode == '3p10XTCR':
    
    if trns_ref is None:
        #trns_ref = f'{pwd}/data/IG_C_human.fa'
        trns_ref = f'{pwd}/data/TRab_C_mouse.fa'
    """"""
    print('\n\n alignment to C gene with minimap2 \n\n')
    if split:
        
        subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
        fqs=[f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')]

        pool = Pool(int(cores))
        results = pool.map(align_trns, range(int(cores)))
        pool.close()
        
        args=[]
        for i in range(int(cores)): args.append((f'part_00{i+1}',f'{outdir}/split'))

        pool = Pool(int(cores))
        results = pool.starmap(utils.landmark_extract_slideseq_XCR, args)
        pool.close()

        subprocess.call(f'cat {outdir}/split/*_decon*.gz > {outdir}/{sample}_deconcat.fastq.gz',shell=True)
        
        #subprocess.call(f'rm -r {outdir}/split/',shell=True)

    else:
        #subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
        utils.decon_3p10XTCR(sample,outdir)
        #subprocess.call(f'rm {outdir}/*.sam',shell=True)
        
    
    print('\n\n align VDJ with MiXCR and extract clones \n\n')

    if not human: species='mmu'; print(species)
    else: species='hsa'; print(species)

    subprocess.call([ f'{pwd}/scripts/mixcr.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_VDJ.fastq.gz',species, cores ])
    
    """

    print('\n\n align BC-UMI to a reference of barcodes with STAR  \n\n')

    utils.write_bc_slideseq(sample,outdir,barcodes)

    subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])

    subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
           f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])

    print('\n\n generate clone-barcode-UMI table \n\n')

    clones,cloneID=utils.clone_filt_slideseq(sample,outdir)

    utils.process_matching_slideseq_XCR(sample,outdir,cloneID)
    """