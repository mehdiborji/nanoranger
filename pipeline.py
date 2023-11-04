import subprocess
import argparse
import utils
import os
from multiprocessing import Pool
import scanpy as sc

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
parser.add_argument('--xpecies', type=str)
#parser.add_argument('--human', default=False, action='store_true')

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
#human = args.human
xpecies = args.xpecies

pwd=os.path.dirname(os.path.abspath(__file__))

#print('pwd =',pwd)

if not os.path.exists(outdir): os.makedirs(outdir)

def align_trns(i):
    
    print('working on alignment of ', fqs[i])
    
    sam_file=f'{outdir}/split/part_{i+1}_trns.sam'
    if os.path.isfile(sam_file):
        print('alignment file,', sam_file,' exists, will not realign')
    else:
        print('alignment file,', sam_file,' does not exist, will align')
        subprocess.call([ f'{pwd}/scripts/align_trns.sh',cores, trns_ref,fqs[i],f'{outdir}/split',f'part_{i+1}'])
        
def align_bcs(i):
    
    sam_file=f'{outdir}/split/part_{i+1}_matching.sam'
    
    if os.path.isfile(sam_file):
        print('matching file,', sam_file,' exists, will not realign')
    else:
        print('matching file,', sam_file,' does not exist, will align')
        subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/split/part_{i+1}_BCUMI.fasta.gz', 
           f'{outdir}/{sample}_ref/', f'{outdir}/split/part_{i+1}_matching', cores, '-1'])
        
def split_fastq(infile,outdir,cores):
    
    inputfq_name=infile.split('/')[-1]

    inputfq_name=inputfq_name.split('.fastq')[0]

    splitted_file = f'{outdir}/split/{inputfq_name}.part_001.fastq.gz'

    if os.path.isfile(splitted_file):
        print('splitted input fastq exists, will not perform splitting assuming completion of splitting')
    else:
        print('splitted input fastq does not exist, will perform splitting')
        subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])

def split_fastq_unzipped(infile,outdir,cores):
    
    infq_fname = infile.split('/')[-1]

    infq_name = infq_fname.split('.fastq')[0]

    splitted_file = f'{outdir}/split/{infq_name}.part_001.fastq'

    if os.path.isfile(splitted_file):
        print('splitted input fastq exists, will not perform splitting assuming completion of splitting')
    else:
        print('splitted input fastq does not exist, will perform splitting')
        
        if infile.endswith('gz'):
            
            infile_unz = infile.replace('.gz','')
            if os.path.isfile(infile_unz):
                print(infile_unz,'unzipped file exists')
                subprocess.call(['seqkit', 'split2' ,infile_unz, '-p', cores, '-f', '-O', f'{outdir}/split'])
            else:
                
                if os.path.isfile(f'{outdir}/{infq_name}.fastq'):
                    subprocess.call(['seqkit', 'split2' ,f'{outdir}/{infq_name}.fastq', '-p', cores, '-f', '-O', f'{outdir}/split'])
                else:
                    print(infile_unz,'unzipped file and unzipped local does not exist')
                    print(f'pigz -cd {infile} > {outdir}/{infq_name}.fastq')
                    command=f'pigz -cd {infile} > {outdir}/{infq_name}.fastq'
                    subprocess.call(command,shell=True)
                    
                    subprocess.call(['seqkit', 'split2' ,f'{outdir}/{infq_name}.fastq', '-p', cores, '-f', '-O', f'{outdir}/split'])
                
        
                
        #command=f'rm {unR1} {unR2}'
        #subprocess.call(command,shell=True)

if mode == '5p10XGEX':
    if trns_ref is None:
        trns_ref = f'{pwd}/data/panel_MT_trns.fa'
    """"""
    print('\n\n alignment to transcriptome reference and defusing/deconcatenation \n\n')

    if split:
        
        subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
        fqs=sorted([f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')])

        for i in range(int(cores)): align_trns(i)
        
        args=[]
        for i in range(int(cores)): args.append((f'part_{i+1}',f'{outdir}/split'))

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
        trns_ref = f'{pwd}/data/TR_V_human.fa'
    if xpecies is None:
        xpecies = 'hsa'
    """"""
    print('\n\n alignment to transcriptome reference and defusing/deconcatenation \n\n')

    if split:
        
        inputfq_name=inputfq.split('/')[-1]

        inputfq_name=inputfq_name.split('.fastq')[0]

        split_fastq = f'{outdir}/split/{inputfq_name}.part_001.fastq.gz'

        if os.path.isfile(split_fastq):
            print(split_fastq,' exists')
        else:
            print(split_fastq,' does not exist, will perform splitting')
            subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
        
        fqs=sorted([f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')])
        for i in range(int(cores)): align_trns(i)

        args=[]
        
        for i in range(int(cores)): args.append((f'part_{i+1}',f'{outdir}/split'))
        
        decon_file= f'{outdir}/{sample}_deconcat.fastq.gz'
        
        if os.path.isfile(decon_file):
            print(decon_file,' exists')
        else:
            print(decon_file,' does not exist, will perform splitting')

            pool = Pool(int(cores))
            results = pool.starmap(utils.decon_5p10XTCR, args)
            pool.close()

            subprocess.call(f'cat {outdir}/split/*_decon*.gz > {outdir}/{sample}_deconcat.fastq.gz',shell=True)
            subprocess.call(f'cat {outdir}/split/*_BCUMI*.gz > {outdir}/{sample}_BCUMI.fasta.gz',shell=True)
            subprocess.call(f'cat {outdir}/split/*_eds*.csv > {outdir}/{sample}_eds.csv',shell=True)
        
            #subprocess.call(f'rm -r {outdir}/split/',shell=True)

    else:
        subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
        utils.decon_5p10XTCR(sample,outdir)
        #subprocess.call(f'rm {outdir}/*.sam',shell=True)

    subprocess.call([ f'{pwd}/scripts/mixcr.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_deconcat.fastq.gz', xpecies, cores ])
    
    #subprocess.call([ f'{pwd}/scripts/mixcr_asmbl.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_deconcat.fastq.gz',xpecies, cores ])
    
    if barcodes is None:
        barcodes=f'{pwd}/data/737K-august-2016.txt.gz'

    utils.write_bc_5p10X(sample,outdir,barcodes)

    subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])

    subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
           f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])
    
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
        for i in range(int(cores)): args.append((f'part_{i+1}',f'{outdir}/split'))

        pool = Pool(int(cores))
        results = pool.starmap(utils.decon_RTX, args)
        pool.close()

        subprocess.call(f'cat {outdir}/split/*_decon*.gz > {outdir}/{sample}_deconcat.fastq.gz',shell=True)
        
        #subprocess.call(f'rm -r {outdir}/split/',shell=True)

    else:
        subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
        utils.decon_RTX(sample,outdir)
        #subprocess.call(f'rm {outdir}/*.sam',shell=True)

    subprocess.call([ f'{pwd}/scripts/mixcr.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_deconcat.fastq.gz',xpecies, cores ])
    
    
if mode == '3pXCR_slideseq':
    
    if trns_ref is None:
        #trns_ref = f'{pwd}/data/IG_C_human.fa'
        trns_ref = f'{pwd}/data/XR_C_mouse.fa'
    
    print('\n\n alignment to C gene with minimap2 \n\n')
    if split:
        
        split_fastq = f'{outdir}/split/{sample}.part_001.fastq.gz'
        
        if os.path.isfile(split_fastq):
            print(split_fastq,' exists')
        else:
            subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
        
        fqs=sorted([f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')])
        for i in range(int(cores)): align_trns(i)
        
        args=[]
        for i in range(int(cores)): args.append((f'part_{i+1}',f'{outdir}/split'))

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
    

    subprocess.call([ f'{pwd}/scripts/mixcr.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_VDJ.fastq.gz', xpecies, cores ])
    
    #print('\n\n align BC-UMI to a reference of barcodes with STAR  \n\n')

    utils.write_bc_slideseq(sample,outdir,barcodes)

    subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])

    subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
           f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])

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
        for i in range(int(cores)): args.append((f'part_{i+1}',f'{outdir}/split'))

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

    subprocess.call([ f'{pwd}/scripts/mixcr.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_VDJ.fastq.gz',xpecies, cores ])
    
    

    print('\n\n align BC-UMI to a reference of barcodes with STAR  \n\n')

    utils.write_bc_slideseq(sample,outdir,barcodes)

    subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])

    subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
           f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])

    print('\n\n generate clone-barcode-UMI table \n\n')

    clones,cloneID=utils.clone_filt_slideseq(sample,outdir)

    utils.process_matching_slideseq_XCR(sample,outdir,cloneID)
    """"""
    
if mode == '3p10XTCR_nuc':
    if trns_ref is None:
        #trns_ref = f'{pwd}/data/MT_trns.fa'
        trns_ref = f'{pwd}/data/TR_V_human.fa'
    """"""
    print('\n\n alignment to transcriptome reference and defusing/deconcatenation \n\n')

    if split:
        
        split_fastq(infile,outdir,cores)
        
        fqs=sorted([f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('gz')])
        for i in range(int(cores)): align_trns(i)

        args=[]
        
        for i in range(int(cores)): args.append((f'part_{i+1}',f'{outdir}/split'))
        
        pool = Pool(int(cores))
        results = pool.starmap(utils.decon_3p10XTCR_nuc, args)
        pool.close()

        subprocess.call(f'cat {outdir}/split/*_VDJ*.gz > {outdir}/{sample}_VDJ.fastq.gz',shell=True)
        subprocess.call(f'cat {outdir}/split/*_BCUMI*.gz > {outdir}/{sample}_BCUMI.fasta.gz',shell=True)
        subprocess.call(f'cat {outdir}/split/*_eds*.csv > {outdir}/{sample}_eds.csv',shell=True)
        
        #subprocess.call(f'rm -r {outdir}/split/',shell=True)

    else:
        #subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
        utils.decon_3p10XTCR_nuc(sample,outdir)
        #subprocess.call(f'rm {outdir}/*.sam',shell=True)
    
    subprocess.call([ f'{pwd}/scripts/mixcr.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_VDJ.fastq.gz',xpecies, cores ])
    
    #subprocess.call([ f'{pwd}/scripts/mixcr_asmbl.sh', f'{outdir}/{sample}', f'{outdir}/{sample}_deconcat.fastq.gz',xpecies, cores ])
    
    if barcodes is None:
        barcodes=f'{pwd}/data/3M-february-2018.txt.gz'
    
    utils.write_bc_3p10XTCR_nuc(sample,outdir,barcodes)

    subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])

    subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
           f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])
    """"""
    print('\n\n generate clone-barcode-UMI table \n\n')
    
    utils.clone_filt_5p10X(sample,outdir)
    
    utils.process_matching_3p10XTCR_nuc(sample,outdir)
    """"""
    
if mode == '3p10XGEX_PacBio':
    if trns_ref is None:
        trns_ref = f'{pwd}/data/panel_MT_trns.fa'
    
    print('\n\n alignment to transcriptome reference and defusing/deconcatenation \n\n')

    if split:
        
        split_fastq_unzipped(infile,outdir,cores)
        #subprocess.call(['seqkit', 'split2' ,infile, '-p', cores, '-f', '-O', f'{outdir}/split'])
        fqs=sorted([f'{outdir}/split/'+f for f in os.listdir(f'{outdir}/split') if f.endswith('fastq')])
        
        print('\n ----alignment to transcriptome---- \n')
        for i in range(int(cores)): align_trns(i)
        
        args=[]
        for i in range(int(cores)): args.append((f'part_{i+1}',f'{outdir}/split'))
        
        print('\n ----deconcatenation---- \n')
        pool = Pool(int(cores))
        results = pool.starmap(utils.decon_3p10XGEX, args)
        pool.close()

        #subprocess.call(f'cat {outdir}/split/*_decon*.gz > {outdir}/{sample}_deconcat.fastq.gz',shell=True)
        #subprocess.call(f'cat {outdir}/split/*_BCUMI*.gz > {outdir}/{sample}_BCUMI.fasta.gz',shell=True)
        #subprocess.call(f'rm -r {outdir}/split/',shell=True)

    else:
        subprocess.call([ f'{pwd}/scripts/align_trns.sh', cores, trns_ref, infile, outdir, sample])
        utils.decon_3p10XGEX(sample,outdir)
        #subprocess.call(f'rm {outdir}/*.sam',shell=True)

    print('\n\n alignment of BC-UMI candidates to a reference of barcodes with STAR  \n\n')
    
    if barcodes is None:
        whitelist = f'{pwd}/data/3M-february-2018.txt.gz'

    utils.write_bc_3p10XGEX(sample, outdir, whitelist)
    
    subprocess.call([ f'{pwd}/scripts/barcode_ref.sh', f'{outdir}/{sample}_bcreads.fasta', f'{outdir}/{sample}_ref/'])
    
    for i in range(int(cores)): align_bcs(i)
    
    args=[]
    
    for i in range(int(cores)): args.append((f'part_{i+1}',f'{outdir}/split'))
    
    #if not os.path.exists(f'{outdir}/split/barcodes'): os.makedirs(f'{outdir}/split/barcodes')

    pool = Pool(int(cores))
    results = pool.starmap(utils.process_matching_3p10XGEX, args)
    pool.close()
    
    """
   

    #subprocess.call([ f'{pwd}/scripts/barcode_align.sh', f'{outdir}/{sample}_BCUMI.fasta.gz', 
    #       f'{outdir}/{sample}_ref/', f'{outdir}/{sample}_matching', cores, '-1'])
    
    print('\n\n alignment to genome and generation of BC-UMI-Transcript tagged BAM \n\n')

    if genome_ref is None:
        genome_ref = f'{pwd}/data/MT_chr.fa'

    subprocess.call([ f'{pwd}/scripts/align_genome.sh', cores, genome_ref, f'{outdir}/{sample}_deconcat.fastq.gz', outdir, sample])
    
    
    utils.process_matching_5p10X(sample,outdir)

    subprocess.call(f'samtools index -@{cores} {outdir}/{sample}_genome_tagged.bam',shell=True)
    subprocess.call(f'rm -r {outdir}/{sample}_ref',shell=True)
    subprocess.call(f'rm -r {outdir}/{sample}_matching_*',shell=True)
    subprocess.call(f'rm -r {outdir}/{sample}_bcreads*',shell=True)
    subprocess.call(f'rm -r {outdir}/{sample}_genome.*',shell=True)
    """