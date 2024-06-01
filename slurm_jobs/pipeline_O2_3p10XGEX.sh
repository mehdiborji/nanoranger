#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=40G
#SBATCH -t 0:30:00
#SBATCH -p short
#SBATCH -o pipeline_3p10XGEX_job_%A.out
#SBATCH --account=chen_fec176

echo 'inputfq =' $1
echo 'outdir =' $2
echo 'sample =' $3
echo 'trans_ref =' $4
echo 'genome_ref =' $5

python ~/nanoranger/pipeline.py --c 16 --i $1 --o $2 --e $3 --s --m 3p10XGEX --t $4 --g $5
#~/refs/CAR.fa
#--g ~/refs/GRCh38.primary_assembly.genome.fa.gz 
#--g ~/refs/chr2.fasta.gz
