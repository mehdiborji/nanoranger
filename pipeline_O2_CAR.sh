#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 1:00:00
#SBATCH -p short
#SBATCH -o poreranger_job_%A.out

module load gcc/9.2.0
module load python/3.8.12

echo 'inputfq =' $1
echo 'outdir =' $2
echo 'sample =' $3
echo 'cores =' $4

python3.8 ./pipeline.py --s \
--i $1 --o $2 --e $3 --c $4 \
--t data/CAR.fa \
--g data/CAR.fa

### for CAR transcripts amplified with CD247 primer use CAR_CD247.fa as reference
#python3.8 ./pipeline.py --s \
#--i $1 --o $2 --e $3 --c $4 \
#--t data/CAR_CD247.fa \
#--g data/CAR_CD247.fa
