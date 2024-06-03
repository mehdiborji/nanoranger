#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=2G
#SBATCH -t 0:10:00
#SBATCH -p priority
#SBATCH -o poreranger_5p10XTCR_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'time =' $2
echo 'dev_basecall =' $3
echo 'total_reads =' $4
echo 'out_file =' $5

if [ "$3" = "dev_basecall" ] ; then
    python ~/nanoranger/scripts/store_nanopore_stats.py --indir $1 --time "$2" --dev_basecall --total_reads $4 --out_file $5

else
    python ~/nanoranger/scripts/store_nanopore_stats.py --indir $1 --time "$2" --total_reads $4 --out_file $5
fi

