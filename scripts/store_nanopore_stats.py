import pysam
from tqdm import tqdm
from datetime import datetime
import numpy as np
import csv
import os
import argparse
from concurrent.futures import ProcessPoolExecutor
import glob
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument("--indir", type=str)
parser.add_argument("--dev_basecall", default=False, action="store_true")
parser.add_argument('--time', type=lambda s: datetime.strptime(s, '%Y-%m-%d %H:%M:%S'))
parser.add_argument("--total_reads", type=int, default=1e3)
parser.add_argument("--out_file", type=str)

args = parser.parse_args()
indir = args.indir
dev_basecall = args.dev_basecall
input_datetime = args.time
total_reads = args.total_reads
out_file = args.out_file

print('dev_basecall status',dev_basecall)

# cd /n/scratch/users/m/meb521/20240520_2209_MC-113060_FAU97225/other_reports
# sbatch ~/nanoranger/slurm_jobs/pipeline_O2_collect_stats.sh fastq_pass "2024-05-20 22:09:00" dev_basecall 50000 all_stats_fastq_pass.csv

def process_fastq(fastq_file):
    t_ref = input_datetime

    dir_name = os.path.dirname(fastq_file)
    base_name = os.path.basename(fastq_file)
    csv_name = base_name.split(".fastq")[0]

    output_file = os.path.join(dir_name, f"{csv_name}.csv")

    #print(output_file)

    with pysam.FastxFile(fastq_file) as fh, open(
        output_file, "w", newline=""
    ) as csvfile:
        i = 0
        writer = csv.writer(csvfile)
        """
        writer.writerow(
            [
                "Length",
                "Mean Quality",
                "Read Number",
                "Channel Number",
                "Time Read Reference",
            ]
        )  # Writing headers
        """
        for entry in fh:
            
            i += 1

            seq = entry.sequence
            comms = entry.comment.split(" ")
            
            if dev_basecall:
                
                #print(comms[1:4])
                r_num = int(comms[1].split("=")[1])
                ch_num = int(comms[2].split("=")[1])
                t_read = datetime.strptime(comms[3].split("=")[1].split(".")[0], "%Y-%m-%dT%H:%M:%S")
            else:
                r_num = int(comms[2].split("=")[1])
                ch_num = int(comms[3].split("=")[1])
                t_read = datetime.strptime(comms[4].split("=")[1], "%Y-%m-%dT%H:%M:%SZ")
                
            t_read_ref = (t_read - t_ref).total_seconds()
            mean_q = np.mean(entry.get_quality_array())
            writer.writerow(
                [len(seq), mean_q, r_num, ch_num, t_read_ref]
            )  # Writing data to CSV
            if i > total_reads:
                break


fastq_files = glob.glob(f"{indir}/*.fastq.gz")

#for f in fastq_files:
#    print(f)
#process_fastq(fastq_files[0])

with ProcessPoolExecutor(max_workers=20) as executor:
    executor.map(process_fastq, fastq_files)

subprocess.call(f"cat {indir}/*.csv > {out_file}",shell=True)
subprocess.call(f"cat {out_file} | wc -l", shell=True)