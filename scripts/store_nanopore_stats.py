import pysam
from tqdm import tqdm
from datetime import datetime
import numpy as np
import csv
import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--indir", type=str)

args = parser.parse_args()
indir = args.indir


def process_fastq(fastq_file):

    t_ref = datetime(2024, 4, 23, 10, 20, 0)

    dir_name = os.path.dirname(fastq_file)
    base_name = os.path.splitext(os.path.basename(f))[0]

    # Construct the output file path
    output_file = os.path.join(dir_name, f"{base_name}.csv")

    with pysam.FastxFile(fastq_file) as fh, open(output_file, "w", newline="") as csvfile:

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
        for entry in tqdm(fh):
            i += 1

            seq = entry.sequence
            comms = entry.comment.split(" ")
            r_num = int(comms[2].split("=")[1])
            ch_num = int(comms[3].split("=")[1])
            t_read = datetime.strptime(comms[4].split("=")[1], "%Y-%m-%dT%H:%M:%SZ")
            t_read_ref = (t_read - t_ref).total_seconds() / 3600
            mean_q = np.mean(entry.get_quality_array())
            writer.writerow(
                [len(seq), mean_q, r_num, ch_num, t_read_ref]
            )  # Writing data to CSV
            if i > 400000:
                break

from concurrent.futures import ProcessPoolExecutor
import glob

fastq_files = glob.glob(f'{indir}/*.fastq.gz')  # Replace with your actual path

print(fastq_files)

with ProcessPoolExecutor() as executor:
    executor.map(process_fastq, fastq_files)