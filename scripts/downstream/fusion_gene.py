# read out reads mapping to gene from tagged bamfile alongside start and end of alignment

import pandas as pd
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--bam", type=str, required=True, help="tagged bam file for processing"
)
parser.add_argument("--output", type=str, required=True, help="output file")

args = parser.parse_args()

samfile = pysam.AlignmentFile(args.bam, "r")
alignments = []
for read in samfile.fetch():
    alignments.append(
        [
            read.reference_start,
            read.reference_end,
            read.reference_name,
            read.qlen,
            read.rlen,
            read.get_tag("CB"),
            read.get_tag("UB"),
        ]
    )

seqs = pd.DataFrame(alignments)

seqs.columns = [
    "ref_start",
    "ref_end",
    "ref_name",
    "query_length",
    "read_length",
    "BC",
    "UMI",
]

seqs.to_csv(args.output, index=None)
