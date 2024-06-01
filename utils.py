import pysam
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import edlib
import subprocess
import json

print("nanoranger packages loaded")

linker = "TCTTCAGCGTTCCCGAGA"
ad_seq = [("N", "A"), ("N", "T"), ("N", "G"), ("N", "C")]


def rev(seq):
    return str(Seq(seq).reverse_complement())


def seq_counter(seq_dict, seq_instance):
    if seq_dict.get(seq_instance) is None:
        seq_dict[seq_instance] = 1
    else:
        seq_dict[seq_instance] += 1


def quad_dict_store(quad_dict, quad_key, quad_items):
    if quad_dict.get(quad_key) is None:
        quad_dict[quad_key] = [quad_items]
    else:
        quad_dict[quad_key].extend([quad_items])


def sort_cnt(arr):
    npcnt = np.array(np.unique(arr, return_counts=True)).T
    dfcnt = pd.DataFrame(npcnt)
    dfcnt[1] = dfcnt[1].astype("int")
    dfcnt = dfcnt.sort_values(by=1, ascending=False)
    return dfcnt


def decon_RTX(sample, outdir):
    rclipV = 100

    file = f"{outdir}/{sample}_trns.sam"

    f1 = open(f"{outdir}/{sample}_deconcat.fastq", "w")

    samfile = pysam.AlignmentFile(f"{file}", "r")

    for read in samfile.fetch():
        read.get_tag("AS")
        rlen = read.rlen
        seq = read.seq
        flag = read.flag
        qstrt = read.query_alignment_start
        qend = read.query_alignment_end
        # rname=read.qname;trans=read.reference_name.split('|')[4] #for mouse (unedited)
        trans = read.reference_name.split("|")[0]  # for human (edited)

        if flag == 16 or flag == 2064:
            qstrt_mod = rlen - qend
            qend_mod = rlen - qstrt
        else:
            qstrt_mod = qstrt
            qend_mod = qend

        if rlen - qend > rclipV:
            sub_end = qend + rclipV
        else:
            sub_end = rlen

        sub_seq = seq[qstrt:sub_end]
        sub_qual = read.qual[qstrt:sub_end]

        clip = sub_end - qend
        # subs.append(clip)
        # newnamef=f'>{read.qname}_{puck}_{sub_strt}_{sub_end}_{flag}_{trans}'
        # dd=end_qu
        newnamef = f"{read.qname}_{qstrt_mod}_{qend_mod}_{flag}_{trans}"
        if len(sub_seq) > 100 and clip > 40:
            f1.write(f"@{newnamef}\n")
            f1.write(f"{sub_seq}\n")
            f1.write("+\n")
            f1.write(f"{sub_qual}\n")
        # reads.append([rlen,qlen,qstrt,qend,qstrt_mod,qend_mod,flag,ref_len])
        # trns.append([rname,transs])
        # tot+=1
        # if tot>10000000000:break
    samfile.close()
    f1.close()
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_deconcat.fastq"])


def decon_5p10XGEX(sample, outdir):
    tot = 0
    reads = []
    trns = []
    eds = []
    lclip = 80
    rclip = 20
    r_hang = 0
    const = "CGCTCTTCCGATCT" + 26 * "N" + "TTTCTTATATG"  # motif for 5' 10x

    file = f"{outdir}/{sample}_trns.sam"

    f1 = open(f"{outdir}/{sample}_deconcat.fastq", "w")
    f2 = open(f"{outdir}/{sample}_BCUMI.fasta", "w")

    samfile = pysam.AlignmentFile(f"{file}", "r")

    for read in samfile.fetch():
        read.get_tag("AS")
        qlen = read.qlen
        rlen = read.rlen
        seq = read.seq
        qstrt = read.query_alignment_start
        qend = read.query_alignment_end
        rname = read.qname
        trans = read.reference_name
        flag = read.flag
        ref_len = read.reference_length
        ref_s = read.reference_start
        ref_e = read.reference_end
        ref_e - ref_s

        if qstrt > lclip:
            beg_qu = seq[qstrt - lclip : qstrt + rclip]
        else:
            beg_qu = seq[: qstrt + rclip]

        ed = edlib.align(const, beg_qu, "HW", "locations", 6, ad_seq)

        if ed["editDistance"] > -1 and ed["editDistance"] < 7:
            start = ed["locations"][-1][0]
            end = ed["locations"][-1][1]
            bcumi = beg_qu[start:end]
            if qstrt > lclip:
                start = lclip - start
            else:
                start = qstrt - start
            eds.append([start, len(bcumi), ed["editDistance"]])

            sub_end = qend + r_hang
            sub_strt = qstrt  # -start
            sub_seq = seq[sub_strt:sub_end]
            sub_qual = read.qual[sub_strt:sub_end]

            if flag == 16 or flag == 2064:
                qstrt_mod = rlen - qend
                qend_mod = rlen - qstrt
            else:
                qstrt_mod = qstrt
                qend_mod = qend

            reads.append([rlen, qlen, ref_len, qstrt_mod, qend_mod, flag])
            trns.append([rname, trans])

            newnamef = f"{read.qname}_{qstrt_mod}_{qend_mod}_{flag}_{trans}"
            f1.write(f"@{newnamef}\n")
            f1.write(f"{sub_seq}\n")
            f1.write("+\n")
            f1.write(f"{sub_qual}\n")

            f2.write(f">{newnamef}\n")
            f2.write(f"{bcumi}\n")
        tot += 1
        if tot % 10000 == 0:
            print(tot, f" records from {sample} processed")

    f1.close()
    f2.close()
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_deconcat.fastq"])
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_BCUMI.fasta"])


def decon_5p10XTCR(sample, outdir):
    tot = 0
    eds = []
    lclip = 200  # lclip large to account for missing 5' UTR in reference transcripts,
    rclip = 20  # rclip to account for accidental alignment of BC-UMI-TSO to transcripts (rare)

    # lclipV = 150;
    # must generate custom transcripts and report these unannotated regions

    # rclipV = 100;

    # lclipV: clip into 5' direction of of V gene, this can be technically 0

    lclipV = 40

    # rclipV: clip into direction of C gene, this will be covering CDR3 junction
    # and J gene region and additionally parts of C gene, optimize to clip on adapaterss

    rclipV = 80

    const = "CGCTCTTCCGATCT" + 26 * "N" + "TTTCTTATATG"

    file = f"{outdir}/{sample}_trns.sam"

    f1 = open(f"{outdir}/{sample}_deconcat.fastq", "w")
    f2 = open(f"{outdir}/{sample}_BCUMI.fasta", "w")

    samfile = pysam.AlignmentFile(f"{file}", "r")

    for read in samfile.fetch():
        read.get_tag("AS")
        qlen = read.qlen
        rlen = read.rlen
        seq = read.seq
        qstrt = read.query_alignment_start
        qend = read.query_alignment_end
        trans = read.reference_name
        flag = read.flag
        ref_s = read.reference_start
        ref_e = read.reference_end
        ref_e - ref_s

        if qstrt > lclip:
            beg_qu = seq[qstrt - lclip : qstrt + rclip]
        else:
            beg_qu = seq[: qstrt + rclip]

        """
        # select only lclipV bases from 3' end of V to the 5' direction
        if qlen>lclipV:
            sub_strt = qend-lclipV
        else:
            sub_strt = qstrt
            
        # select only rclipV bases from 3' end of V to the 3' direction

        if rlen-qend>rclipV:
            sub_end = qend+rclipV
        else:
            sub_end = rlen
        sub_seq=seq[sub_strt:sub_end]
        sub_qual=read.qual[sub_strt:sub_end]
        """
        if qlen > lclipV:
            sub_strt = qend - lclipV
        else:
            sub_strt = qstrt

        if rlen - qend > rclipV:
            sub_end = qend + rclipV
        else:
            sub_end = rlen

        sub_seq = read.seq[sub_strt:sub_end]
        sub_qual = read.qual[sub_strt:sub_end]

        if flag == 16 or flag == 2064:
            qstrt_mod = rlen - qend
            qend_mod = rlen - qstrt
        else:
            qstrt_mod = qstrt
            qend_mod = qend

        newnamef = f"{read.qname}_{qstrt_mod}_{qend_mod}_{flag}_{trans}"
        f1.write(f"@{newnamef}\n")
        f1.write(f"{sub_seq}\n")
        f1.write("+\n")
        f1.write(f"{sub_qual}\n")

        ed = edlib.align(const, beg_qu, "HW", "locations", 6, ad_seq)

        if ed["editDistance"] > -1 and ed["editDistance"] < 7:
            start = ed["locations"][-1][0]
            end = ed["locations"][-1][1]
            bcumi = beg_qu[start:end]
            if qstrt > lclip:
                start = lclip - start
            else:
                start = qstrt - start
            eds.append([start, trans, len(bcumi), ed["editDistance"]])

            # sub_end = qend;

            # reads.append([rlen,qlen,ref_len,qstrt_mod,qend_mod,flag])
            # trns.append([rname,trans])

            f2.write(f">{newnamef}\n")
            f2.write(f"{bcumi}\n")
        tot += 1
        if tot % 10000 == 0:
            print(tot, f" records from {sample} processed")

    f1.close()
    f2.close()
    eds = pd.DataFrame(np.array(eds))
    eds.to_csv(f"{outdir}/{sample}_eds.csv")
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_deconcat.fastq"])
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_BCUMI.fasta"])


def decon_3p10XTCR(sample, outdir):
    tot = 0
    eds = []
    file = f"{outdir}/{sample}_trns.sam"

    samfile = pysam.AlignmentFile(f"{file}", "r", threads=8)

    const = 6 * "A" + 28 * "N" + "AGATCGGAAGAGCGTCGTGT"

    lclip = 350
    r_search = 150
    f1 = open(f"{outdir}/{sample}_VDJ.fastq", "w")
    f2 = open(f"{outdir}/{sample}_BCUMI.fasta", "w")
    print(f"{outdir}/{sample}_VDJ.fastq")
    for read in samfile.fetch():
        qlen = read.qlen
        rlen = read.rlen
        seq = read.seq
        flag = read.flag
        qstrt = read.query_alignment_start
        qend = read.query_alignment_end
        rname = read.qname[-10:]
        trans = read.reference_name.split("-")[0]
        ref_s = read.reference_start
        ref_e = read.reference_end
        ref_e - ref_s

        rclip = 100

        if rlen - qend > r_search:
            end_qu = seq[qend : qend + r_search]
        else:
            end_qu = seq[qend:]

        sub_end = qstrt + rclip  # for including part of C gene
        if qstrt > lclip:
            sub_strt = qstrt - lclip
        else:
            sub_strt = 0

        sub_seq = seq[sub_strt:sub_end]
        sub_qual = read.qual[sub_strt:sub_end]

        ed = edlib.align(const, end_qu, "HW", "locations", 5, ad_seq)
        dist = ed["editDistance"]
        eds.append(dist)
        newname = f"{rname}_q{qlen}_d{dist}_s{sub_strt}_e{sub_end}_f{flag}_{trans}"

        if dist > -1 and dist < 6 and len(sub_seq) > 100 and qlen > 100:
            f1.write(f"@{newname}\n")
            f1.write(f"{sub_seq}\n")
            f1.write("+\n")
            f1.write(f"{sub_qual}\n")
            bcumi = rev(end_qu[ed["locations"][0][0] : ed["locations"][0][1]])[14:]
            f2.write(f">{newname}\n")
            f2.write(f"{bcumi}\n")
        tot += 1
        if tot % 50000 == 0:
            print(tot, " records processed")
    samfile.close()
    f1.close()
    f2.close()

    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_VDJ.fastq"])
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_BCUMI.fasta"])

    sort_cnt(eds).to_csv(f"{outdir}/{sample}_eds.csv")


def decon_3pXCR_slideseq(sample, outdir):
    tot = 0
    c_hangs = []
    c_eds = []
    polyAs = []
    newnames = []
    file = f"{outdir}/{sample}_trns.sam"

    subs = []
    samfile = pysam.AlignmentFile(f"{file}", "r", threads=8)
    const = rev(linker)

    rclip = 80  # selection length into 3' direction of 5' softclip of C gene (for including part of C gene)
    lclip = 200  # selection length into 5' direction of 5' softclip of C gene  (for including VDJ sequence) potential optimization is to clip on PCR adapter sequence

    rclipA = 22
    lclipA = 16  # nul

    r_search = 200  # search length into 3' direction of 3' softclip of C gene (to search for BC-UMI)

    f1 = open(f"{outdir}/{sample}_VDJ.fastq", "w")
    print(f"{outdir}/{sample}_VDJ.fastq")
    for read in samfile.fetch():
        accept = False
        read.get_tag("AS")
        rlen = read.rlen
        seq = read.seq
        flag = read.flag
        qstrt = read.query_alignment_start
        qend = read.query_alignment_end
        trans = read.reference_name

        ref_s = read.reference_start
        ref_e = read.reference_end
        span = ref_e - ref_s

        if span > 400:
            accept = True

        if rlen - qend > r_search:
            end_qu = seq[qend : qend + r_search]
        else:
            end_qu = seq[qend:]

        sub_end = qstrt + rclip
        if qstrt > lclip:
            sub_strt = qstrt - lclip
        else:
            sub_strt = 0
        # sub_seq=seq[sub_strt:sub_end]
        # sub_qual=read.qual[sub_strt:sub_end]
        sub_seq = seq[sub_strt:sub_end]
        sub_qual = read.qual[sub_strt:sub_end]

        subs.append(len(sub_seq))
        newnamef = f">{read.qname}_{sample}_{sub_strt}_{sub_end}_{flag}_{trans}"
        dd = end_qu

        if len(sub_seq) > 100 and accept:
            f1.write(f"@{read.qname}_{sample}_{sub_strt}_{sub_end}_{flag}_{trans}\n")
            f1.write(f"{sub_seq}\n")
            f1.write("+\n")
            f1.write(f"{sub_qual}\n")

            for i in range(int(len(dd) / 20)):
                w = dd[20 * i : 20 * i + 40]
                ed = edlib.align(const, w, "HW", "locations", 2)

                if ed["editDistance"] > -1 and ed["editDistance"] < 4:
                    # print(ed)
                    start = ed["locations"][0][0] + 20 * i
                    end = ed["locations"][0][1] + 20 * i
                    if start - rclipA < 0:
                        upstart = 0
                    else:
                        upstart = start - rclipA
                    upend = end + lclipA
                    select = rev(dd[upstart:upend])
                    polyA = dd[: upstart + 5]
                    c_hangs.append(select)
                    polyAs.append(polyA)
                    c_eds.append(ed["editDistance"])
                    newnames.append(newnamef)
                    break
        tot += 1
        if tot % 50000 == 0:
            print(tot, " records processed")
    samfile.close()
    f1.close()

    pd.DataFrame([newnames, c_eds]).T.to_csv(
        f"{outdir}/{sample}_eds_names.csv", index=None
    )

    f2 = open(f"{outdir}/{sample}_BCUMI.fasta", "w")
    f3 = open(f"{outdir}/{sample}_polyA.fasta", "w")

    for i in range(len(newnames)):
        accept = False
        if len(c_hangs[i]) > 45 and len(polyAs[i]) < 70:
            accept = True
        if len(c_hangs[i]) > 45 and len(polyAs[i]) > 70:
            if polyAs[i].count("A") / len(polyAs[i]) > 0.5:
                accept = True
        if accept:
            f2.write(f"{newnames[i]}\n")
            f2.write(f"{c_hangs[i]}\n")
            f3.write(f"{newnames[i]}\n")
            f3.write(f"{polyAs[i]}\n")
    f2.close()
    f3.close()

    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_VDJ.fastq"])
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_BCUMI.fasta"])
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_polyA.fasta"])
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_eds_names.csv"])


def clone_filt_slideseq(sample, outdir):
    # for slideseq

    clones = pd.read_table(f"{outdir}/{sample}_clones.txt.gz")
    # clones=clones[clones.cloneCount>1].copy()
    # clones=clones[clones.topChains==clones.chains].copy()
    clones = clones[
        [
            "chains",
            "aaSeqImputedCDR3",
            "cloneCount",
            "cloneId",
            "nSeqImputedCDR3",
            "allVHitsWithScore",
            "allDHitsWithScore",
            "allJHitsWithScore",
        ]
    ]
    # clones=clones[clones.chains.isin(('TRA','TRB'))].copy()
    cloneID = pd.read_table(f"{outdir}/{sample}_cloneID.txt.gz")
    # cloneID=cloneID[(cloneID.chains.isin(('TRA','TRB')))&(cloneID.cloneId.isin(clones.cloneId))].copy()

    cloneID = cloneID[cloneID.cloneId.isin(clones.cloneId)].copy()

    """

    ncdr=sort_cnt(clones.aaSeqImputedCDR3)
    repncdr=ncdr[ncdr[1]>1][0].to_list()
    reclone=clones[clones.aaSeqImputedCDR3.isin(repncdr)].sort_values(by=['aaSeqImputedCDR3','cloneCount'],ascending=False)

    maps={}
    for rep in repncdr:
        dd=reclone[reclone.aaSeqImputedCDR3==rep].index.tolist()
        for idx in dd[1:]:
            maps[idx]=dd[0]

    cloneID.cloneId=cloneID.cloneId.apply(lambda x: maps[x] if x in maps.keys() else x).copy()
    clones.drop_duplicates(subset='aaSeqImputedCDR3',keep='first',inplace=True)
    """
    cloneID.set_index("descrsR1", inplace=True)
    cloneID = cloneID[["chains", "cloneId"]]
    cloneID.to_csv(f"{outdir}/{sample}_cloneID_filtered.csv")
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_cloneID_filtered.csv"])
    print("clone filtering finished")
    return (clones, cloneID)


def clone_filt_5p10X(sample, outdir):
    clones = pd.read_table(f"{outdir}/{sample}_clones.txt.gz")
    clones = clones[clones.cloneCount > 1].copy()
    # clones=clones[clones.topChains==clones.chains].copy()
    clones = clones[
        [
            "chains",
            "aaSeqImputedCDR3",
            "cloneCount",
            "cloneId",
            "nSeqImputedCDR3",
            "allVHitsWithScore",
            "allDHitsWithScore",
            "allJHitsWithScore",
        ]
    ]
    # clones=clones[clones.chains.isin(('TRA','TRB'))].copy()
    cloneID = pd.read_table(f"{outdir}/{sample}_cloneID.txt.gz")
    # cloneID=cloneID[(cloneID.chains.isin(('TRA','TRB')))&(cloneID.cloneId.isin(clones.cloneId))].copy()

    cloneID = cloneID[cloneID.cloneId.isin(clones.cloneId)].copy()

    ncdr = sort_cnt(clones.aaSeqImputedCDR3)
    repncdr = ncdr[ncdr[1] > 1][0].to_list()
    reclone = clones[clones.aaSeqImputedCDR3.isin(repncdr)].sort_values(
        by=["aaSeqImputedCDR3", "cloneCount"], ascending=False
    )

    maps = {}
    for rep in repncdr:
        dd = reclone[reclone.aaSeqImputedCDR3 == rep].index.tolist()
        for idx in dd[1:]:
            maps[idx] = dd[0]

    cloneID.cloneId = cloneID.cloneId.apply(
        lambda x: maps[x] if x in maps.keys() else x
    ).copy()
    clones.drop_duplicates(subset="aaSeqImputedCDR3", keep="first", inplace=True)
    cloneID.set_index("descrsR1", inplace=True)
    cloneID = cloneID[["chains", "cloneId"]]

    print("clone filtering finished")
    clones.to_csv(f"{outdir}/{sample}_clones_filtered.csv")
    cloneID.to_csv(f"{outdir}/{sample}_cloneID_filtered.csv")
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_cloneID_filtered.csv"])
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_clones_filtered.csv"])


def write_bc_slideseq(sample, outdir, bc_file):
    left = 56 - 41
    right = 56 - 32
    barcodes = pd.read_table(bc_file, names=["bc"])

    if "BeadBarcodes" in bc_file:
        bcs = np.unique(barcodes.bc.apply(lambda x: "".join(x.split(","))).to_list())

    if "matched" in bc_file:
        bcs = np.unique(barcodes.bc.apply(lambda x: x.split("-")[0]).to_list())

    bc32base = ["N" * left + b[:8] + linker + b[8:] + "N" * right for b in bcs]

    # with open(f'{outdir}/{sample}_bcreads_pad_{left}_{right}.fasta', 'w') as f:
    with open(f"{outdir}/{sample}_bcreads.fasta", "w") as f:
        for i, b in enumerate(bc32base):
            f.write(f">{bcs[i]}\n")
            f.write(f"{bc32base[i]}\n")


def write_bc_5p10X(sample, outdir, bc_file):
    bcs = pd.read_table(bc_file, names=["bc"])
    bcs = bcs.bc.apply(lambda x: x.split("-")[0]).to_list()
    left = 30
    right = 40
    bc_pad = ["N" * left + b + "N" * right for b in bcs]

    with open(f"{outdir}/{sample}_bcreads.fasta", "w") as f:
        for i, b in enumerate(bc_pad):
            f.write(f">{bcs[i]}\n")
            f.write(f"{bc_pad[i]}\n")


def process_matching_slideseq_XCR(sample, outdir, cloneID):
    tot = 0
    readIDs = []
    all_AS = []
    bcs = []
    umis = []
    bad_bc = []
    samfile = pysam.AlignmentFile(f"{outdir}/{sample}_matching.sam", "r", threads=8)
    print("matching bam processing")
    for read in samfile.fetch():
        tot += 1
        AS = read.get_tag("AS")
        all_AS.append([AS, read.flag])
        if AS >= 28 and read.flag == 0:
            # reads.append(read)
            name = read.query_name

            bc = read.reference_name
            # print(read,bc,'\n')
            try:
                pairs = np.array(read.aligned_pairs)
                umi_s = dict(zip(pairs[:, 1], pairs[:, 0]))[
                    47
                ]  # left pad value (15) +32 = 47
                umi = read.query[umi_s : umi_s + 9]
                # print(read,bc,umi,'\n')
                if len(umi) == 9:
                    umis.append(umi)
                    bcs.append(bc)
                    readIDs.append(name)
            except:
                bad_bc.append(bc)
                # pass
        if tot % 50000 == 0:
            print(tot, " records processed")
    print("number of short UMI reads = ", len(bad_bc))
    set(bcs)

    all_AS = np.array(all_AS)
    scores = sort_cnt(all_AS[all_AS[:, 1] == 0][:, 0])
    scores.columns = ["score", "count"]

    plt.rcParams["figure.figsize"] = (5, 4)
    plt.rcParams["ytick.right"] = True
    plt.rcParams["ytick.labelright"] = True
    sns.barplot(data=scores[scores.score > 20], x="score", y="count")
    plt.savefig(f"{outdir}/{sample}_barcode_scores.pdf", bbox_inches="tight")
    print("clone_bcumi producing")
    bcumiID = np.array([readIDs, bcs, umis]).T
    bcumiID_df = pd.DataFrame(bcumiID, columns=["ID", "bc", "umi"])
    bcumiID_df.set_index("ID", inplace=True)
    merged_IDs = pd.merge(
        cloneID, bcumiID_df, how="inner", left_index=True, right_index=True
    )
    merged_IDs = merged_IDs.sort_values(by=["cloneId", "bc", "umi"])
    merged_IDs.to_csv(f"{outdir}/{sample}_clone_bcumi.csv", index=None)

    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_clone_bcumi.csv"])


def process_matching_5p10X(sample, outdir):
    tot = 0
    all_AS = []
    reads = []
    bad_bc = []
    read_bcumi_dic = {}
    bcumi_dic = {}

    samfile = pysam.AlignmentFile(f"{outdir}/{sample}_matching.sam", "rb", threads=8)

    for read in samfile.fetch():
        tot += 1
        AS = read.get_tag("AS")
        all_AS.append([AS, read.flag])
        if AS >= 14 and read.flag == 0:
            # print(read)
            name = read.query_name
            bc = read.reference_name
            seq = read.query
            try:
                pairs = np.array(read.aligned_pairs)
                pair_dic = dict(zip(pairs[:, 1], pairs[:, 0]))
                # rstart.append(pair_dic[54])
                umi = seq[pair_dic[46] : pair_dic[46] + 10]  # left pad +16 = 30 +16
            except:
                umi = "N"
            try:
                bcumi_dic[bc].append(umi)
            except:
                bcumi_dic[bc] = [umi]
            if len(umi) < 10:
                bad_bc.append(bc)
            else:
                read_bcumi_dic[name] = [bc, umi]
                # print(bc,umi,seq,'\n')
            # reads.append([read.seq,umi,bc]) #for debuging loc of umi bc
        if tot % 10000 == 0:
            print(tot, "barcode candidates processed")

    print("number of short UMI reads = ", len(bad_bc))

    all_AS = np.array(all_AS)

    scores = sort_cnt(all_AS[all_AS[:, 1] == 0][:, 0])
    scores.columns = ["score", "count"]
    scores.to_csv(f"{outdir}/{sample}_barcode_scores.csv", index=None)
    plt.rcParams["figure.figsize"] = (5, 3)
    plt.rcParams["ytick.right"] = True
    plt.rcParams["ytick.labelright"] = True
    sns.barplot(data=scores[scores.score > 8], x="score", y="count")
    plt.savefig(f"{outdir}/{sample}_barcode_scores.pdf", bbox_inches="tight")
    plt.close()

    read_keys = list(read_bcumi_dic.keys())
    bcumi_dic = {}
    read_bc_umi_trns_dic = {}
    for i, read in enumerate(read_keys):
        bcumi = read_bcumi_dic[read]
        # rname=read.split('_')[0]
        # rname='_'.join(read.split('_')[:-1])
        trns = read.split("_")[4]
        # trns_org=read.split('_')[-2]
        # if trns==trns_org:
        CB = bcumi[0]
        UB = bcumi[1]
        try:
            bcumi_dic[CB].append(UB)
        except:
            bcumi_dic[CB] = [UB]
        read_bc_umi_trns_dic[read] = [CB, UB, trns]
        if i % 20000 == 0 and i > 0:
            print(i, "Read-BC-UMI-Transcript tuples saved")
        # if i>10000: break

    umis = []
    reads = []
    for i, bc in enumerate(bcumi_dic.keys()):
        umis.append(np.unique(bcumi_dic[bc]).shape[0])
        reads.append(len(bcumi_dic[bc]))
        # if i%1000==0: print(i,'barcodes processed')

    bcumi_dedup = pd.DataFrame([bcumi_dic.keys(), umis, reads]).T
    bcumi_dedup.columns = ["bc", "umi_cnt", "read_cnt"]
    bcumi_dedup.umi_cnt = bcumi_dedup.umi_cnt.astype("int")
    bcumi_dedup.read_cnt = bcumi_dedup.read_cnt.astype("int")
    # bcumi_dedup=bcumi_dedup.sort_values(by='read_cnt',ascending=False)
    bcumi_dedup = bcumi_dedup.sort_values(by="umi_cnt", ascending=False)
    bcumi_dedup.set_index("bc", inplace=True)
    bcumi_dedup["dup_rate"] = bcumi_dedup.read_cnt / bcumi_dedup.umi_cnt

    bcumi_dedup[bcumi_dedup.umi_cnt > 0].to_csv(
        f"{outdir}/{sample}_bcumi_dedup.csv"
    )  # [bcumi_dedup.umi_cnt>100].shape

    # logcnt=np.log10(bcumi_dedup.read_cnt)
    # logcntumi=np.log10(bcumi_dedup.umi_cnt)
    # data=pd.concat([logcntumi,logcnt],axis=1)
    # data=data[data.umi_cnt>0]
    # plt.rcParams['figure.figsize'] = (5, 3)
    # sns.kdeplot(data=data);
    # plt.savefig(f'{outdir}/{sample}_bcumi.pdf',bbox_inches='tight');
    # plt.close()

    plt.rcParams["figure.figsize"] = (5, 5)
    # plt.plot(x=np.log10(np.arange(1,len(bcumi_dedup.umi_cnt)+1)),y=np.log10(bcumi_dedup.umi_cnt));
    plt.plot(
        np.log10(np.arange(1, len(bcumi_dedup.umi_cnt) + 1)),
        np.log10(bcumi_dedup.umi_cnt),
    )
    plt.ylabel("log10 UMI counts")
    plt.xlabel("log10 cell rank")
    plt.title("library knee plot")
    plt.savefig(f"{outdir}/{sample}_knee.pdf", bbox_inches="tight")
    plt.close()
    # subprocess.call([ 'pigz', '-f',  f'{outdir}/{sample}_clone_bcumi.csv'])

    sam_tag = f"{outdir}/{sample}_genome_tagged.bam"
    sam = f"{outdir}/{sample}_genome.bam"

    samfile = pysam.AlignmentFile(sam, "rb", threads=8)
    tagged_bam = pysam.AlignmentFile(sam_tag, "wb", template=samfile, threads=8)

    i = 0
    all_trns = []
    for read in samfile.fetch():
        name = read.qname  # .split('_')[0]
        # names.append(name)
        if read_bc_umi_trns_dic.get(name) is not None and read.flag < 20:
            i += 1
            bcumi_trns = read_bc_umi_trns_dic[name]
            read.set_tag("CB", bcumi_trns[0])
            read.set_tag("UB", bcumi_trns[1])
            read.set_tag("XT", bcumi_trns[2])
            all_trns.append(bcumi_trns[2])
            tagged_bam.write(read)
            if i % 20000 == 0:
                print(i, "genomic records saved")

    tagged_bam.close()
    samfile.close()

    trns_df = sort_cnt(all_trns)
    trns_df.to_csv(f"{outdir}/{sample}_trns_ct.csv", index=None)


def process_matching_5p10XTCR(sample, outdir):
    tot = 0
    all_AS = []
    reads = []
    bad_bc = []
    read_bcumi_dic = {}
    bcumi_dic = {}
    readIDs = []

    samfile = pysam.AlignmentFile(f"{outdir}/{sample}_matching.sam", "rb", threads=8)

    for read in samfile.fetch():
        tot += 1
        AS = read.get_tag("AS")
        all_AS.append([AS, read.flag])
        if AS >= 14 and read.flag == 0:
            # print(read)
            name = read.query_name
            bc = read.reference_name
            seq = read.query
            try:
                pairs = np.array(read.aligned_pairs)
                pair_dic = dict(zip(pairs[:, 1], pairs[:, 0]))
                # rstart.append(pair_dic[54])
                umi = seq[pair_dic[46] : pair_dic[46] + 10]  # left pad +16 = 30 +16
            except:
                umi = "N"
            try:
                bcumi_dic[bc].append(umi)
            except:
                bcumi_dic[bc] = [umi]
            if len(umi) < 10:
                bad_bc.append(bc)
            else:
                read_bcumi_dic[name] = [bc, umi]
                readIDs.append([name, bc, umi])
                # print(bc,umi,seq,'\n')
                # reads.append([read.seq,umi,bc]) #for debuging loc of umi bc

        if tot % 10000 == 0:
            print(tot, "barcode candidates processed")

    print("number of short UMI reads = ", len(bad_bc))

    all_AS = np.array(all_AS)

    scores = sort_cnt(all_AS[all_AS[:, 1] == 0][:, 0])
    scores.columns = ["score", "count"]
    scores.to_csv(f"{outdir}/{sample}_barcode_scores.csv", index=None)
    plt.rcParams["figure.figsize"] = (5, 3)
    plt.rcParams["ytick.right"] = True
    plt.rcParams["ytick.labelright"] = True
    sns.barplot(data=scores[scores.score > 8], x="score", y="count")
    plt.savefig(f"{outdir}/{sample}_barcode_scores.pdf", bbox_inches="tight")
    plt.close()

    read_keys = list(read_bcumi_dic.keys())
    bcumi_dic = {}
    read_bc_umi_trns_dic = {}
    for i, read in enumerate(read_keys):
        bcumi = read_bcumi_dic[read]
        # rname=read.split('_')[0]
        # rname='_'.join(read.split('_')[:-1])
        trns = read.split("_")[4]
        # trns_org=read.split('_')[-2]
        # if trns==trns_org:
        CB = bcumi[0]
        UB = bcumi[1]
        try:
            bcumi_dic[CB].append(UB)
        except:
            bcumi_dic[CB] = [UB]
        read_bc_umi_trns_dic[read] = [CB, UB, trns]
        if i % 20000 == 0 and i > 0:
            print(i, "Read-BC-UMI-Transcript tuples saved")
        # if i>10000: break

    umis = []
    reads = []
    for i, bc in enumerate(bcumi_dic.keys()):
        umis.append(np.unique(bcumi_dic[bc]).shape[0])
        reads.append(len(bcumi_dic[bc]))
        # if i%1000==0: print(i,'barcodes processed')

    bcumi_dedup = pd.DataFrame([bcumi_dic.keys(), umis, reads]).T
    bcumi_dedup.columns = ["bc", "umi_cnt", "read_cnt"]
    bcumi_dedup.umi_cnt = bcumi_dedup.umi_cnt.astype("int")
    bcumi_dedup.read_cnt = bcumi_dedup.read_cnt.astype("int")
    # bcumi_dedup=bcumi_dedup.sort_values(by='read_cnt',ascending=False)
    bcumi_dedup = bcumi_dedup.sort_values(by="umi_cnt", ascending=False)
    bcumi_dedup.set_index("bc", inplace=True)
    bcumi_dedup["dup_rate"] = bcumi_dedup.read_cnt / bcumi_dedup.umi_cnt

    bcumi_dedup[bcumi_dedup.umi_cnt > 0].to_csv(
        f"{outdir}/{sample}_bcumi_dedup.csv"
    )  # [bcumi_dedup.umi_cnt>100].shape
    # subprocess.call([ 'pigz', '-f',  f'{outdir}/{sample}_clone_bcumi.csv'])

    # logcnt=np.log10(bcumi_dedup.read_cnt)
    # logcntumi=np.log10(bcumi_dedup.umi_cnt)
    # data=pd.concat([logcntumi,logcnt],axis=1)
    # data=data[data.umi_cnt>0]
    # plt.rcParams['figure.figsize'] = (5, 3)
    # sns.kdeplot(data=data);
    # plt.savefig(f'{outdir}/{sample}_bcumi.pdf',bbox_inches='tight');
    # plt.close()

    # plt.rcParams['figure.figsize'] = (5, 5)
    # plt.plot(x=np.log10(np.arange(1,len(bcumi_dedup.umi_cnt)+1)),y=np.log10(bcumi_dedup.umi_cnt));
    plt.plot(
        np.log10(np.arange(1, len(bcumi_dedup.umi_cnt) + 1)),
        np.log10(bcumi_dedup.umi_cnt),
    )
    plt.ylabel("log10 UMI counts")
    plt.xlabel("log10 cell rank")
    plt.title("library knee plot")
    plt.savefig(f"{outdir}/{sample}_knee_UMI.pdf", bbox_inches="tight")
    plt.close()

    # plt.rcParams['figure.figsize'] = (5, 5)
    # plt.plot(x=np.log10(np.arange(1,len(bcumi_dedup.umi_cnt)+1)),y=np.log10(bcumi_dedup.umi_cnt));
    bcumi_dedup = bcumi_dedup.sort_values(by="read_cnt", ascending=False)
    plt.plot(
        np.log10(np.arange(1, len(bcumi_dedup.read_cnt) + 1)),
        np.log10(bcumi_dedup.read_cnt),
    )
    plt.ylabel("log10 read counts")
    plt.xlabel("log10 cell rank")
    plt.title("library knee plot")
    plt.savefig(f"{outdir}/{sample}_knee_reads.pdf", bbox_inches="tight")
    plt.close()

    print("clone_bcumi producing")

    cloneID = pd.read_csv(f"{outdir}/{sample}_cloneID_filtered.csv.gz", index_col=0)
    bcumiID = np.array(readIDs)
    # return(readIDs)

    bcumiID_df = pd.DataFrame(bcumiID, columns=["ID", "bc", "umi"])
    bcumiID_df.set_index("ID", inplace=True)
    merged_IDs = pd.merge(
        bcumiID_df, cloneID, how="inner", left_index=True, right_index=True
    )
    merged_IDs = merged_IDs.sort_values(by=["cloneId", "bc", "umi"])
    merged_IDs.to_csv(f"{outdir}/{sample}_clone_bcumi.csv", index=None)

    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_clone_bcumi.csv"])


def decon_3p10XTCR_nuc(sample, outdir):
    tot = 0
    short_bc = 0
    file = f"{outdir}/{sample}_trns.sam"

    samfile = pysam.AlignmentFile(f"{file}", "r", threads=8)

    # const=6*'A'+28*'N'+'AGATCGGAAGAGCGTCGTGT'

    const = "AGATCGGAAGAGCGTCGTGT"

    r_search = 2000  # search length into 3' direction of 3' softclip of V gene (to search for BC-UMI)

    f1 = open(f"{outdir}/{sample}_VDJ.fastq", "w")
    f2 = open(f"{outdir}/{sample}_BCUMI.fasta", "w")
    print(f"{outdir}/{sample}_VDJ.fastq")
    for read in samfile.fetch():
        rlen = read.rlen
        seq = read.seq
        flag = read.flag
        qstrt = read.query_alignment_start
        qend = read.query_alignment_end
        read.qname[-10:]
        trans = read.reference_name.split("-")[0]
        ref_s = read.reference_start
        ref_e = read.reference_end
        ref_e - ref_s

        rclip = 100  # PART INTO 3' OF V GENE (SHOULD INCLUDE CDR3 AND PARTIAL C GENE or C-J INTRON)

        if rlen - qend > r_search:
            end_qu = seq[qend : qend + r_search]
        else:
            end_qu = seq[qend:]

        # sub_end = qstrt+rclip # for including part of C gene (when aligning to C genes)
        # sub_end = qend+rclip # for including CDR3 J and part of C-J-intron or C gene

        if rlen - qend > rclip:
            sub_end = qend + rclip
        else:
            sub_end = rlen

        # if qstrt>lclip: sub_strt = qstrt-lclip
        # else: sub_strt = 0
        sub_strt = qstrt

        sub_seq = seq[sub_strt:sub_end]
        sub_qual = read.qual[sub_strt:sub_end]

        rclip_truseq = 40
        lclip_truseq = 10
        step = 200

        # subs.append(len(sub_seq))
        newnamef = f">{read.qname}_{sample}_{sub_strt}_{sub_end}_{flag}_{trans}"

        if len(sub_seq) > 100:
            f1.write(f"@{read.qname}_{sample}_{sub_strt}_{sub_end}_{flag}_{trans}\n")
            f1.write(f"{sub_seq}\n")
            f1.write("+\n")
            f1.write(f"{sub_qual}\n")

            number_steps = int(len(end_qu) / step)
            # print(number_steps,len(end_qu))

            for i in range(number_steps + 1):
                w = end_qu[step * i : step * (i + 1) + 70]
                # print(w)
                ed = edlib.align(const, w, "HW", "locations", 2)

                # print(i,ed)
                # print(w)

                if ed["editDistance"] > -1 and ed["editDistance"] < 3:
                    # print(ed)
                    start = ed["locations"][0][0] + 200 * i
                    end = ed["locations"][0][1] + 200 * i

                    if start - rclip_truseq < 0:
                        pass
                    else:
                        start - rclip_truseq
                    end + lclip_truseq
                    # bcumi=rev(end_qu[upstart:upend])

                    bcumi = rev(end_qu[start - 35 : end - 12])
                    # polyA=dd[:upstart+5]
                    # c_hangs.append(select)
                    # polyAs.append(polyA)
                    # print('>bcumui=',bcumi)
                    # eds.append([ed['editDistance'],bcumi])

                    # eds.append(ed['editDistance'])
                    # newnames.append(newnamef)
                    if len(bcumi) > 30:
                        f2.write(f"{newnamef}\n")
                        f2.write(f"{bcumi}\n")
                    else:
                        short_bc += 1
                        # print('short_BCUMI')
                    break
            # print('\n')

        # ed=edlib.align(const, end_qu,'HW','locations',5,ad_seq)
        # dist=ed['editDistance']
        # eds.append(dist)
        # newname=f'{rname}_q{qlen}_d{dist}_s{sub_strt}_e{sub_end}_f{flag}_{trans}'
        #
        # if dist>-1 and dist<6 and len(sub_seq)>100 and qlen >100:
        #    f1.write(f'@{newname}\n')
        #    f1.write(f'{sub_seq}\n')
        #    f1.write('+\n')
        #    f1.write(f'{sub_qual}\n')
        #    #bcumi=rev(end_qu[ed['locations'][0][0]:ed['locations'][0][1]])[14:]
        #    f2.write(f'>{newname}\n')
        #    f2.write(f'{bcumi}\n')
        tot += 1
        if tot % 50000 == 0:
            print(tot, " records processed")
        # if tot>5000:break
    samfile.close()
    f1.close()
    f2.close()

    print(short_bc)

    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_VDJ.fastq"])
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_BCUMI.fasta"])

    # sort_cnt(eds).to_csv(f'{outdir}/{sample}_eds.csv')
    # return(eds)


def write_bc_3p10XTCR_nuc(sample, outdir, bc_file):
    import scanpy as sc

    adata = sc.read_10x_h5(bc_file)

    sc.pp.filter_cells(adata, min_genes=4)

    bcs = [a[:16] for a in adata.obs.index]

    left = 16
    right = 28
    bc_pad = ["N" * left + b + "N" * right for b in bcs]
    # with open(f'./outputs/737k_pad_{left}_{}.fasta', 'w') as f:
    with open(f"{outdir}/{sample}_bcreads.fasta", "w") as f:
        for i, b in enumerate(bc_pad):
            f.write(f">{bcs[i]}\n")
            f.write(f"{bc_pad[i]}\n")


def process_matching_3p10XTCR_nuc(sample, outdir):
    tot = 0
    all_AS = []
    reads = []
    bad_bc = []
    read_bcumi_dic = {}
    bcumi_dic = {}
    readIDs = []

    samfile = pysam.AlignmentFile(f"{outdir}/{sample}_matching.sam", "rb", threads=8)

    for read in samfile.fetch():
        tot += 1
        AS = read.get_tag("AS")
        all_AS.append([AS, read.flag])
        if AS >= 14 and read.flag == 0:
            # print(read)
            name = read.query_name
            bc = read.reference_name
            seq = read.query
            try:
                pairs = np.array(read.aligned_pairs)
                pair_dic = dict(zip(pairs[:, 1], pairs[:, 0]))
                # rstart.append(pair_dic[54])
                umi = seq[pair_dic[32] : pair_dic[32] + 12]  # left pad +16 = 16 +16
            except:
                umi = "N"
            try:
                bcumi_dic[bc].append(umi)
            except:
                bcumi_dic[bc] = [umi]
            if len(umi) < 12:
                bad_bc.append(bc)
            else:
                read_bcumi_dic[name] = [bc, umi]
                readIDs.append([name, bc, umi])
                # print(bc,umi,seq,'\n')
                # reads.append([read.seq,umi,bc]) #for debuging loc of umi bc

        if tot % 10000 == 0:
            print(tot, "barcode candidates processed")

    print("number of short UMI reads = ", len(bad_bc))

    all_AS = np.array(all_AS)

    scores = sort_cnt(all_AS[all_AS[:, 1] == 0][:, 0])
    scores.columns = ["score", "count"]
    scores.to_csv(f"{outdir}/{sample}_barcode_scores.csv", index=None)
    plt.rcParams["figure.figsize"] = (5, 3)
    plt.rcParams["ytick.right"] = True
    plt.rcParams["ytick.labelright"] = True
    sns.barplot(data=scores[scores.score > 8], x="score", y="count")
    plt.savefig(f"{outdir}/{sample}_barcode_scores.pdf", bbox_inches="tight")
    plt.close()

    read_keys = list(read_bcumi_dic.keys())
    bcumi_dic = {}
    read_bc_umi_trns_dic = {}
    for i, read in enumerate(read_keys):
        bcumi = read_bcumi_dic[read]
        # rname=read.split('_')[0]
        # rname='_'.join(read.split('_')[:-1])
        trns = read.split("_")[4]
        # trns_org=read.split('_')[-2]
        # if trns==trns_org:
        CB = bcumi[0]
        UB = bcumi[1]
        try:
            bcumi_dic[CB].append(UB)
        except:
            bcumi_dic[CB] = [UB]
        read_bc_umi_trns_dic[read] = [CB, UB, trns]
        if i % 20000 == 0 and i > 0:
            print(i, "Read-BC-UMI-Transcript tuples saved")
        # if i>10000: break

    umis = []
    reads = []
    for i, bc in enumerate(bcumi_dic.keys()):
        umis.append(np.unique(bcumi_dic[bc]).shape[0])
        reads.append(len(bcumi_dic[bc]))
        # if i%1000==0: print(i,'barcodes processed')

    bcumi_dedup = pd.DataFrame([bcumi_dic.keys(), umis, reads]).T
    bcumi_dedup.columns = ["bc", "umi_cnt", "read_cnt"]
    bcumi_dedup.umi_cnt = bcumi_dedup.umi_cnt.astype("int")
    bcumi_dedup.read_cnt = bcumi_dedup.read_cnt.astype("int")
    # bcumi_dedup=bcumi_dedup.sort_values(by='read_cnt',ascending=False)
    bcumi_dedup = bcumi_dedup.sort_values(by="umi_cnt", ascending=False)
    bcumi_dedup.set_index("bc", inplace=True)
    bcumi_dedup["dup_rate"] = bcumi_dedup.read_cnt / bcumi_dedup.umi_cnt

    bcumi_dedup[bcumi_dedup.umi_cnt > 0].to_csv(
        f"{outdir}/{sample}_bcumi_dedup.csv"
    )  # [bcumi_dedup.umi_cnt>100].shape
    # subprocess.call([ 'pigz', '-f',  f'{outdir}/{sample}_clone_bcumi.csv'])

    # logcnt=np.log10(bcumi_dedup.read_cnt)
    # logcntumi=np.log10(bcumi_dedup.umi_cnt)
    # data=pd.concat([logcntumi,logcnt],axis=1)
    # data=data[data.umi_cnt>0]
    # plt.rcParams['figure.figsize'] = (5, 3)
    # sns.kdeplot(data=data);
    # plt.savefig(f'{outdir}/{sample}_bcumi.pdf',bbox_inches='tight');
    # plt.close()

    # plt.rcParams['figure.figsize'] = (5, 5)
    # plt.plot(x=np.log10(np.arange(1,len(bcumi_dedup.umi_cnt)+1)),y=np.log10(bcumi_dedup.umi_cnt));
    plt.plot(
        np.log10(np.arange(1, len(bcumi_dedup.umi_cnt) + 1)),
        np.log10(bcumi_dedup.umi_cnt),
    )
    plt.ylabel("log10 UMI counts")
    plt.xlabel("log10 cell rank")
    plt.title("library knee plot")
    plt.savefig(f"{outdir}/{sample}_knee_UMI.pdf", bbox_inches="tight")
    plt.close()

    # plt.rcParams['figure.figsize'] = (5, 5)
    # plt.plot(x=np.log10(np.arange(1,len(bcumi_dedup.umi_cnt)+1)),y=np.log10(bcumi_dedup.umi_cnt));
    bcumi_dedup = bcumi_dedup.sort_values(by="read_cnt", ascending=False)
    plt.plot(
        np.log10(np.arange(1, len(bcumi_dedup.read_cnt) + 1)),
        np.log10(bcumi_dedup.read_cnt),
    )
    plt.ylabel("log10 read counts")
    plt.xlabel("log10 cell rank")
    plt.title("library knee plot")
    plt.savefig(f"{outdir}/{sample}_knee_reads.pdf", bbox_inches="tight")
    plt.close()

    print("clone_bcumi producing")

    cloneID = pd.read_csv(f"{outdir}/{sample}_cloneID_filtered.csv.gz", index_col=0)
    bcumiID = np.array(readIDs)

    bcumiID_df = pd.DataFrame(bcumiID, columns=["ID", "bc", "umi"])
    bcumiID_df.set_index("ID", inplace=True)
    merged_IDs = pd.merge(
        bcumiID_df, cloneID, how="inner", left_index=True, right_index=True
    )
    merged_IDs = merged_IDs.sort_values(by=["cloneId", "bc", "umi"])
    merged_IDs.to_csv(f"{outdir}/{sample}_clone_bcumi.csv", index=None)

    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_clone_bcumi.csv"])


def decon_3p10XGEX(sample, outdir):
    tot = 0
    short_BC = 0
    file = f"{outdir}/{sample}_trns.sam"
    samfile = pysam.AlignmentFile(f"{file}", "r")
    bc_count_dict = {}
    bc_count_json = f"{outdir}/{sample}_bc_count.json"
    if os.path.isfile(bc_count_json):
        print(bc_count_json, " exists, skip")
        return

    # const=6*'A'+28*'N'+'AGATCGGAAGAGCGTCGTGT'

    const = "AGATCGGAAGAGCGTCGTGT"
    r_search = 700  # search length into 3' direction of 3' softclip of V gene (to search for BC-UMI)

    l_search = 500

    rclip = 1  # to keep in deconcat file

    lclip = 1  # to keep in deconcat file

    f1 = open(f"{outdir}/{sample}_deconcat.fastq", "w")
    f2 = open(f"{outdir}/{sample}_BCUMI.fasta", "w")

    for read in samfile.fetch():
        rlen = read.rlen
        seq = read.seq
        flag = read.flag
        qstrt = read.query_alignment_start
        qend = read.query_alignment_end
        trans = read.reference_name  # .split('-')[0]
        ref_s = read.reference_start
        ref_e = read.reference_end
        ref_e - ref_s

        # PART TO 3' OF V GENE (SHOULD INCLUDE CDR3 AND PARTIAL C GENE)

        if rlen - qend > r_search:
            end_qu = seq[qend - 70 : qend + r_search]
        else:
            end_qu = seq[qend - 70 :]

        if qstrt < l_search:
            seq[:qstrt]
        else:
            seq[qstrt - l_search : qstrt]

        if rlen - qend > rclip:
            sub_end = qend + rclip
        else:
            sub_end = rlen

        if qstrt < lclip:
            sub_strt = 0
        else:
            sub_strt = qstrt - lclip

        sub_seq = seq[sub_strt:sub_end]  # subset only the selected parts to save
        sub_qual = read.qual[sub_strt:sub_end]

        if flag == 16 or flag == 2064:
            qstrt_mod = rlen - qend
            qend_mod = rlen - qstrt
        else:
            qstrt_mod = qstrt
            qend_mod = qend

        newnamef = f"{read.qname}_{qstrt_mod}_{qend_mod}_{flag}_{trans}"

        if len(sub_seq) > 50:
            f1.write(f"@{newnamef}\n")
            f1.write(f"{sub_seq}\n")
            f1.write("+\n")
            f1.write(f"{sub_qual}\n")

            step = 200
            number_steps = int(len(end_qu) / step)

            # print(number_steps,len(end_qu))

            for i in range(number_steps + 1):
                w = end_qu[step * i : step * (i + 1) + 70]
                # print(w)
                ed = edlib.align(const, w, "HW", "locations", 3)

                if ed["editDistance"] > -1 and ed["editDistance"] < 4:
                    # print(ed)
                    start = ed["locations"][0][0] + 200 * i
                    ed["locations"][0][1] + 200 * i

                    bcumi = rev(end_qu[start - 16 - 12 - 4 : start + 3])

                    seq_counter(bc_count_dict, bcumi[3 : 3 + 16])

                    if len(bcumi) > 30:
                        f2.write(f">{newnamef}\n")
                        f2.write(f"{bcumi}\n")
                    else:
                        short_BC += 1
                    break
                # if ed['editDistance']>=3 or ed['editDistance']<0:
                #    eds.append([i,ed['editDistance']])
            # print('\n')

        tot += 1
        if tot % 4e5 == 0:
            print(tot, " records processed")

        if tot > 1e9:
            break
    samfile.close()
    f1.close()
    f2.close()

    print("number of short BCUMIs = ", short_BC)

    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_deconcat.fastq"])
    subprocess.call(["pigz", "-f", f"{outdir}/{sample}_BCUMI.fasta"])

    with open(bc_count_json, "w") as json_file:
        json.dump(bc_count_dict, json_file)

    # read_cnt=pd.Series(bc_count_dict)
    # read_cnt.to_csv(f'{outdir}/{sample}_bc_count.csv')
    # pd.DataFrame(eds).to_csv(f'{outdir}/{sample}_eds.csv')
    # sort_cnt(eds).to_csv(f'{outdir}/{sample}_eds.csv')


def write_bc_3p10XGEX(sample, outdir, barcodes):
    bcreads_fasta = f"{outdir}/{sample}_bcreads.fasta"

    if os.path.isfile(bcreads_fasta):
        print(bcreads_fasta, " exists, skip")
        return

    jsons = sorted(
        [f for f in os.listdir(f"{outdir}/split/") if f.endswith("bc_count.json")]
    )

    bc_agg = {}
    for i in range(len(jsons)):
        with open(f"{outdir}/split/{jsons[i]}", "r") as json_file:
            bc_count = json.load(json_file)
            for k in bc_count:
                if bc_agg.get(k) is not None:
                    bc_agg[k] += bc_count[k]
                else:
                    bc_agg[k] = bc_count[k]

    read_cnt = pd.Series(bc_agg)
    read_cnt.name = "read_count"
    print(read_cnt.shape)
    read_cnt.to_csv(f"{outdir}/{sample}_bc_read_count.csv")
    # raw_bcs = pd.read_csv(f'{outdir}/{sample}_bc_read_count.csv',index_col=0)

    raw_bcs = read_cnt[read_cnt > 20].copy()

    print(raw_bcs.shape)

    white_bcs = pd.read_table(barcodes, names=["bc"])
    # print(white_bcs.shape)

    shared = set(white_bcs["bc"]) & set(raw_bcs.index)
    bcs = raw_bcs[raw_bcs.index.isin(shared)].index
    print(len(bcs))

    # 3 on Truseq, 4 on polyT
    left = 3 + 1
    right = 12 + 4 + 1
    bc_pad = ["N" * left + b + "N" * right for b in bcs]

    with open(bcreads_fasta, "w") as f:
        for i, b in enumerate(bc_pad):
            f.write(f">{bcs[i]}\n")
            f.write(f"{bc_pad[i]}\n")


def process_matching_3p10XGEX(sample, outdir):
    tot = 0
    all_AS = []
    bad_bc = []

    samfile = pysam.AlignmentFile(f"{outdir}/{sample}_matching.sam", "r")

    quad_dict = {}
    quads_json = f"{outdir}/{sample}_quads.json"

    if os.path.isfile(quads_json):
        print(quads_json, " exists, skip")
        return

    for read in samfile.fetch():
        tot += 1
        AS = read.get_tag("AS")
        all_AS.append([AS, read.flag])
        if AS >= 14 and read.flag == 0:
            # print(read)
            name = read.query_name

            trns = "_".join(name.split("/ccs_")[-1].split("_")[3:])

            if "|" in trns:
                trns = trns.split("|")[-1].split("_")[0]

            bc = read.reference_name
            seq = read.query
            umi_start = 4 + 16  # left pad + 16 = 4 + 16
            umi_length = 12
            try:
                pairs = np.array(read.aligned_pairs)
                pair_dic = dict(zip(pairs[:, 1], pairs[:, 0]))
                umi = seq[
                    pair_dic[umi_start] : pair_dic[umi_start] + umi_length
                ]  # left pad + 16 = 4 + 16
            except:
                umi = "N"

            if len(umi) < 12:
                bad_bc.append(bc)
            else:
                quad_dict_store(quad_dict, bc, [umi, trns])

        if tot % 400000 == 0:
            print(tot, "barcode candidates processed")
        if tot > 1e20:
            break

    print("number of short UMI reads = ", len(bad_bc))

    with open(quads_json, "w") as json_file:
        json.dump(quad_dict, json_file)

    all_AS = np.array(all_AS)

    scores = sort_cnt(all_AS[all_AS[:, 1] == 0][:, 0])
    scores.columns = ["score", "count"]
    scores.to_csv(f"{outdir}/{sample}_barcode_scores.csv", index=None)


def make_count_mtx_3p10XGEX(sample, outdir):
    mtx_file = f"{outdir}/{sample}_gex.mtx.gz"

    if os.path.isfile(mtx_file):
        print(mtx_file, " exists, skip")
        return

    jsons = sorted(
        [f for f in os.listdir(f"{outdir}/split/") if f.endswith("quads.json")]
    )

    data_agg = {}

    """
    will implement parallel version
    """

    for i in tqdm(range(len(jsons))):
        with open(f"{dir_split}{jsons[i]}", "r") as json_file:
            data_sub = json.load(json_file)
            print(jsons[i], len(data_sub))
            for k in data_sub:
                if data_agg.get(k) is not None:
                    data_agg[k].extend(data_sub[k])
                else:
                    data_agg[k] = data_sub[k]
