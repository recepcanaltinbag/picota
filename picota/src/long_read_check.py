#!/usr/bin/env python3
import os
import subprocess
import pysam
from pathlib import Path
import matplotlib.pyplot as plt
import re
from Bio import SeqIO
# -----------------------------
# PARAMETRELER
# -----------------------------
RAW_READS = "/media/lin-bio/back2/picota_test/SRR12917043/sra_files/SRR12917043_1.fastq"  # Ham nanopore okuma dosyası (tüm klasörler için aynı ise buraya yaz)
THREADS = 24
IDENTITY_THRESHOLD = 0.90  # %90 eşleşme
BASE_DIR = Path("/media/lin-bio/back2/picota_out/blasted/SRR12917032_k105")

MERGED_FASTA = BASE_DIR / "merged.fasta"
MERGED_BAM = BASE_DIR / "merged.sorted.bam"

# -----------------------------
# Yardımcı Fonksiyonlar
# -----------------------------
def merge_fastas(base_dir, merged_path):
    with open(merged_path, "w") as fout:
        for fasta_path in base_dir.rglob("*.fasta"):
            if fasta_path == MERGED_FASTA:
                continue  # merged fasta’yı atla
            folder = fasta_path.parent.name
            fasta_name = fasta_path.stem
            with open(fasta_path) as f:
                seq = []
                header = ""
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if seq:  # önceki sequence’i yaz
                            fout.write("".join(seq) + "\n")
                            seq = []
                        header = line[1:].replace(" ", "_")
                        fout.write(f">{folder}_{fasta_name}|{header}\n")
                    else:
                        seq.append(line)
                if seq:  # dosya sonundaki seq’i yaz
                    fout.write("".join(seq) + "\n")


def run_minimap2(ref_fasta, reads_fastq, out_bam):
    sam_path = out_bam.with_suffix(".sam")
    subprocess.run([
        "minimap2", "-ax", "map-ont", "--secondary=no",
        "-t", str(THREADS), str(ref_fasta), str(reads_fastq)
    ], stdout=open(sam_path, "w"), check=True)
    subprocess.run(["samtools", "view", "-bS", str(sam_path)], stdout=open(out_bam, "wb"), check=True)
    subprocess.run(["samtools", "sort", "-o", str(out_bam), str(out_bam)], check=True)
    subprocess.run(["samtools", "index", str(out_bam)], check=True)
    sam_path.unlink()

def extract_matching_reads_by_contig(bam_path, identity_threshold):
    bam = pysam.AlignmentFile(bam_path, "rb")
    contig_reads = dict()
    for aln in bam:
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        matches = sum(length for op,length in aln.cigartuples if op==0)
        aln_len = aln.query_alignment_length
        identity = matches / aln_len if aln_len>0 else 0
        if identity < identity_threshold:
            continue
        contig = aln.reference_name
        contig_reads.setdefault(contig, []).append({
            "read_id": aln.query_name,
            "start": aln.reference_start,
            "end": aln.reference_end,
            "len": aln.query_alignment_length,
            "aln": aln  # <-- işte buraya ekliyoruz
        })
    return contig_reads


from pysam import AlignedSegment

def get_softclip_lengths(aln: AlignedSegment):
    cigar = aln.cigartuples or []
    left_clip = cigar[0][1] if cigar and cigar[0][0] == 4 else 0
    right_clip = cigar[-1][1] if cigar and cigar[-1][0] == 4 else 0
    return left_clip, right_clip

def check_circularity(matching_reads, contig_length, softclip_frac=0.1):
    """
    Circular read tespiti:
    - Split alignment wrap-around
    - Tek parça wrap-around (end < start)
    - Çift taraflı softclip
    """
    circ_reads = []

    # Read ID'ye göre grupla
    reads_by_id = {}
    for r in matching_reads:
        reads_by_id.setdefault(r["read_id"], []).append(r)

    for read_id, parts in reads_by_id.items():
        if len(parts) == 1:
            r = parts[0]
            read_len = r["len"]

            # 1) Tek parça wrap-around
            if r["end"] < r["start"] and read_len > contig_length * 0.8:
                circ_reads.append(r)

            # 2) Başta başlıyor, sonda bitiyor
            elif r["start"] < contig_length * 0.1 and r["end"] > contig_length * 0.9:
                circ_reads.append(r)

            # 3) Çift taraflı softclip kontrolü
            elif isinstance(r["aln"], pysam.AlignedSegment):
                left_clip, right_clip = get_softclip_lengths(r["aln"])
                near_start = r["start"] < contig_length * 0.1
                near_end = r["end"] > contig_length * 0.9
                left_big = left_clip > read_len * softclip_frac
                right_big = right_clip > read_len * softclip_frac

                if left_big and near_end and right_big and near_start:
                    circ_reads.append(r)

        else:
            # Çok parçalı alignment wrap-around kontrolü
            for r1 in parts:
                for r2 in parts:
                    if r1 is r2:
                        continue
                    if r1["end"] > contig_length * 0.9 and r2["start"] < contig_length * 0.1:
                        total_cov = r1["len"] + r2["len"]
                        if total_cov > contig_length * 0.8:
                            circ_reads.extend(parts)
                            break

    return circ_reads



def plot_coverage(matching_reads, circular_reads, contig_length, output_png):
    coverage = [0]*contig_length
    circ_positions = set()
    for r in matching_reads:
        for pos in range(r["start"], min(r["end"], contig_length)):
            coverage[pos]+=1
        if r in circular_reads:
            for pos in range(r["start"], min(r["end"], contig_length)):
                circ_positions.add(pos)
    plt.figure(figsize=(10,4))
    plt.plot(range(contig_length), coverage, color="teal", label="Tüm eşleşmeler")
    if circular_reads:
        circ_cov = [coverage[i] if i in circ_positions else 0 for i in range(contig_length)]
        plt.plot(range(contig_length), circ_cov, color="red", linewidth=2, label="Circular adayları")
    plt.xlabel("Pozisyon (bp)")
    plt.ylabel("Coverage")
    plt.title("Coverage Plot")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_png)
    plt.close()

def plot_wraparound(circular_reads, contig_length, output_png):
    if not circular_reads:
        return
    plt.figure(figsize=(10,2))
    for r in circular_reads:
        plt.plot([r["start"], r["end"]], [1,1], color="red", linewidth=3)
    plt.xlim(0, contig_length)
    plt.yticks([])
    plt.xlabel("Dairesel Referans Pozisyonu (bp)")
    plt.title("Wrap-around hizalamaları (Circular adaylar)")
    plt.tight_layout()
    plt.savefig(output_png)
    plt.close()

def contig_length_from_merged(fasta_path):
    lengths = dict()
    with open(fasta_path) as f:
        seq = []
        header = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    lengths[header] = len("".join(seq))
                    seq = []
                header = line.strip()[1:]
            else:
                seq.append(line.strip())
        if seq:
            lengths[header] = len("".join(seq))
    return lengths

def extract_reads_fastq(txt_file, raw_fastq, out_fastq):
    matching_ids=set()
    with open(txt_file) as f:
        for line in f:
            matching_ids.add(line.split("\t")[0])
    with open(raw_fastq) as inf, open(out_fastq,"w") as outf:
        for rec in SeqIO.parse(inf,"fastq"):
            if rec.id in matching_ids:
                SeqIO.write(rec, outf, "fastq")

# -----------------------------
# ANA ÇALIŞMA
# -----------------------------
#merge_fastas(BASE_DIR, MERGED_FASTA)
#run_minimap2(MERGED_FASTA, RAW_READS, MERGED_BAM)

def safe_filename(name, max_len=60):
    name = re.sub(r'[^A-Za-z0-9_\-]', '_', name)
    if len(name) > max_len:
        name = name[:max_len]
    return name

contig_lengths = contig_length_from_merged(MERGED_FASTA)
contig_reads = extract_matching_reads_by_contig(MERGED_BAM, IDENTITY_THRESHOLD)

for contig, reads in contig_reads.items():
    folder_fasta, contig_id = contig.split("|",1)
    folder = BASE_DIR / folder_fasta.split("_")[0]
    os.makedirs(folder, exist_ok=True)


    fasta_name = "_".join(folder_fasta.split("_")[1:])
    contig_len = contig_lengths[contig]
    circular_reads = check_circularity(reads, contig_len)
    print(f"{contig} -> {len(reads)} eşleşme, {len(circular_reads)} circular aday")

    # Güvenli dosya isimleri
    safe_fasta_name = safe_filename(fasta_name)
    safe_contig_id = safe_filename(contig_id)
    coverage_png = folder/f"coverage_{safe_contig_id}.png"
    wrap_png = folder/f"wraparound_{safe_contig_id}.png"
    txt_file = folder/f"matching_reads_{safe_contig_id}.txt"
    fastq_out = folder/f"matched_reads_{safe_contig_id}.fastq"

    # Çıktılar

    plot_coverage(reads, circular_reads, contig_len, coverage_png)


    plot_wraparound(circular_reads, contig_len, wrap_png)
    with open(txt_file,"w") as out:
        for r in reads:
            circ_flag = "CIRCULAR" if r in circular_reads else ""
            out.write(f"{r['read_id']}\t{r['start']}\t{r['end']}\t{r['len']}\t{circ_flag}\n")
    extract_reads_fastq(txt_file, RAW_READS, fastq_out)

print("\nBitti ✅")