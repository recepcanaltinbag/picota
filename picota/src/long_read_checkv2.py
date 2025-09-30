import subprocess
import pysam
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
import math
import numpy as np


# === 1. Mapping ===
def run_minimap2(ref_fasta, fastq_file, bam_out="mapping.bam", threads=4, tmp_dir="tmp_mapping"):
    """Minimap2 ile mapping yapar ve sorted BAM döner. 
    Temp dosyalar ayrı klasörde tutulur. Var olanlar tekrar kullanılmaz."""

    os.makedirs(tmp_dir, exist_ok=True)

    sam_out = os.path.join(tmp_dir, os.path.basename(bam_out).replace(".bam", ".sam"))
    bam_out = os.path.join(tmp_dir, os.path.basename(bam_out))
    sorted_bam = bam_out.replace(".bam", "_sorted.bam")
    bai_out = sorted_bam + ".bai"

    # Eğer sorted bam + index varsa hiç çalıştırma
    if os.path.exists(sorted_bam) and os.path.exists(bai_out):
        print(f"[OK] {sorted_bam} ve {bai_out} bulundu, mapping atlandı.")
        return sorted_bam

    # Minimap2
    if not os.path.exists(sam_out):
        cmd = [
            "minimap2", "-ax", "map-ont", "-t", str(threads),
            ref_fasta, fastq_file
        ]
        print(f"[RUN] Minimap2: {' '.join(cmd)}")
        with open(sam_out, "w") as f:
            subprocess.run(cmd, stdout=f, check=True)

    # SAM -> BAM
    if not os.path.exists(bam_out):
        subprocess.run(["samtools", "view", "-bS", sam_out, "-o", bam_out], check=True)

    # BAM -> sorted BAM
    if not os.path.exists(sorted_bam):
        subprocess.run(["samtools", "sort", bam_out, "-o", sorted_bam], check=True)

    # BAM index
    if not os.path.exists(bai_out):
        subprocess.run(["samtools", "index", sorted_bam], check=True)

    return sorted_bam




# === 2. Circular read bulucu (Nanopore-friendly, SA tag düzeltmesi) ===
def find_circular_reads(bam_file, ref_name, ref_len, fastq_file, outdir, tolerance=0.50):
    """BAM dosyasından circular candidate read'leri bulur, split read'leri ve SA tag'leri dikkate alır"""
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Read sequence’leri saklamak için
    fastq_dict = {rec.id: rec for rec in SeqIO.parse(fastq_file, "fastq")}

    reads = {}
    for read in bam.fetch(ref_name):
        if read.query_name not in reads:
            reads[read.query_name] = []
        # alignment bilgilerini tuple olarak saklayalım: (start, end)
        reads[read.query_name].append((read.reference_start, read.reference_end))

        # Supplementary alignments (SA tag) varsa ekle
        if read.has_tag("SA"):
            sa_tag = read.get_tag("SA")
            for sa_entry in sa_tag.split(";"):
                if not sa_entry:
                    continue
                sa_fields = sa_entry.split(",")
                sa_start = int(sa_fields[1])
                sa_cigar = sa_fields[3]
                # CIGAR'dan M uzunluğunu al
                import re
                match = re.findall(r'(\d+)M', sa_cigar)
                sa_end = sa_start + sum(map(int, match)) if match else sa_start
                reads[read.query_name].append((sa_start, sa_end))

    circular_reads = []
    for qname, aligns in reads.items():
        read_len = max(end for start, end in aligns) - min(start for start, end in aligns)

        # 1. read uzunluğu referans ile uyumlu mu?
        if not (ref_len * (1 - tolerance) <= read_len <= ref_len * (1 + tolerance)):
            continue

        # 2. alignment parçaları ≥2 mi?
        if len(aligns) < 2:
            continue

        # 3. referansın başı ve sonu kapsanıyor mu? (parçaların birleşik kapsaması)
        min_start = min(start for start, end in aligns)
        max_end = max(end for start, end in aligns)
        if not (min_start < 100 and max_end > ref_len - 100):
            continue

        circular_reads.append((qname, read_len))

    # FASTQ olarak kaydet
    fastq_out = os.path.join(outdir, "circular_reads.fastq")
    with open(fastq_out, "w") as out_f:
        for qname, _ in circular_reads:
            if qname in fastq_dict:
                SeqIO.write(fastq_dict[qname], out_f, "fastq")

    print(f"✔ Circular read FASTQ kaydedildi: {fastq_out}")
    return circular_reads



def plot_linear(bam_file, ref_name, read_name, ref_len, outdir):
    bam = pysam.AlignmentFile(bam_file, "rb")
    aligns = [aln for aln in bam.fetch(ref_name) if aln.query_name == read_name]

    plt.figure(figsize=(12, 2))
    plt.hlines(1, 0, ref_len, colors="black", linestyles="--", label="Reference")

    for aln in aligns:
        color = "blue" if not aln.is_reverse else "red"
        # SA tag varsa split
        segs = [(aln.reference_start, aln.reference_end)]
        if aln.has_tag("SA"):
            sa_tag = aln.get_tag("SA")
            for sa in sa_tag.split(";"):
                if sa:
                    parts = sa.split(",")
                    sa_start = int(parts[1])
                    sa_strand = parts[2]
                    sa_cigar = parts[3]
                    # Cigar’dan soft/hard clipping çıkar
                    match_len = sum([int(n) for n in ''.join([c if c.isdigit() else ' ' for c in sa_cigar]).split()])
                    sa_end = sa_start + match_len
                    segs.append((sa_start, sa_end))
        for start, end in segs:
            plt.hlines(1, start, end, colors=color, linewidth=5)

    plt.title(f"Linear alignment of {read_name}")
    plt.xlabel("Reference position (bp)")
    plt.yticks([])
    plt.tight_layout()

    out_png = os.path.join(outdir, f"{read_name}_linear.png")
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"✔ Linear plot saved: {out_png}")


# === Circular plot ===
def plot_circular(bam_file, ref_name, read_name, ref_len, outdir):
    bam = pysam.AlignmentFile(bam_file, "rb")
    aligns = [aln for aln in bam.fetch(ref_name) if aln.query_name == read_name]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': 'polar'})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2.0)
    ax.set_yticks([])
    ax.set_xticks([])

    # Her alignment parçasını çiz
    for aln in aligns:
        color = "blue" if not aln.is_reverse else "red"
        segs = [(aln.reference_start, aln.reference_end)]
        if aln.has_tag("SA"):
            sa_tag = aln.get_tag("SA")
            for sa in sa_tag.split(";"):
                if sa:
                    parts = sa.split(",")
                    sa_start = int(parts[1])
                    sa_cigar = parts[3]
                    match_len = sum([int(n) for n in ''.join([c if c.isdigit() else ' ' for c in sa_cigar]).split()])
                    sa_end = sa_start + match_len
                    segs.append((sa_start, sa_end))
        for start, end in segs:
            theta = np.linspace((start / ref_len) * 2 * np.pi, (end / ref_len) * 2 * np.pi, 100)
            ax.plot(theta, np.ones_like(theta), color=color, linewidth=4)

    plt.title(f"Circular alignment of {read_name}")
    out_png = os.path.join(outdir, f"{read_name}_circular.png")
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"✔ Circular plot saved: {out_png}")







# === 4. Main ===
if __name__ == "__main__":
    ref_fasta = "/media/lin-bio/back2/picota_out/blasted/SRR12917049_k105/199_Cycle_513-len3350-.fasta"
    fastq_file = "/media/lin-bio/back2/picota_test/testCycle513.fastq"

    # temp output klasörü
    outdir = "/media/lin-bio/back2/picota_IS26_test/circle_Control"
    os.makedirs(outdir, exist_ok=True)

    # 1. Mapping
    sorted_bam = run_minimap2(ref_fasta, fastq_file)

    # 2. Referans uzunluğu ve ismi
    ref_record = next(SeqIO.parse(ref_fasta, "fasta"))
    ref_len = len(ref_record)
    ref_name = ref_record.id
    print(ref_len)
    # 3. Circular reads bul
    circular_reads = find_circular_reads(sorted_bam, ref_name, ref_len, fastq_file, outdir)

    print("\nCircular candidate reads (isim, uzunluk):")
    for r in circular_reads:
        print(r)
        # 4. Görselleştir (lineer + circular)
        plot_linear(sorted_bam, ref_name, r[0], ref_len, outdir)
        plot_circular(sorted_bam, ref_name, r[0], ref_len, outdir)

    print(f"\n✔ Tüm çıktılar '{outdir}' klasörüne kaydedildi ✅")




'''
if __name__ == "__main__":

    main(
        reference="/media/lin-bio/back2/picota_out/blasted/SRR12917049_k105/199_Cycle_513-len3350-.fasta",
        reads="/media/lin-bio/back2/picota_test/SRR12917036/sra_files/SRR12917036_1.fastq",
        output="/media/lin-bio/back2/picota_out/longreadcheck/circular_reads.txt",
        fastq_out="/media/lin-bio/back2/picota_out/longreadcheck/circular_reads.fastq",
        image_out="/media/lin-bio/back2/picota_out/longreadcheck/circular_reads.png",
        min_cov=0.9
    )
'''