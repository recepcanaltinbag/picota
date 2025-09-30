import subprocess
import pysam
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
import math
import numpy as np
import re 

def run_minimap2(ref_fasta, fastq_file, bam_out="mapping.bam", threads=4, run_dir=None):
    """Minimap2 ile mapping yapar ve sorted BAM döner.
    Temp dosyalar ayrı klasörde tutulur. Var olanlar tekrar kullanılmaz."""
    
    if run_dir is None:
        run_dir = "tmp_mapping"
    os.makedirs(run_dir, exist_ok=True)

    sam_out = os.path.join(run_dir, os.path.basename(bam_out).replace(".bam", ".sam"))
    bam_out = os.path.join(run_dir, os.path.basename(bam_out))
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
def find_circular_reads_old(bam_file, ref_name, ref_len, fastq_file, outdir, tolerance=0.1):
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
        read_len = len(fastq_dict[qname].seq)
        #read_len = max(end for start, end in aligns) - min(start for start, end in aligns)

        # 1. read uzunluğu referans ile uyumlu mu?
        if (read_len < ref_len * (1 - tolerance)) or (read_len > ref_len * (1 + tolerance)):
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


def plot_linear(bam_file, ref_name, ref_len, outdir):
    bam = pysam.AlignmentFile(bam_file, "rb")

    plt.figure(figsize=(12, 6))  # yüksekliği artır
    plt.hlines(0, 0, ref_len, colors="black", linestyles="--", label="Reference")

    read_dict = {}  # read_id -> y_offset
    y_base = 1
    for aln in bam.fetch(ref_name):
        read_id = aln.query_name
        if read_id not in read_dict:
            read_dict[read_id] = y_base
            y_base += 1  # her yeni read farklı y konumu

        y_pos = read_dict[read_id]

        color = "blue" if not aln.is_reverse else "red"
        segs = [(aln.reference_start, aln.reference_end)]

        if aln.has_tag("SA"):
            sa_tag = aln.get_tag("SA")
            for sa in sa_tag.split(";"):
                if sa:
                    parts = sa.split(",")
                    sa_start = int(parts[1])
                    sa_cigar = parts[3]
                    match_len = sum([int(n) for n in re.findall(r'(\d+)M', sa_cigar)])
                    sa_end = sa_start + match_len
                    segs.append((sa_start, sa_end))

        for start, end in segs:
            plt.hlines(y_pos, start, end, colors=color, linewidth=3)

    plt.title(f"Linear alignment of all reads")
    plt.xlabel("Reference position (bp)")
    plt.ylabel("Reads")
    plt.yticks(list(read_dict.values()), list(read_dict.keys()))
    plt.tight_layout()

    out_png = os.path.join(outdir, f"all_reads_linear.png")
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"✔ Linear plot of all reads saved: {out_png}")



def plot_circular(bam_file, ref_name, read_name, ref_len, outdir):
    """
    Circular alignment plot. 
    - read_name: BAM'deki query_name (sadece read_id kısmı ile eşleşmeli)
    - ref_name: BAM header'daki contig adı
    """
    bam = pysam.AlignmentFile(bam_file, "rb")

    # BAM'deki query_name sadece read_id içerecek şekilde ayarla
    read_id_only = read_name.split("_")[0]
    aligns = [aln for aln in bam.fetch(ref_name) if aln.query_name == read_id_only]

    if not aligns:
        print(f"[WARN] {read_name} için alignment bulunamadı.")
        return

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': 'polar'})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2.0)
    ax.set_yticks([])
    ax.set_xticks([])

    # Her alignment parçasını çiz
    for aln_idx, aln in enumerate(aligns):
        color = "blue" if not aln.is_reverse else "red"

        # Ana alignment
        segs = [(aln.reference_start, aln.reference_end)]

        # SA tag varsa parçaları ekle
        if aln.has_tag("SA"):
            sa_tag = aln.get_tag("SA")
            for sa in sa_tag.split(";"):
                if sa:
                    parts = sa.split(",")
                    sa_start = int(parts[1])
                    sa_cigar = parts[3]
                    # Sadece M'leri hesapla
                    matches = re.findall(r'(\d+)M', sa_cigar)
                    match_len = sum(map(int, matches))
                    sa_end = sa_start + match_len
                    segs.append((sa_start, sa_end))

        # Her segmenti çiz
        for seg_idx, (start, end) in enumerate(segs):
            theta = np.linspace((start / ref_len) * 2 * np.pi,
                                (end / ref_len) * 2 * np.pi, 100)
            # Farklı segmentleri hafifçe farklı radius ile çiz
            radius = 1 + 0.1 * (aln_idx + seg_idx)
            ax.plot(theta, np.ones_like(theta) * radius, color=color, linewidth=4)

    plt.title(f"Circular alignment of {read_name}")
    out_png = os.path.join(outdir, f"{read_name}_circular.png")
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"✔ Circular plot saved: {out_png}")






# === 0. Importlar ===
import os
from Bio import SeqIO
import pysam
import re

# === 1. Fonksiyonlar ===
def load_fastq(fastq_file):
    return {rec.id: rec for rec in SeqIO.parse(fastq_file, "fastq")}

def load_reads(bam_file, ref_name):
    """
    BAM dosyasından read alignment bilgilerini yükler.
    Primary ve supplementary (SA) alignments için query koordinatlarını hesaplar.
    Returns:
        reads: dict[read_name] = list of tuples (ref_start, ref_end, query_start, query_end)
    """
    import re
    import pysam

    bam = pysam.AlignmentFile(bam_file, "rb")
    reads = {}

    for read in bam.fetch(ref_name):
        if read.query_name not in reads:
            reads[read.query_name] = []

        # --- Primary alignment ---
        reads[read.query_name].append((
            read.reference_start,
            read.reference_end,
            read.query_alignment_start,
            read.query_alignment_end
        ))

        # --- Supplementary alignments (SA tag) ---
        if read.has_tag("SA"):
            sa_tag = read.get_tag("SA")
            for sa_entry in sa_tag.split(";"):
                if not sa_entry:
                    continue
                sa_fields = sa_entry.split(",")
                sa_ref_start = int(sa_fields[1])
                sa_cigar = sa_fields[3]

                # Query alignment uzunluğunu CIGAR’dan hesapla
                q_start = 0
                q_len = 0
                matches = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
                for l, t in matches:
                    l = int(l)
                    if t in ("M", "I", "=", "X"):  # query ilerleyen
                        q_len += l
                sa_ref_end = sa_ref_start + sum(int(l) for l, t in matches if t in ("M", "D", "N", "=", "X"))

                reads[read.query_name].append((
                    sa_ref_start,
                    sa_ref_end,
                    q_start,
                    q_len
                ))

    return reads


def is_circular_ordered(aligns, ref_len, transposon_len,
                        cargo_start_win=500, cargo_end_win=500, tolerance=50):
    aligns = sorted(aligns, key=lambda x: x[2])  # query sırasına göre
    cargo_end_flag = False
    transposon_flag = False
    cargo_start_flag = False

    for ref_start, ref_end, q_start, q_end in aligns:
        # 1. Cargo sonu
        if not cargo_end_flag and ref_end > ref_len - cargo_end_win:
            cargo_end_flag = True
            continue
        # 2. Transposon (tam kapsama)
        if cargo_end_flag and not transposon_flag:
            if ref_start <= 0 + tolerance and ref_end >= transposon_len - tolerance:
                transposon_flag = True
                continue
        # 3. Cargo başı
        if cargo_end_flag and transposon_flag and not cargo_start_flag:
            if ref_start < cargo_start_win:
                cargo_start_flag = True

    return cargo_end_flag and transposon_flag and cargo_start_flag

def find_circular_reads(bam_file, ref_name, ref_len, fastq_file, transposon_len,
                        cargo_start_win=500, cargo_end_win=500, tolerance=50):
    fastq_dict = load_fastq(fastq_file)
    reads = load_reads(bam_file, ref_name)
    circular_reads = []
    for qname, aligns in reads.items():
        if is_circular_ordered(aligns, ref_len, transposon_len,
                               cargo_start_win, cargo_end_win, tolerance):
            read_len = len(fastq_dict[qname].seq)
            circular_reads.append((qname, read_len))
    return circular_reads


def debug_all_reads_aligns(bam_file, ref_name, ref_len, transposon_len,
                           cargo_start_win=1000, cargo_end_win=1000, tolerance=100):
    """
    Tüm read'leri debug için yazdırır.
    Merge edilmiş alignments ile ref üzerindeki konumları gösterir.
    """
    reads = load_reads(bam_file, ref_name)

    for qname, aligns in reads.items():
        # Merge overlapping alignments
        merged = merge_overlapping_aligns(aligns)

        print(f"\n[DEBUG] Read: {qname} (merged {len(merged)} blok)")

        cargo_end_flag = False
        transposon_flag = False
        cargo_start_flag = False

        for ref_start, ref_end, q_start, q_end in merged:
            region_hit = None
            # Cargo sonu
            if not cargo_end_flag and ref_end > ref_len - cargo_end_win:
                cargo_end_flag = True
                region_hit = "cargo_end"
            # Transposon full
            elif cargo_end_flag and not transposon_flag:
                if ref_start <= 0 + tolerance and ref_end >= transposon_len - tolerance:
                    transposon_flag = True
                    region_hit = "transposon_full"
            # Cargo başı
            elif cargo_end_flag and transposon_flag and not cargo_start_flag:
                if ref_start < cargo_start_win:
                    cargo_start_flag = True
                    region_hit = "cargo_start"

            print(f"  Align: ref({ref_start}-{ref_end}), query({q_start}-{q_end}), region_hit={region_hit}")

        print(f"  ✅ Flags: cargo_end={cargo_end_flag}, transposon_full={transposon_flag}, cargo_start={cargo_start_flag}")


def merge_overlapping_aligns(aligns):
    """
    Overlapping alignment bloklarını birleştirir (ref koordinatına göre)
    sadece query koordinatı olanlar işleme alınır
    """
    # query koordinatı olmayanları (0-0) filtrele
    aligns = [a for a in aligns if a[2] != 0 and a[3] != 0]
    if not aligns:
        return []

    # ref_start ile sırala
    aligns = sorted(aligns, key=lambda x: x[0])


    merged = []
    cur_start, cur_end, cur_qstart, cur_qend = aligns[0]

    for ref_start, ref_end, q_start, q_end in aligns[1:]:
        if ref_start <= cur_end:  # overlap varsa
            cur_end = max(cur_end, ref_end)
            cur_qend = max(cur_qend, q_end)
        else:
            merged.append((cur_start, cur_end, cur_qstart, cur_qend))
            cur_start, cur_end, cur_qstart, cur_qend = ref_start, ref_end, q_start, q_end

    merged.append((cur_start, cur_end, cur_qstart, cur_qend))
    return merged


def is_circular_ordered_debug_merged(aligns, ref_len, transposon_len,
                                     cargo_start_win=250, cargo_end_win=250, tolerance=200):
    """
    Merge edilmiş alignments ile ardışık cargo_end -> transposon_full -> cargo_start kontrolü
    Debug log verir
    """
    merged = merge_overlapping_aligns(aligns)

    cargo_end_flag = False
    transposon_flag = False
    cargo_start_flag = False
    debug_info = []

    for ref_start, ref_end, q_start, q_end in merged:
        region_hit = None
        # 1. Cargo sonu
        if not cargo_end_flag and ref_end > ref_len - cargo_end_win:
            cargo_end_flag = True
            region_hit = "cargo_end"
        # 2. Transposon (tam kapsama)
        elif cargo_end_flag and not transposon_flag:
            if ref_start <= 0 + tolerance and ref_end >= transposon_len - tolerance:
                transposon_flag = True
                region_hit = "transposon_full"
        # 3. Cargo başı
        elif cargo_end_flag and transposon_flag and not cargo_start_flag:
            if ref_start < cargo_start_win:
                cargo_start_flag = True
                region_hit = "cargo_start"

        debug_info.append((ref_start, ref_end, q_start, q_end, region_hit))

    circular = cargo_end_flag and transposon_flag and cargo_start_flag
    return circular, debug_info


def find_circular_reads_debug_merged(bam_file, ref_name, ref_len, fastq_file,
                                     transposon_len, cargo_start_win=500, cargo_end_win=500, tolerance=100):
    fastq_dict = load_fastq(fastq_file)
    reads = load_reads(bam_file, ref_name)
    circular_reads = []

    for qname, aligns in reads.items():
        circular, debug_info = is_circular_ordered_debug_merged(
            aligns, ref_len, transposon_len, cargo_start_win, cargo_end_win, tolerance
        )

        if circular:
            read_len = len(fastq_dict[qname].seq)
            circular_reads.append((qname, read_len))

            # Sadece circular read için debug
            print(f"\n[DEBUG] Circular Read: {qname}")
            for r in debug_info:
                print(f"  Align: ref({r[0]}-{r[1]}), query({r[2]}-{r[3]}), region_hit={r[4]}")
            print(f"  ✅ Circular read bulundu! uzunluk: {read_len} bp")

    return circular_reads



def debug_circular_read_regions(read_name, aligns, ref_len, transposon_len,
                                cargo_start_win=500, cargo_end_win=500, tolerance=50):
    """
    Merge edilmiş alignments ile read üzerinde blokları __region(block)__ şeklinde gösterir
    Sadece tüm flag’ler True ise yazdırır
    """
    merged = merge_overlapping_aligns(aligns)

    blocks = []
    cargo_end_flag = False
    transposon_flag = False
    cargo_start_flag = False

    for ref_start, ref_end, q_start, q_end in merged:
        region_hit = None
        # Cargo sonu
        if not cargo_end_flag and ref_end > ref_len - cargo_end_win:
            cargo_end_flag = True
            region_hit = "cargo_end"
        # Transposon (tam kapsama)
        elif cargo_end_flag and not transposon_flag:
            if ref_start <= 0 + tolerance and ref_end >= transposon_len - tolerance:
                transposon_flag = True
                region_hit = "transposon"
        # Cargo başı
        elif cargo_end_flag and transposon_flag and not cargo_start_flag:
            if ref_start < cargo_start_win:
                cargo_start_flag = True
                region_hit = "cargo_start"

        block_str = f"{region_hit}({ref_start}-{ref_end})" if region_hit else f"none({ref_start}-{ref_end})"
        blocks.append(block_str)

    # Sadece tüm flag’ler True ise debug yazdır
    if cargo_end_flag and transposon_flag and cargo_start_flag:
        print(f"\n[DEBUG] Circular Read: {read_name}")
        print("  " + "__".join(blocks))
        print(f"Flags: cargo_start={cargo_start_flag}, transposon_full={transposon_flag}, cargo_end={cargo_end_flag}")
        return True

    return False

def find_circular_readsv2(bam_file, ref_name, ref_len, fastq_file, transposon_len,
                        cargo_start_win=500, cargo_end_win=500, tolerance=50):
    fastq_dict = load_fastq(fastq_file)
    reads = load_reads(bam_file, ref_name)
    circular_reads = []

    for qname, aligns in reads.items():
        circular = debug_circular_read_regions(
            qname, aligns, ref_len, transposon_len,
            cargo_start_win, cargo_end_win, tolerance
        )
        if circular:
            read_len = len(fastq_dict[qname].seq)
            circular_reads.append((qname, read_len))
    return circular_reads


# === 4. Main ===
if __name__ == "__main__":
    base_dir = "/media/lin-bio/back2/picota_IS26_test"
    outdir = os.path.join(base_dir, "circle_Control")
    for root, _, files in os.walk(outdir):
        for f in files:
            if f.endswith(".sam"):
                os.remove(os.path.join(root,f))

    os.makedirs(outdir, exist_ok=True)

    print(f"[INFO] Çıktılar {outdir} klasörüne kaydedilecek.")

    # Referanslara göre transposon başlangıçları
    transposon_starts = {
        "Tn4352-M20306": 820,
        "Tn6010-EU370913": 820,
        "Tn6023-GU562437.2": 820,
        "Tn6309-KX710094": 820,
        "Tn6925-CP076822": 820,
        "TnPMLUA4-KC964607.1": 820

        # istediğin kadar ref ekle
    }

    fasta_files = [f for f in os.listdir(base_dir) if f.endswith(".fasta")]
    if not fasta_files:
        print("[WARN] Hiç fasta dosyası bulunamadı!")
        exit(1)

    for ref_fasta in fasta_files:
        ref_path = os.path.join(base_dir, ref_fasta)
        ref_record = next(SeqIO.parse(ref_path, "fasta"))
        ref_len = len(ref_record)
        ref_name = os.path.splitext(ref_fasta)[0]

        print(f"\n[INFO] Referans: {ref_fasta} ({ref_name}, {ref_len} bp)")

        # transposon uzunluğunu al, yoksa default 1200
        transposon_len = transposon_starts.get(ref_name, 1200)

        sra_dir = os.path.join(base_dir, ref_name, "sra_files")
        if not os.path.isdir(sra_dir):
            print(f"[WARN] {sra_dir} bulunamadı, geçiliyor.")
            continue

        fastq_files = [f for f in os.listdir(sra_dir) if f.endswith(".fastq")]
        if not fastq_files:
            print(f"[WARN] {sra_dir} içinde fastq yok, geçiliyor.")
            continue

        for fastq_file in fastq_files:
            fq_path = os.path.join(sra_dir, fastq_file)
            print(f"\n[STEP] Mapping başlatılıyor → {fastq_file}")

            run_dir = os.path.join(outdir, ref_name + "_" + os.path.splitext(fastq_file)[0])
            os.makedirs(run_dir, exist_ok=True)

            try:
                sorted_bam = run_minimap2(
                    ref_path,
                    fq_path,
                    bam_out=f"{os.path.splitext(fastq_file)[0]}_{ref_name}.bam",
                    threads=4,
                    run_dir=run_dir
                )
            except Exception as e:
                print(e)
            '''
            circular_reads = find_circular_readsv2(
                sorted_bam,
                ref_name,
                ref_len,
                fq_path,
                transposon_len=transposon_len
            )
            
            if not circular_reads:
                print("[INFO] Circular read bulunamadı.")
                continue

            print(f"[INFO] {len(circular_reads)} circular read bulundu.")

            plot_linear(sorted_bam, ref_name, ref_len, run_dir)
            for r in circular_reads:
                read_id = r[0]
                custom_id = f"{read_id}_{ref_name}"
                print(f"   ↳ {read_id} için grafik çiziliyor ({custom_id})")
                plot_circular(sorted_bam, ref_name, custom_id, ref_len, run_dir)
            '''
    print(f"\n✔ Tüm işler tamamlandı. Çıktılar: {outdir}")




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