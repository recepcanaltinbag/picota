import re
from collections import defaultdict
import matplotlib.pyplot as plt
import os

def parse_annotated_file(file_path):
    reads = {}
    read_lengths = {}
    current_read = None
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_read = line[1:]
                reads[current_read] = []
            else:
                m = re.search(r"\((transposon|cargo)\).*Read (\d+)-(\d+) \| ReadLen=(\d+)", line)
                if m:
                    block_type = m.group(1)
                    start, end, rlen = int(m.group(2)), int(m.group(3)), int(m.group(4))
                    reads[current_read].append((block_type, start, end))
                    read_lengths[current_read] = rlen
    return reads, read_lengths

import re

def parse_annotated_files(file_list):
    """
    Birden fazla annotated read dosyasını okur.
    - Aynı read_id varsa blokları birleştirir.
    - Read length son dosyadan alınır.
    """
    reads_all = {}
    read_lengths_all = {}

    for file_path in file_list:
        with open(file_path, "r") as f:
            current_read = None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    current_read = line[1:]
                    if current_read not in reads_all:
                        reads_all[current_read] = []
                else:
                    m = re.search(r"\((transposon|cargo)\).*Read (\d+)-(\d+) \| ReadLen=(\d+)", line)
                    if m:
                        block_type = m.group(1)
                        start, end, rlen = int(m.group(2)), int(m.group(3)), int(m.group(4))
                        reads_all[current_read].append((block_type, start, end))
                        read_lengths_all[current_read] = rlen  # son dosyadaki read length

    return reads_all, read_lengths_all


def check_adjacent(coords, tolerance=100):
    for (prev_s, prev_e), (next_s, next_e) in zip(coords, coords[1:]):
        if abs(next_s - prev_e) > tolerance:
            return False
    return True

def find_patterns(reads, read_lengths, tolerance=100):
    patterns = {
        "CT": ["cargo", "transposon"],
        "TC": ["transposon", "cargo"],
        "TCT": ["transposon", "cargo", "transposon"],
        "CTC": ["cargo", "transposon", "cargo"],
        "CTCT": ["cargo", "transposon", "cargo", "transposon"],
        "TCTC": ["transposon", "cargo", "transposon", "cargo"],
        "TCTCT": ["transposon", "cargo", "transposon", "cargo", "transposon"],
        "CTCTC": ["cargo", "transposon", "cargo", "transposon", "cargo"],
        "TCTCTC": ["transposon", "cargo", "transposon", "cargo", "transposon", "cargo"],
        "TCTCTCT": ["transposon", "cargo", "transposon", "cargo", "transposon", "cargo", "transposon"],
    }

    srr_counts = defaultdict(lambda: {p: 0 for p in patterns})
    total = {p: 0 for p in patterns}
    pattern_readlens = {p: [] for p in patterns}
    pattern_positions = {p: [] for p in patterns}  # <-- burayı ekledim

    for read_id, blocks in reads.items():
        srr_id = read_id.split(".")[0]
        n = len(blocks)
        readlen = read_lengths.get(read_id, None)
        if not readlen or readlen <= 0:
            continue

        for pname, ptypes in patterns.items():
            k = len(ptypes)
            for i in range(n - k + 1):
                window = blocks[i:i+k]
                types = [b[0] for b in window]
                coords = [(b[1], b[2]) for b in window]

                if types == ptypes and check_adjacent(coords, tolerance):
                    srr_counts[srr_id][pname] += 1
                    total[pname] += 1
                    pattern_readlens[pname].append(readlen)

                    # pattern pozisyonu: pattern'in center'ı (baş ve son blok arasında)
                    pattern_start = coords[0][0]
                    pattern_end = coords[-1][1]
                    center_norm = ((pattern_start + pattern_end) / 2.0) / readlen
                    # güvenlik: 0..1 aralığının dışına çıkarsa clamp
                    center_norm = max(0.0, min(1.0, center_norm))
                    pattern_positions[pname].append(center_norm)

    return srr_counts, total, pattern_readlens, pattern_positions


def find_patterns_by_ref(reads, read_lengths, tolerance=100):
    patterns = {
        "CT": ["cargo", "transposon"],
        "TC": ["transposon", "cargo"],
        "TCT": ["transposon", "cargo", "transposon"],
        "CTC": ["cargo", "transposon", "cargo"],
        "CTCT": ["cargo", "transposon", "cargo", "transposon"],
        "TCTC": ["transposon", "cargo", "transposon", "cargo"],
        "TCTCT": ["transposon", "cargo", "transposon", "cargo", "transposon"],
        "CTCTC": ["cargo", "transposon", "cargo", "transposon", "cargo"],
        "TCTCTC": ["transposon", "cargo", "transposon", "cargo", "transposon", "cargo"],
        "TCTCTCT": ["transposon", "cargo", "transposon", "cargo", "transposon", "cargo", "transposon"],
    }

    # srr_counts[srr_id][ref_name][pattern]
    srr_counts = defaultdict(lambda: defaultdict(lambda: {p: 0 for p in patterns}))
    total = defaultdict(lambda: {p: 0 for p in patterns})

    for read_id, blocks in reads.items():
        srr_id = read_id.split(".")[0]
        readlen = read_lengths.get(read_id, None)
        if not readlen or readlen <= 0:
            continue

        n = len(blocks)
        for pname, ptypes in patterns.items():
            k = len(ptypes)
            for i in range(n - k + 1):
                window = blocks[i:i+k]
                types = [b[0] for b in window]
                coords = [(b[1], b[2]) for b in window]

                # referans adı: blokta yoksa unknown
                ref_name = window[0][3] if len(window[0]) > 3 else "unknown"

                if types == ptypes and check_adjacent(coords, tolerance):
                    srr_counts[srr_id][ref_name][pname] += 1
                    total[ref_name][pname] += 1

    return srr_counts, total


def plot_histograms(pattern_readlens, outdir="figures", prefix="global"):
    os.makedirs(outdir, exist_ok=True)
    for pname, lens in pattern_readlens.items():
        if not lens:
            continue
        plt.figure(figsize=(6,4))
        plt.hist(lens, bins=30)
        plt.title(f"{prefix} - {pname} ReadLen distribution")
        plt.xlabel("ReadLen")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{prefix}_{pname}_hist.png"), dpi=300)
        plt.close()

def plot_scatter_positions(pattern_readlens, pattern_positions, outdir="figures", prefix="global"):
    os.makedirs(outdir, exist_ok=True)
    for pname in pattern_readlens.keys():
        lens = pattern_readlens[pname]
        pos = pattern_positions[pname]
        if not lens:
            continue
        plt.figure(figsize=(7,5))
        plt.scatter(lens, pos, alpha=0.6, s=30)
        plt.xlabel("ReadLen")
        plt.ylabel("Normalized pattern center (0=baş, 1=son)")
        plt.title(f"{prefix} - {pname} position vs ReadLen")
        plt.ylim(-0.05, 1.05)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{prefix}_{pname}_pos_vs_len.png"), dpi=300)
        plt.close()

def plot_scatter_positions_logscale(pattern_readlens, pattern_positions, outdir="figures", prefix="global"):
    """
    Scatter plot:
    X ekseni: ReadLen logaritmik ölçek
    Y ekseni: pattern center normalize edilmiş (0=baş, 1=son)
    """
    os.makedirs(outdir, exist_ok=True)
    for pname in pattern_readlens.keys():
        lens = pattern_readlens[pname]
        pos = pattern_positions[pname]

        # Boy eşleşmesi ve geçersiz değerleri filtrele
        filtered = [(l, p) for l, p in zip(lens, pos) if l and p is not None and l > 0]
        if not filtered:
            continue
        lens_filtered, pos_filtered = zip(*filtered)

        plt.figure(figsize=(7,5))
        plt.scatter(lens_filtered, pos_filtered, alpha=0.6, s=30, color="royalblue", edgecolor="k")
        plt.xscale("log")  # logaritmik X ekseni
        plt.xlabel("ReadLen (log scale)")
        plt.ylabel("Normalized pattern center (0=baş, 1=son)")
        plt.title(f"{prefix} - {pname} ReadLen vs Position (log scale)")
        plt.ylim(-0.05, 1.05)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{prefix}_{pname}_pos_vs_len_log.png"), dpi=300)
        plt.close()



def analyze_blocks(file1, file2, figure_out_dir, figure_name):

    reads, read_lengths = parse_annotated_files([file1, file2])
    srr_counts, total, pattern_readlens, pattern_positions = find_patterns(reads, read_lengths)
    output_file = os.path.join(figure_out_dir, "srr_summary.txt")

    with open(output_file, "w") as f:
        # SRR bazında
        f.write("=== SRR  ===\n")
        print("=== SRR  ===")
        for s, d in srr_counts.items():
            total_patterns = sum(d.values())
            if total_patterns > 0:
                line = f"{s}: " + ", ".join(f"{k}={v}" for k,v in d.items() if v>0)
                print(line)
                f.write(line + "\n")

        # Genel toplam
        f.write("\n=== Total ===\n")
        print("\n=== Total ===")
        print(total)
        f.write(str(total) + "\n")

    # global histogramlar (ReadLen dağılımı)
    plot_histograms(pattern_readlens,figure_out_dir, figure_name)
    # ReadLen vs normalized position scatter'ları
    plot_scatter_positions(pattern_readlens, pattern_positions, figure_out_dir, figure_name)
    plot_scatter_positions_logscale(pattern_readlens, pattern_positions, figure_out_dir, figure_name)
    print("\n[INFO] figures/ klasörüne histogram ve position scatter'ları kaydedildi.")
