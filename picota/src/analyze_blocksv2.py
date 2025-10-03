import re
from collections import defaultdict
import matplotlib.pyplot as plt
import os

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

def check_adjacent(coords, tolerance=1000, small_overlap=300):
    """
    Ardışıklık kontrolü:
    - Gap tolerans içinde ise adjacent
    - Küçük overlap (ör. -100bp) varsa adjacent
    - Büyük overlap varsa adjacent sayma
    """
    for (prev_s, prev_e), (next_s, next_e) in zip(coords, coords[1:]):
        gap = next_s - prev_e

        if gap >= 0:  
            # Normal boşluk
            if gap > tolerance:
                return False
        else:  
            # Overlap durumu
            if abs(gap) > small_overlap:  
                return False

    return True


def check_overlap(prev, next_, min_overlap=1, max_tolerated=300):
    """
    İki aralık overlap ediyor mu?
    - min_overlap: en az kaç bp üst üste binerse overlap saysın
    - max_tolerated: bu kadar bp'den küçük overlapleri sayma (relax)
    """
    prev_s, prev_e = prev
    next_s, next_e = next_
    overlap_len = min(prev_e, next_e) - max(prev_s, next_s)

    if overlap_len >= min_overlap:
        return overlap_len > max_tolerated  # küçük overlapleri "yok" gibi say
    return False


def find_patterns_with_insertions(reads, read_lengths, tolerance=1000):
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
    pattern_positions = {p: [] for p in patterns}
    insertions = defaultdict(list)

    for read_id, blocks in reads.items():
        srr_id = read_id.split(".")[0]
        n = len(blocks)
        readlen = read_lengths.get(read_id, None)
        if not readlen or readlen <= 0:
            continue

        # --- Transposon-in-cargo tespiti ---
        overlap_indices = set()
        for i, (b_type, b_start, b_end) in enumerate(blocks):
            if b_type == "transposon":
                # Önceki cargo ile kontrol
                if i > 0 and blocks[i-1][0] == "cargo" and check_overlap((blocks[i-1][1], blocks[i-1][2]), (b_start, b_end)):
                    insertions[read_id].append(("cargo-inserted-transposon", i-1, i))
                    overlap_indices.update([i-1, i])
                # Sonraki cargo ile kontrol
                if i < n-1 and blocks[i+1][0] == "cargo" and check_overlap((b_start, b_end), (blocks[i+1][1], blocks[i+1][2])):
                    insertions[read_id].append(("transposon-in-cargo", i, i+1))
                    overlap_indices.update([i, i+1])

        # --- Pattern tespiti (overlap içeren blokları sayma) ---
        for pname, ptypes in patterns.items():
            k = len(ptypes)
            for i in range(n - k + 1):
                window = blocks[i:i+k]
                types = [b[0] for b in window]
                coords = [(b[1], b[2]) for b in window]

                # Eğer window içinde overlap var ise sayma
                if any((i+j) in overlap_indices for j in range(k)):
                    continue

                if types == ptypes and check_adjacent(coords, tolerance):
                    srr_counts[srr_id][pname] += 1
                    total[pname] += 1
                    pattern_readlens[pname].append(readlen)
                    pattern_start = coords[0][0]
                    pattern_end = coords[-1][1]
                    center_norm = ((pattern_start + pattern_end) / 2.0) / readlen
                    center_norm = max(0.0, min(1.0, center_norm))
                    pattern_positions[pname].append(center_norm)

    return srr_counts, total, pattern_readlens, pattern_positions, insertions

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
    # Reads ve read lengths oku
    reads, read_lengths = parse_annotated_files([file1, file2])
    
    # Pattern + transposon-in-cargo tespiti
    srr_counts, total, pattern_readlens, pattern_positions, insertions = find_patterns_with_insertions(reads, read_lengths)
    
    # Çıktı dosyası oluştur
    os.makedirs(figure_out_dir, exist_ok=True)
    output_file = os.path.join(figure_out_dir, f"srr_summary_{figure_name}.txt")
    
    with open(output_file, "w") as f:
        # SRR bazında pattern sayıları
        f.write("=== SRR ===\n")
        print("=== SRR ===")
        for srr_id, d in srr_counts.items():
            total_patterns = sum(d.values())
            line_parts = [f"{k}={v}" for k, v in d.items() if v > 0]
            
            # Eğer insertions varsa ekle
            ins_events = insertions.get(srr_id, [])
            if ins_events:
                ins_str = "; ".join([f"{e[0]}({e[1]}->{e[2]})" for e in ins_events])
                line_parts.append(f"insertions: {ins_str}")
            
            if line_parts:
                line = f"{srr_id}: " + ", ".join(line_parts)
                print(line)
                f.write(line + "\n")
        
        # Genel toplam
        f.write("\n=== Total Patterns ===\n")
        print("\n=== Total Patterns ===")
        print(total)
        f.write(str(total) + "\n")
        
        # Genel insertions sayısı
        total_insertions = sum(len(v) for v in insertions.values())
        f.write(f"\n=== Total Insertions ===\n{total_insertions}\n")
        print(f"\n=== Total Insertions ===\n{total_insertions}")
    
    # Histogram ve scatter plotlar
    plot_histograms(pattern_readlens, figure_out_dir, figure_name)
    plot_scatter_positions(pattern_readlens, pattern_positions, figure_out_dir, figure_name)
    plot_scatter_positions_logscale(pattern_readlens, pattern_positions, figure_out_dir, figure_name)
    
    print("\n[INFO] figures/ klasörüne histogram ve position scatter'ları kaydedildi.")
    print("[INFO] insertions ve pattern özetleri txt dosyasına kaydedildi.")
