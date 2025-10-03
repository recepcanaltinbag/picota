# picota_analysis.py
import re
import os
from collections import defaultdict, Counter
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
import logging


from src.config_loader import ToleranceConfig
from src.config_loader import Config

logger: logging.Logger = None


# -------------------------
# Parsing
# -------------------------
def parse_annotated_files(file_list: List[str]) -> Tuple[Dict[str, List[Tuple[str,int,int]]], Dict[str,int]]:
    reads_all: Dict[str, List[Tuple[str,int,int]]] = {}
    read_lengths_all: Dict[str, int] = {}

    for file_path in file_list:
        with open(file_path, "r") as fh:
            current_read = None
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    current_read = line[1:]
                    if current_read not in reads_all:
                        reads_all[current_read] = []
                else:
                    m = re.search(r"\((transposon|cargo)\).*Read\s+(\d+)-(\d+)\s*\|\s*ReadLen=(\d+)", line)
                    if m and current_read is not None:
                        block_type = m.group(1)
                        start, end, rlen = int(m.group(2)), int(m.group(3)), int(m.group(4))
                        reads_all[current_read].append((block_type, start, end))
                        read_lengths_all[current_read] = rlen

    return reads_all, read_lengths_all


# -------------------------
# Overlap utilities
# -------------------------
def check_overlap(prev: Tuple[int,int], next_: Tuple[int,int], min_overlap: int = 1, max_tolerated: int = 300) -> bool:
    prev_s, prev_e = prev
    next_s, next_e = next_
    overlap_len = min(prev_e, next_e) - max(prev_s, next_s)
    if overlap_len >= min_overlap:
        return overlap_len > max_tolerated
    return False


# -------------------------
# Gap statistics collection
# -------------------------
def build_gap_statistics(reads: Dict[str, List[Tuple[str,int,int]]]) -> Dict[Tuple[str,str], dict]:
    gap_lists = defaultdict(list)
    for blocks in reads.values():
        for prev, next_ in zip(blocks, blocks[1:]):
            prev_type, prev_s, prev_e = prev
            next_type, next_s, next_e = next_
            gap = next_s - prev_e
            gap_lists[(prev_type, next_type)].append(gap)

    gap_summary = {}
    for k, gaps in gap_lists.items():
        if gaps:
            avg_gap = sum(gaps) / len(gaps)
            mode_gap = Counter(gaps).most_common(1)[0][0]
            max_gap = max(gaps)
            gap_summary[k] = {
                'gaps': gaps,
                'count': len(gaps),
                'mode_gap': mode_gap,
                'avg_gap': avg_gap,
                'max_gap': max_gap,
                'large_count': 0
            }
        else:
            gap_summary[k] = {
                'gaps': [],
                'count': 0,
                'mode_gap': 0,
                'avg_gap': 0.0,
                'max_gap': 0,
                'large_count': 0
            }
    return gap_summary


def annotate_large_counts(gap_summary: Dict[Tuple[str,str], dict], cfg: ToleranceConfig) -> None:
    for k, info in gap_summary.items():
        mode = info['mode_gap']
        threshold = max(cfg.absolute_small_gap + 1, int(mode * 1.2))
        large_count = sum(1 for g in info['gaps'] if g > threshold)
        info['large_count'] = large_count


# -------------------------
# Adjacent / asymmetry decision
# -------------------------
def check_adjacent(coords: List[Tuple[int,int]], types: List[str], gap_summary: Dict[Tuple[str,str], dict], cfg: ToleranceConfig) -> bool:
    if len(coords) < 2:
        return True

    gaps = [next_s - prev_e for (prev_s, prev_e), (next_s, next_e) in zip(coords, coords[1:])]
    for i, gap in enumerate(gaps):
        prev_type = types[i]
        next_type = types[i+1]
        info = gap_summary.get((prev_type, next_type))
        reverse_info = gap_summary.get((next_type, prev_type))

        if gap < 0 and abs(gap) > cfg.small_overlap:
            logger.debug(f"Overlap too large: {prev_type}->{next_type}, overlap={abs(gap)}")
            return False

        if gap <= cfg.absolute_small_gap:
            logger.debug(f"Small absolute gap accepted: {prev_type}->{next_type}, gap={gap}")
            continue

        if not info:
            if gap <= cfg.tolerance:
                logger.debug(f"No stats, gap within tolerance: {prev_type}->{next_type}, gap={gap}")
                continue
            else:
                logger.debug(f"No stats and gap > tolerance: {prev_type}->{next_type}, gap={gap}")
                return False

        total_count = info.get('count', 0)
        large_count = info.get('large_count', 0)

        if abs(gap - info['mode_gap']) <= cfg.mode_gap_tolerance:
            logger.debug(f"Gap near mode accepted: {prev_type}->{next_type}, gap={gap}, mode={info['mode_gap']}")
            continue

        is_large = gap > max(cfg.absolute_small_gap, int(info['mode_gap'] * 1.2))
        if is_large and total_count > 0 and (large_count / total_count) < cfg.min_large_fraction:
            logger.debug(f"Singleton large gap ignored: {prev_type}->{next_type}, gap={gap} (large_count={large_count}, total={total_count})")
            continue

        if reverse_info:
            rev_mode = reverse_info['mode_gap']
            is_reverse_large = gap > max(cfg.absolute_small_gap, int(rev_mode * 1.2))
            if is_large != is_reverse_large:
                logger.warning(f"Asymmetric gap pattern: {prev_type}->{next_type}, gap={gap}, mode={info['mode_gap']}, reverse_mode={rev_mode}")
                continue

        final_boosted = min(cfg.max_boost, max(cfg.boosted_tolerance, int(info['mode_gap'] * 1.5)))
        if gap > final_boosted:
            logger.debug(f"Gap exceeds final_boosted: {prev_type}->{next_type}, gap={gap}, final_boosted={final_boosted} -> reject")
            return False

        logger.debug(f"Gap accepted by default: {prev_type}->{next_type}, gap={gap}")
    return True


# -------------------------
# Pattern finding and insertion detection
# -------------------------
def find_patterns_with_insertions(reads: Dict[str, List[Tuple[str,int,int]]], read_lengths: Dict[str,int], cfg: Config):
    patterns = cfg.patterns
    srr_counts = defaultdict(lambda: {p: 0 for p in patterns})
    total = {p: 0 for p in patterns}
    pattern_readlens = {p: [] for p in patterns}
    pattern_positions = {p: [] for p in patterns}
    insertions = defaultdict(list)

    gap_summary = build_gap_statistics(reads)
    annotate_large_counts(gap_summary, cfg.tolerances)

    logger.info("Gap summary (sample):")
    for k, v in gap_summary.items():
        logger.info(f"  {k}: count={v['count']}, mode={v['mode_gap']}, large_count={v['large_count']}")

    for read_id, blocks in reads.items():
        srr_id = read_id.split(".")[0]
        n = len(blocks)
        readlen = read_lengths.get(read_id)
        if not readlen or readlen <= 0:
            continue

        overlap_indices = set()
        for i, (b_type, b_start, b_end) in enumerate(blocks):
            if b_type == "transposon":
                if i > 0 and blocks[i-1][0] == "cargo":
                    if check_overlap((blocks[i-1][1], blocks[i-1][2]), (b_start, b_end)):
                        insertions[read_id].append(("cargo-inserted-transposon", i-1, i))
                        overlap_indices.update([i-1, i])
                if i < n-1 and blocks[i+1][0] == "cargo":
                    if check_overlap((b_start, b_end), (blocks[i+1][1], blocks[i+1][2])):
                        insertions[read_id].append(("transposon-in-cargo", i, i+1))
                        overlap_indices.update([i, i+1])

        for pname, ptypes in patterns.items():
            k = len(ptypes)
            for i in range(n - k + 1):
                window = blocks[i:i+k]
                types = [b[0] for b in window]
                coords = [(b[1], b[2]) for b in window]
                if any((i + j) in overlap_indices for j in range(k)):
                    continue
                if types == ptypes:
                    if check_adjacent(coords, types, gap_summary, cfg.tolerances):
                        srr_counts[srr_id][pname] += 1
                        total[pname] += 1
                        pattern_readlens[pname].append(readlen)
                        pattern_start = coords[0][0]
                        pattern_end = coords[-1][1]
                        center_norm = ((pattern_start + pattern_end) / 2.0) / readlen
                        center_norm = max(0.0, min(1.0, center_norm))
                        pattern_positions[pname].append(center_norm)

    return srr_counts, total, pattern_readlens, pattern_positions, insertions


# -------------------------
# Plotting helpers
# -------------------------
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
    os.makedirs(outdir, exist_ok=True)
    for pname in pattern_readlens.keys():
        lens = pattern_readlens[pname]
        pos = pattern_positions[pname]
        filtered = [(l, p) for l, p in zip(lens, pos) if l and p is not None and l > 0]
        if not filtered:
            continue
        lens_filtered, pos_filtered = zip(*filtered)

        plt.figure(figsize=(7,5))
        plt.scatter(lens_filtered, pos_filtered, alpha=0.6, s=30, color="royalblue", edgecolor="k")
        plt.xscale("log")
        plt.xlabel("ReadLen (log scale)")
        plt.ylabel("Normalized pattern center (0=baş, 1=son)")
        plt.title(f"{prefix} - {pname} ReadLen vs Position (log scale)")
        plt.ylim(-0.05, 1.05)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{prefix}_{pname}_pos_vs_len_log.png"), dpi=300)
        plt.close()


# -------------------------
# Main analysis entry point
# -------------------------
def analyze_blocks(file1: str, file2: str, figure_out_dir, contig_name, cfg: Config):
    global logger
    logger = logging.getLogger(cfg.logging.logger_name)
    
    reads, read_lengths = parse_annotated_files([file1, file2])
    srr_counts, total, pattern_readlens, pattern_positions, insertions = find_patterns_with_insertions(reads, read_lengths, cfg)

    os.makedirs(figure_out_dir, exist_ok=True)
    summary_path = os.path.join(figure_out_dir, f"srr_summary_{contig_name}.txt")

    with open(summary_path, "w") as fh:
        fh.write("=== SRR ===\n")
        for srr_id, d in srr_counts.items():
            line_parts = [f"{k}={v}" for k,v in d.items() if v>0]
            ins_events = insertions.get(srr_id, [])
            if ins_events:
                ins_str = "; ".join([f"{e[0]}({e[1]}->{e[2]})" for e in ins_events])
                line_parts.append(f"insertions: {ins_str}")
            if line_parts:
                fh.write(f"{srr_id}: " + ", ".join(line_parts) + "\n")

        fh.write("\n=== Total Patterns ===\n")
        fh.write(str(total) + "\n")
        logger.info("=== Total Patterns ===")
        logger.info(str(total))

        total_insertions = sum(len(v) for v in insertions.values())
        fh.write("\n=== Total Insertions ===\n")
        fh.write(str(total_insertions) + "\n")
        logger.info(f"=== Total Insertions === {total_insertions}")

    plot_histograms(pattern_readlens, figure_out_dir, contig_name)
    plot_scatter_positions(pattern_readlens, pattern_positions, figure_out_dir, contig_name)
    plot_scatter_positions_logscale(pattern_readlens, pattern_positions, figure_out_dir, contig_name)

    logger.info(f"Summary written to {summary_path}")
    logger.info(f"Figures written to {figure_out_dir}")




'''
logger = setup_logger("log.txt", level=logging.INFO)

from config_loader import load_config
cfg = load_config("picota/config.yaml")

print(cfg)

analyze_blocks(
    "/media/lin-bio/back2/picota_project_longTestIS26/mapping/SRR11108582/SRR11108582_Cycle_3-len5783-_split.fasta_partial",
    "/media/lin-bio/back2/picota_project_longTestIS26/mapping/SRR11108582/SRR11108582_Cycle_3-len5783-_split.fasta_full",
    "/media/lin-bio/back2/picota_project_longTestIS26/mapping/SRR11108582", 
    "SRR11108582_Cycle_3-len5783",
    cfg
)
'''