"""
ct_cluster.py
-------------
Clusters novel CT sequences across multiple samples using blastn all-vs-all.
Produces a CTCluster tag (CTCluster001, CTCluster002, ...) for each group of
similar CTs, and writes a cross-sample summary TSV.

Algorithm:
  1. Collect all novel CT cycle FASTA sequences from annot/ directories.
  2. Build a blastn database from the combined FASTA.
  3. Run blastn all-vs-all (query = db) with identity/coverage thresholds.
  4. Build a union-find graph from hits → connected components = clusters.
  5. Assign CTCluster tags (largest cluster first).
  6. Write summary TSV: one row per (CTCluster, SRA_ID, CT_Tag).

Requires blastn and makeblastdb in PATH.
"""

import csv
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Union-Find for clustering
# ---------------------------------------------------------------------------
class _UF:
    def __init__(self):
        self._parent = {}

    def find(self, x):
        self._parent.setdefault(x, x)
        if self._parent[x] != x:
            self._parent[x] = self.find(self._parent[x])
        return self._parent[x]

    def union(self, a, b):
        self._parent[self.find(a)] = self.find(b)

    def groups(self):
        g = defaultdict(list)
        for x in self._parent:
            g[self.find(x)].append(x)
        return list(g.values())


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------
def _read_fasta(path: str) -> dict:
    """Return {header: sequence} from a FASTA file."""
    seqs = {}
    header = None
    buf = []
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    seqs[header] = ''.join(buf)
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if header is not None:
        seqs[header] = ''.join(buf)
    return seqs


def _write_fasta(records: list, path: str):
    """Write list of (header, seq) to FASTA."""
    with open(path, 'w') as fh:
        for hdr, seq in records:
            fh.write(f'>{hdr}\n{seq}\n')


# ---------------------------------------------------------------------------
# Main clustering function
# ---------------------------------------------------------------------------
def cluster_novel_cts(
    output_dir: str,
    sample_ids: list,
    enriched_rows: list,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
    logger=None,
) -> str:
    """
    Cluster novel CT sequences from all samples and write a summary TSV.

    Parameters
    ----------
    output_dir : str
        Pipeline output root (contains per-sample subdirs).
    sample_ids : list[str]
        List of short SRA IDs processed in this run.
    enriched_rows : list[dict]
        Combined enriched CSV rows from all samples.
    min_identity : float
        Minimum blastn %identity to consider two CTs similar (default 90).
    min_coverage : float
        Minimum query coverage % (default 80).
    logger :
        Optional logger instance.

    Returns
    -------
    str
        Path to the written cluster summary TSV, or '' on failure.
    """
    def _log(msg):
        if logger:
            logger.info(msg)
        else:
            print(msg)

    if not shutil.which('blastn') or not shutil.which('makeblastdb'):
        _log('  ⚠  blastn/makeblastdb not found — skipping CT clustering')
        return ''

    output_path = Path(output_dir)

    # ── Step 1: collect novel CT sequences ───────────────────────────────────
    # novel CT tags from enriched rows
    novel_tags = {
        r['CT_Tag'] for r in enriched_rows
        if r.get('Category', '').strip().lower() == 'novel'
    }
    if not novel_tags:
        _log('  ⤼ No novel CTs found — skipping clustering')
        return ''

    # Map CT_Tag → (SRA_ID, cycle_id)
    tag_meta = {}
    for r in enriched_rows:
        tag = r.get('CT_Tag', '')
        if tag not in tag_meta:
            tag_meta[tag] = {
                'SRA_ID':   r.get('SRA_ID', ''),
                'CycleID':  r.get('CycleID', ''),
                'CT_Length_bp': r.get('CT_Length_bp', ''),
                'Score':    r.get('Score', ''),
                'IS_Group': r.get('IS_Group', ''),
                'Antibiotic_Class': r.get('Antibiotic_Class', ''),
                'Category': r.get('Category', ''),
            }

    # Collect sequences from per-sample cycle FASTA files
    # File pattern: output_dir/{sra_id}/{sra_id}_cycles.fasta
    ct_sequences = {}   # CT_Tag → (header_in_fasta, sequence)
    for sra_id in sample_ids:
        cycle_fasta = output_path / sra_id / f'{sra_id}_cycles.fasta'
        if not cycle_fasta.exists():
            continue
        seqs = _read_fasta(str(cycle_fasta))
        for hdr, seq in seqs.items():
            # Match header to CT_Tag via CycleID
            for tag, meta in tag_meta.items():
                if meta['SRA_ID'] == sra_id and meta['CycleID'] and meta['CycleID'] in hdr:
                    if tag in novel_tags:
                        ct_sequences[tag] = (hdr, seq)
                    break

    if len(ct_sequences) < 2:
        _log(f'  ⤼ Only {len(ct_sequences)} novel CT sequence(s) — skipping clustering')
        return ''

    _log(f'  Clustering {len(ct_sequences)} novel CT sequences '
         f'(identity≥{min_identity}%, coverage≥{min_coverage}%)')

    # ── Step 2: blastn all-vs-all ─────────────────────────────────────────────
    with tempfile.TemporaryDirectory(dir=str(output_path)) as tmp:
        combined_fa = os.path.join(tmp, 'novel_cts.fasta')
        db_prefix   = os.path.join(tmp, 'novel_cts_db')
        blast_out   = os.path.join(tmp, 'blast_results.tsv')

        records = [(tag, seq) for tag, (_, seq) in ct_sequences.items()]
        _write_fasta(records, combined_fa)

        # Build db
        subprocess.run(
            ['makeblastdb', '-in', combined_fa, '-dbtype', 'nucl', '-out', db_prefix],
            check=True, capture_output=True
        )

        # Run blastn: outfmt 6 = qseqid sseqid pident length qlen slen
        subprocess.run(
            ['blastn', '-query', combined_fa, '-db', db_prefix,
             '-outfmt', '6 qseqid sseqid pident length qlen slen',
             '-out', blast_out,
             '-perc_identity', str(min_identity),
             '-evalue', '1e-10',
             '-num_threads', '4'],
            check=True, capture_output=True
        )

        # ── Step 3: parse hits → union-find ──────────────────────────────────
        uf = _UF()
        for tag in ct_sequences:
            uf.find(tag)   # register all nodes

        with open(blast_out) as fh:
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                q, s, pident, length, qlen, slen = parts
                if q == s:
                    continue
                try:
                    pident_f  = float(pident)
                    length_i  = int(length)
                    qlen_i    = int(qlen)
                    slen_i    = int(slen)
                except ValueError:
                    continue
                # coverage = aligned length / min(qlen, slen)
                cov = length_i / min(qlen_i, slen_i) * 100
                if pident_f >= min_identity and cov >= min_coverage:
                    uf.union(q, s)

    # ── Step 4: assign CTCluster tags ────────────────────────────────────────
    groups = uf.groups()
    # Sort groups: largest first, then by earliest CT_Tag
    groups.sort(key=lambda g: (-len(g), min(g)))

    tag_to_cluster = {}
    for idx, group in enumerate(groups, 1):
        cluster_id = f'CTCluster{idx:03d}'
        for tag in group:
            tag_to_cluster[tag] = cluster_id

    # Assign 'Unique' to singletons (no similar CT in other samples)
    for tag, cluster_id in tag_to_cluster.items():
        members = [t for t, c in tag_to_cluster.items() if c == cluster_id]
        sra_ids  = {tag_meta[t]['SRA_ID'] for t in members if t in tag_meta}
        if len(sra_ids) == 1 and len(members) == 1:
            tag_to_cluster[tag] = 'Unique'

    # ── Step 5: write summary TSV ─────────────────────────────────────────────
    summary_path = output_path / 'picota_ct_clusters.tsv'
    cols = ['CTCluster', 'CT_Tag', 'SRA_ID', 'CycleID',
            'CT_Length_bp', 'Score', 'IS_Group', 'Antibiotic_Class']

    # Count how many distinct SRA IDs share each cluster
    cluster_sra_counts = defaultdict(set)
    for tag, cluster_id in tag_to_cluster.items():
        meta = tag_meta.get(tag, {})
        cluster_sra_counts[cluster_id].add(meta.get('SRA_ID', ''))

    rows_out = []
    for tag, cluster_id in sorted(tag_to_cluster.items(),
                                   key=lambda x: (x[1], x[0])):
        meta = tag_meta.get(tag, {})
        n_samples = len(cluster_sra_counts[cluster_id])
        rows_out.append({
            'CTCluster':       cluster_id,
            'N_Samples':       n_samples,
            'CT_Tag':          tag,
            'SRA_ID':          meta.get('SRA_ID', ''),
            'CycleID':         meta.get('CycleID', ''),
            'CT_Length_bp':    meta.get('CT_Length_bp', ''),
            'Score':           meta.get('Score', ''),
            'IS_Group':        meta.get('IS_Group', ''),
            'Antibiotic_Class': meta.get('Antibiotic_Class', ''),
        })

    with open(summary_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['CTCluster', 'N_Samples'] + cols[1:],
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(rows_out)

    n_clusters = len({r['CTCluster'] for r in rows_out if r['CTCluster'] != 'Unique'})
    n_unique   = sum(1 for r in rows_out if r['CTCluster'] == 'Unique')
    _log(f'  ✓ {n_clusters} cluster(s), {n_unique} unique novel CT(s)')
    _log(f'  Summary → {summary_path}')
    return str(summary_path)
