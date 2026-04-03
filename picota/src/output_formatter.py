"""
output_formatter.py – Enriched CSV generator for PICOTA results.

Reads picota_final_tab (TSV) and produces picota_enriched.csv where:
  - Each CT gets a unique tag (CT001, CT002, ...)
  - Antibiotic classes are expanded to one row per class
  - IS family/group are annotated
  - Key metadata (CT length, IS length, score, organism) is captured
"""

import csv
import os
import re

# ---------------------------------------------------------------------------
# IS-element family lookup
# Keys are case-insensitive prefixes; values are (IS_Group, IS_Family).
# Based on ISFinder/ISFDB taxonomy.
# ---------------------------------------------------------------------------
_IS_PREFIX_MAP = [
    ("IS1 ",    ("IS1",    "IS1")),
    ("IS1$",    ("IS1",    "IS1")),
    ("IS2",     ("IS3",    "IS3")),
    ("IS3",     ("IS3",    "IS3")),
    ("IS4",     ("IS4",    "IS4")),
    ("IS5",     ("IS5",    "IS5")),
    ("IS6",     ("IS6",    "IS6")),
    ("IS10",    ("IS4",    "IS10")),
    ("IS21",    ("IS21",   "IS21")),
    ("IS26",    ("IS6",    "IS26")),
    ("IS30",    ("IS30",   "IS30")),
    ("IS66",    ("IS66",   "IS66")),
    ("IS91",    ("IS91",   "IS91")),
    ("IS110",   ("IS110",  "IS110")),
    ("IS200",   ("IS200",  "IS200/IS605")),
    ("IS256",   ("IS256",  "IS256")),
    ("IS481",   ("IS481",  "IS481")),
    ("IS630",   ("IS630",  "IS630")),
    ("IS701",   ("IS701",  "IS701")),
    ("IS982",   ("IS982",  "IS982")),
    ("IS1111",  ("IS1111", "IS1111")),
    ("IS1380",  ("IS1380", "IS1380")),
    ("IS1595",  ("IS1595", "IS1595")),
    ("IS3000",  ("ISAs1",  "ISAs1")),
    ("ISEcp",   ("IS1380", "ISEcp")),
    ("ISCR",    ("IS91",   "ISCR")),
    ("ISApl",   ("IS256",  "ISApl")),
    ("ISSen",   ("IS1",    "ISSen")),
    ("Tn",      ("Composite", "Composite")),
]

# Compiled patterns for antibiotic class inference from gene/product names
_AMR_CLASS_PATTERNS = [
    (re.compile(r'\bbla\w*', re.I),             "Beta-lactam"),
    (re.compile(r'\baac[(\[]', re.I),            "Aminoglycoside"),
    (re.compile(r'\baph[(\[]', re.I),            "Aminoglycoside"),
    (re.compile(r'\bant[(\[]', re.I),            "Aminoglycoside"),
    (re.compile(r'\baad\w*', re.I),              "Aminoglycoside"),
    (re.compile(r'\barm\w*', re.I),              "Aminoglycoside"),
    (re.compile(r'\brmtA|rmtB|rmtC|rmtD|rmtE|rmtF|rmtG|rmtH\b', re.I), "Aminoglycoside"),
    (re.compile(r'\btet[A-Z]?\b', re.I),         "Tetracycline"),
    (re.compile(r'\berm\w*', re.I),              "Macrolide-Lincosamide-Streptogramin"),
    (re.compile(r'\bmph\w*', re.I),              "Macrolide"),
    (re.compile(r'\bmsr\w*', re.I),              "Macrolide"),
    (re.compile(r'\bcfr\w*', re.I),              "Phenicol-Oxazolidinone"),
    (re.compile(r'\bcat\b|\bcml\b|\bflo\b', re.I), "Phenicol"),
    (re.compile(r'\bsul[0-9]?\b', re.I),         "Sulfonamide"),
    (re.compile(r'\bdfr\w*', re.I),              "Trimethoprim"),
    (re.compile(r'\bqnr[A-Z]?\b', re.I),         "Fluoroquinolone"),
    (re.compile(r'\bacrA|acrB|mexA|mexB|oqxA|oqxB\b', re.I), "Fluoroquinolone"),
    (re.compile(r'\bvan[A-Z]\b', re.I),           "Glycopeptide"),
    (re.compile(r'\bcol\w*|mcr\w*', re.I),        "Colistin"),
    (re.compile(r'\bmerA|merB|merC\b', re.I),     "Mercury"),
    (re.compile(r'\bqac\w*|emrA|emrB\b', re.I),  "Quaternary-ammonium"),
    (re.compile(r'\blnu\w*', re.I),               "Lincosamide"),
    (re.compile(r'\bvgb\w*|vat\w*\b', re.I),      "Streptogramin"),
    (re.compile(r'\bopa\w*|oxa\w*\b', re.I),      "Beta-lactam"),
    (re.compile(r'\bcps\w*|cfi\w*\b', re.I),      "Cephalosporin"),
]


def infer_is_family(is_name: str):
    """Return (IS_Group, IS_Family) for a given IS element name."""
    name = is_name.strip()
    # Sort by prefix length descending so more-specific entries win (e.g. IS26 before IS2)
    for prefix, (group, family) in sorted(_IS_PREFIX_MAP, key=lambda x: len(x[0]), reverse=True):
        pat = prefix.rstrip('$').rstrip()
        if name.upper().startswith(pat.upper()):
            return group, family
    # Fallback: try to extract IS prefix with digits
    m = re.match(r'(IS\d+)', name, re.I)
    if m:
        return m.group(1).upper(), m.group(1).upper()
    return "Unknown", "Unknown"


def infer_antibiotic_class(gene_or_product: str) -> str:
    """Return the antibiotic class for a gene/product name, or 'Other'."""
    text = gene_or_product.strip()
    for pattern, cls in _AMR_CLASS_PATTERNS:
        if pattern.search(text):
            return cls
    return "Other"


def parse_ct_length(cycle_id: str) -> int:
    """Extract CT length from cycle ID like 'Cycle_1-len8500-comp4-'."""
    m = re.search(r'-len(\d+)-', cycle_id)
    return int(m.group(1)) if m else 0


def parse_max_is_length(is_coords_str: str) -> int:
    """Return the max IS element length from semicolon-separated 'start-end' coords."""
    max_len = 0
    for coord in is_coords_str.split(';'):
        parts = coord.strip().split('-')
        if len(parts) == 2:
            try:
                max_len = max(max_len, abs(int(parts[1]) - int(parts[0])))
            except ValueError:
                pass
    return max_len


# ---------------------------------------------------------------------------
# Column names for picota_final_tab
# ---------------------------------------------------------------------------
_TAB_COLS = [
    'CycleID', 'SRAID', 'kmer',
    'score0', 'score1', 'score2',
    'NumIS', 'ISproducts', 'IScoords',
    'NumAnt', 'Antproducts', 'Antcoords',
    'NumXeno', 'Xenoproducts', 'Xenocoords',
    'NumCompTN', 'CompTN', 'CompTNscoords',
]

# Output columns for enriched CSV
ENRICHED_COLS = [
    'CT_Tag', 'Category', 'CycleID', 'SRA_ID',
    'CT_Length_bp', 'Score', 'NumIS', 'IS_Group', 'IS_Family', 'IS_Length_bp',
    'IS_Names', 'NumAMR', 'Antibiotic_Class', 'Resistance_Gene',
    'NumXeno', 'Xenobiotic_Functions',
    'NumCompTN', 'Known_CompTN',
]


def build_enriched_rows(tab_path: str, ct_tag_offset: int = 0) -> list:
    """
    Parse *picota_final_tab* and return a list of enriched row dicts.

    Each detected CT produces one or more rows (one per antibiotic class).
    If a CT has no AMR hit, it still gets one row with Antibiotic_Class='None'.

    Parameters
    ----------
    tab_path : str
        Path to picota_final_tab.
    ct_tag_offset : int
        Start numbering CT tags at ct_tag_offset + 1 (for multi-sample merging).
    """
    rows = []
    if not os.path.exists(tab_path):
        return rows

    ct_counter = ct_tag_offset

    with open(tab_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for record in reader:
            ct_counter += 1
            ct_tag = f"CT{ct_counter:03d}"

            cycle_id   = record.get('CycleID', '')
            sra_id     = record.get('SRAID', '')
            score_str  = record.get('score0', '0')
            try:
                score = round(float(score_str), 3)
            except ValueError:
                score = 0.0

            ct_length   = parse_ct_length(cycle_id)
            num_is      = int(record.get('NumIS', 0) or 0)
            is_products = record.get('ISproducts', '')
            is_coords   = record.get('IScoords', '')
            num_ant     = int(record.get('NumAnt', 0) or 0)
            ant_products = record.get('Antproducts', '')
            num_xeno    = int(record.get('NumXeno', 0) or 0)
            xeno_products = record.get('Xenoproducts', '')
            num_tn      = int(record.get('NumCompTN', 0) or 0)
            comptn      = record.get('CompTN', 'Novel')

            category = 'Novel' if comptn.strip() == 'Novel' else 'Known'

            # IS family/group (use first IS hit for primary annotation)
            is_names_list = [n.strip() for n in is_products.split(';') if n.strip()]
            if is_names_list:
                is_group, is_family = infer_is_family(is_names_list[0])
                is_length_bp = parse_max_is_length(is_coords)
            else:
                is_group, is_family, is_length_bp = 'None', 'None', 0

            # Antibiotic classes – one row per unique class
            ant_genes_list = [g.strip() for g in ant_products.split(';') if g.strip()]
            if ant_genes_list:
                class_gene_pairs = []
                seen_classes = set()
                for gene in ant_genes_list:
                    cls = infer_antibiotic_class(gene)
                    # Keep one representative gene per class
                    if cls not in seen_classes:
                        class_gene_pairs.append((cls, gene))
                        seen_classes.add(cls)
            else:
                class_gene_pairs = [('None', 'None')]

            base = {
                'CT_Tag':             ct_tag,
                'Category':           category,
                'CycleID':            cycle_id,
                'SRA_ID':             sra_id,
                'CT_Length_bp':       ct_length,
                'Score':              score,
                'NumIS':              num_is,
                'IS_Group':           is_group,
                'IS_Family':          is_family,
                'IS_Length_bp':       is_length_bp,
                'IS_Names':           is_products,
                'NumAMR':             num_ant,
                'NumXeno':            num_xeno,
                'Xenobiotic_Functions': xeno_products,
                'NumCompTN':          num_tn,
                'Known_CompTN':       comptn,
            }

            for cls, gene in class_gene_pairs:
                row = dict(base)
                row['Antibiotic_Class'] = cls
                row['Resistance_Gene']  = gene
                rows.append(row)

    return rows


def write_enriched_csv(tab_path: str, out_csv_path: str, ct_tag_offset: int = 0) -> int:
    """
    Generate enriched CSV from *picota_final_tab*.

    Returns the number of CTs written (rows before class-expansion).
    """
    rows = build_enriched_rows(tab_path, ct_tag_offset)
    if not rows:
        return 0

    with open(out_csv_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=ENRICHED_COLS)
        writer.writeheader()
        writer.writerows(rows)

    # Count unique CT tags
    unique_cts = len({r['CT_Tag'] for r in rows})
    return unique_cts
