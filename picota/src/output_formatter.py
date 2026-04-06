"""
output_formatter.py – Enriched CSV generator for PICOTA results.

Reads picota_final_tab (TSV) and produces picota_enriched.csv where:
  - Each CT gets a unique tag (CT001, CT002, ...)
  - Antibiotic classes are expanded to one row per class
  - IS family/group are annotated
  - Key metadata (CT length, IS length, score, organism) is captured
"""

import csv
import json
import os
import re
import urllib.request
import urllib.parse
from pathlib import Path

# ---------------------------------------------------------------------------
# IS-element family lookup — ISFinder/ISFDB superfamily taxonomy
# (Siguier et al. 2006, NAR; updated from ISFinder 2023)
#
# Strategy (applied in order):
#   1. Strip variant suffix (IS26_1 → IS26, ISEcp1B → ISEcp1B)
#   2. Exact lookup in _IS_EXACT dict (case-insensitive key)
#   3. Regex-prefix rules (longest first) for named series not in exact dict
#   4. Fallback: return IS+digits as both group and family
#
# This avoids the classic prefix-collision errors:
#   IS10 starts with IS1  → would incorrectly land in IS1 superfamily
#   IS26 starts with IS2  → would incorrectly land in IS3 superfamily
#   IS21 starts with IS2  → same trap
# ---------------------------------------------------------------------------

# Exact base-name → (IS_Group / superfamily, IS_Family) mapping.
# Keys are lowercase. IS_Group = ISFinder superfamily designation.
_IS_EXACT: dict = {
    # IS1 superfamily
    "is1":    ("IS1",          "IS1"),
    "is1a":   ("IS1",          "IS1"),
    "is1b":   ("IS1",          "IS1"),
    "issen":  ("IS1",          "ISSen"),
    "issen1": ("IS1",          "ISSen"),
    "issen2": ("IS1",          "ISSen"),
    # IS3 superfamily (IS2 is an IS3-family element, not IS1!)
    "is2":    ("IS3",          "IS2"),
    "is3":    ("IS3",          "IS3"),
    "is150":  ("IS3",          "IS150"),
    "is407":  ("IS3",          "IS407"),
    "is51":   ("IS3",          "IS51"),
    "is679":  ("IS3",          "IS679"),
    # IS4 superfamily
    "is4":    ("IS4",          "IS4"),
    "is10":   ("IS4",          "IS10"),
    "is231":  ("IS4",          "IS231"),
    "is1151": ("IS4",          "IS1151"),
    # IS5 superfamily
    "is5":    ("IS5",          "IS5"),
    "is427":  ("IS5",          "IS427"),
    "is903":  ("IS5",          "IS903"),
    "is1031": ("IS5",          "IS1031"),
    # IS6 superfamily  ← IS26 lives here, NOT IS3
    "is6":    ("IS6",          "IS6"),
    "is15":   ("IS6",          "IS15"),
    "is26":   ("IS6",          "IS26"),
    "is257":  ("IS6",          "IS257"),
    "is1936": ("IS6",          "IS1936"),
    # IS21 superfamily
    "is21":   ("IS21",         "IS21"),
    "is408":  ("IS21",         "IS408"),
    # IS30 superfamily
    "is30":   ("IS30",         "IS30"),
    "is1655": ("IS30",         "IS1655"),
    # IS66 superfamily
    "is66":   ("IS66",         "IS66"),
    "is1133": ("IS66",         "IS1133"),
    # IS91 superfamily
    "is91":   ("IS91",         "IS91"),
    "iscr":   ("IS91",         "ISCR"),
    "iscr1":  ("IS91",         "ISCR"),
    "iscr2":  ("IS91",         "ISCR"),
    "iscr3":  ("IS91",         "ISCR"),
    # IS110 superfamily
    "is110":  ("IS110",        "IS110"),
    "is492":  ("IS110",        "IS492"),
    # IS200/IS605 superfamily
    "is200":  ("IS200/IS605",  "IS200"),
    "is605":  ("IS200/IS605",  "IS605"),
    "is608":  ("IS200/IS605",  "IS608"),
    # IS256 superfamily
    "is256":  ("IS256",        "IS256"),
    "is1294": ("IS256",        "IS1294"),
    "is1301": ("IS256",        "IS1301"),
    "isapl1": ("IS30",         "ISApl1"),   # ISFinder: IS30 superfamily
    "isapl2": ("IS3",          "IS150"),    # ISFinder: IS3/IS150
    "isapl3": ("IS3",          "IS150"),
    # IS481 superfamily
    "is481":  ("IS481",        "IS481"),
    # IS630 superfamily
    "is630":  ("IS630",        "IS630"),
    "is869":  ("IS630",        "IS869"),
    # IS701 superfamily
    "is701":  ("IS701",        "IS701"),
    # IS982 superfamily
    "is982":  ("IS982",        "IS982"),
    # IS1111 superfamily
    "is1111": ("IS1111",       "IS1111"),
    # IS1380 superfamily  ← ISEcp1 lives here
    "is1380": ("IS1380",       "IS1380"),
    "isecp":  ("IS1380",       "ISEcp"),
    "isecp1": ("IS1380",       "ISEcp"),
    "isecp1b":("IS1380",       "ISEcp"),
    # IS1595 superfamily
    "is1595": ("IS1595",       "IS1595"),
    "isncy":  ("IS1595",       "ISNcy"),
    # ISAs1 superfamily
    "isas1":  ("ISAs1",        "ISAs1"),
    "is3000": ("ISAs1",        "ISAs1"),
    # Composite / Tn
    "tn":     ("Composite",    "Composite"),
}

# Regex prefix rules for IS *series* variants not individually listed above.
# Checked after exact lookup fails. Sorted longest-pattern first (greedy).
# Format: (compiled_re, IS_Group/superfamily, IS_Family)
_IS_REGEX_RULES = sorted([
    (re.compile(r'^IS26\b',   re.I), "IS6",          "IS26"),
    (re.compile(r'^IS10\b',   re.I), "IS4",          "IS10"),
    (re.compile(r'^IS21\b',   re.I), "IS21",         "IS21"),
    (re.compile(r'^IS30\b',   re.I), "IS30",         "IS30"),
    (re.compile(r'^IS66\b',   re.I), "IS66",         "IS66"),
    (re.compile(r'^IS91\b',   re.I), "IS91",         "IS91"),
    (re.compile(r'^ISCR\b',   re.I), "IS91",         "ISCR"),
    (re.compile(r'^IS110\b',  re.I), "IS110",        "IS110"),
    (re.compile(r'^IS256\b',  re.I), "IS256",        "IS256"),
    (re.compile(r'^IS481\b',  re.I), "IS481",        "IS481"),
    (re.compile(r'^IS630\b',  re.I), "IS630",        "IS630"),
    (re.compile(r'^IS701\b',  re.I), "IS701",        "IS701"),
    (re.compile(r'^IS982\b',  re.I), "IS982",        "IS982"),
    (re.compile(r'^IS1111\b', re.I), "IS1111",       "IS1111"),
    (re.compile(r'^IS1380\b', re.I), "IS1380",       "IS1380"),
    (re.compile(r'^IS1595\b', re.I), "IS1595",       "IS1595"),
    (re.compile(r'^ISEcp\b',  re.I), "IS1380",       "ISEcp"),
    (re.compile(r'^ISApl\b',  re.I), "IS256",        "ISApl"),
    (re.compile(r'^ISAs1\b',  re.I), "ISAs1",        "ISAs1"),
    (re.compile(r'^IS6\b',    re.I), "IS6",          "IS6"),
    (re.compile(r'^IS5\b',    re.I), "IS5",          "IS5"),
    (re.compile(r'^IS4\b',    re.I), "IS4",          "IS4"),
    (re.compile(r'^IS3\b',    re.I), "IS3",          "IS3"),
    (re.compile(r'^IS2\b',    re.I), "IS3",          "IS2"),    # IS2 ∈ IS3 superfamily
    (re.compile(r'^IS1\b',    re.I), "IS1",          "IS1"),
], key=lambda x: len(x[0].pattern), reverse=True)

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


# ---------------------------------------------------------------------------
# ISFinder taxonomy table — loaded lazily from isfinder_taxonomy.tsv
# Generated by scripts/build_isfinder_taxonomy.py from the bundled IS.fna.
# Covers ~5570 IS elements with ISFinder-sourced superfamily annotations.
# ---------------------------------------------------------------------------
_isfinder_table: dict = {}   # lowercase IS_Name → (IS_Family, IS_Group)
_isfinder_loaded: bool = False


def _load_isfinder_table() -> dict:
    """
    Load isfinder_taxonomy.tsv into _isfinder_table (once per process).

    TSV columns (from build_isfinder_taxonomy.py / IS.fna headers):
      IS_Name   IS_Family(=superfamily)   IS_Group(=specific group)

    PICOTA output semantics:
      IS_Group  = superfamily  → TSV IS_Family  (e.g. "IS6"  for IS26)
      IS_Family = element name → TSV IS_Group if not "unknown", else IS_Name

    Stored tuple: (IS_Group/superfamily, IS_Family/element)
    """
    global _isfinder_table, _isfinder_loaded
    if _isfinder_loaded:
        return _isfinder_table

    tsv_path = Path(__file__).resolve().parent.parent / 'DBs' / 'ISes' / 'isfinder_taxonomy.tsv'
    if not tsv_path.exists():
        _isfinder_loaded = True
        return _isfinder_table

    with open(tsv_path, encoding='utf-8', errors='replace') as fh:
        next(fh, None)  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name    = parts[0].strip()
            tsv_fam = parts[1].strip()   # superfamily (IS3, IS6 …)
            tsv_grp = parts[2].strip()   # specific group (IS2, IS26 …) or == superfamily

            if not name:
                continue

            # IS_Group  = superfamily (= tsv_fam)
            # IS_Family = specific element name:
            #   use tsv_grp when it's more specific than the superfamily,
            #   otherwise fall back to the element's own name.
            is_group = tsv_fam
            if tsv_grp and tsv_grp.lower() not in ('unknown', tsv_fam.lower()):
                is_family = tsv_grp
            else:
                is_family = name

            _isfinder_table[name.lower()] = (is_group, is_family)

    _isfinder_loaded = True
    return _isfinder_table


def infer_is_family(is_name: str) -> tuple:
    """
    Return (IS_Group, IS_Family) for a given IS element name.

    Lookup order — authoritative first:
      1. ISFinder taxonomy TSV (isfinder_taxonomy.tsv, ~5570 entries from IS.fna)
         — tried with full name, then base name (variant suffix stripped)
      2. Curated _IS_EXACT dict — covers canonical elements missing from IS.fna
         (IS10, ISCR1-3, IS1380 …) and common aliases
      3. Regex-prefix rules — handles series variants (IS26_1, ISCR3 …)
      4. IS+digits fallback

    Examples:
      IS26       → ("IS6",   "IS26")    # IS6 superfamily, not IS3
      IS26_1     → ("IS6",   "IS26")    # variant suffix stripped
      IS10       → ("IS4",   "IS10")    # IS4 superfamily (via _IS_EXACT)
      IS21       → ("IS21",  "IS21")    # own superfamily
      ISEcp1     → ("IS1380","IS1380")  # via ISFinder TSV
      ISCR3      → ("IS91",  "ISCR")   # via regex rule
      IS2        → ("IS3",   "IS2")    # IS2 ∈ IS3 superfamily
    """
    name = is_name.strip()
    if not name:
        return ("Unknown", "Unknown")

    # Strip database source prefixes added by PICOTA's BLAST pipeline
    # e.g. "ISFinder_IS26" → "IS26", "ISFinder_ISEcp1" → "ISEcp1"
    for _prefix in ('ISFinder_', 'ISFDB_', 'TNcentral_'):
        if name.startswith(_prefix):
            name = name[len(_prefix):]
            break

    # Strip accession suffixes: "IS26-MH257753" → "IS26", "IS1216E-KR349520.1" → "IS1216E"
    # (done later via base split, but also handle hyphen-accession pattern here)

    table = _load_isfinder_table()

    # Step 1a — full name in ISFinder TSV
    entry = table.get(name.lower())
    if entry:
        return entry  # already (IS_Group/superfamily, IS_Family/element)

    # Step 1b — base name in ISFinder TSV (strip variant suffix: _1, _B, -like)
    base = re.split(r'[_\-:]', name)[0]
    if base.lower() != name.lower():
        entry = table.get(base.lower())
        if entry:
            return entry

    # Step 2 — curated exact dict (canonical elements + aliases)
    entry2 = _IS_EXACT.get(base.lower()) or _IS_EXACT.get(name.lower())
    if entry2:
        return entry2

    # Step 3 — regex prefix rules (longest pattern first)
    for pattern, group, family in _IS_REGEX_RULES:
        if pattern.match(base):
            return (group, family)

    # Step 4 — IS+digits fallback
    m = re.match(r'(IS\d+)', name, re.I)
    if m:
        num = m.group(1).upper()
        return (num, num)

    # Step 5 — Tn (composite transposon, not an IS element per se)
    if re.match(r'^Tn\d', name, re.I):
        return ("Composite", "Composite")

    return ("Unknown", "Unknown")


def infer_antibiotic_class(gene_or_product: str) -> str:
    """Return the antibiotic class for a gene/product name, or 'Other'."""
    text = gene_or_product.strip()
    for pattern, cls in _AMR_CLASS_PATTERNS:
        if pattern.search(text):
            return cls
    return "Other"


_aro_data = None
_MODULE_DIR = Path(__file__).resolve().parent.parent  # picota/picota/


def load_aro_data():
    global _aro_data
    if _aro_data is not None:
        return _aro_data

    aro_index_path = _MODULE_DIR / 'DBs' / 'Antibiotics' / 'card-data' / 'aro_index.tsv'
    if not aro_index_path.exists():
        _aro_data = {}
        return _aro_data

    aro_data = {}
    with open(aro_index_path, 'r', encoding='utf-8', errors='replace') as fh:
        header = fh.readline().strip().split('\t')
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            row = dict(zip(header, parts))
            aro_id = row.get('ARO Accession', '').strip()
            if aro_id:
                aro_data[aro_id.upper()] = {
                    'ARO_Name': row.get('ARO Name', '').strip(),
                    'Drug_Class': row.get('Drug Class', '').strip(),
                    'Model_Name': row.get('Model Name', '').strip()
                }
    _aro_data = aro_data
    return aro_data


def get_aro_metadata(resistance_gene: str):
    """Return aro_id, aro_name, drug_class for a Resistance_Gene value."""
    if 'ARO::' in resistance_gene:
        aro_id = resistance_gene.split('ARO::')[-1].strip()
        if aro_id and not aro_id.upper().startswith('ARO:'):
            aro_id = 'ARO:' + aro_id
    elif resistance_gene.upper().startswith('ARO:'):
        aro_id = resistance_gene.strip().upper()
    else:
        aro_id = None

    if not aro_id:
        return None, None, None

    aro_data = load_aro_data()
    entry = aro_data.get(aro_id.upper())
    if entry:
        return aro_id.upper(), entry.get('ARO_Name', ''), entry.get('Drug_Class', '')
    return aro_id.upper(), None, None


def query_sra_organism(sra_id: str):
    cache_path = _MODULE_DIR / 'DBs' / 'cache' / 'sra_taxonomy_cache.json'
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = {}
    if cache_path.exists():
        try:
            with open(cache_path, 'r', encoding='utf-8') as fh:
                cache = json.load(fh)
        except Exception:
            cache = {}

    if not sra_id:
        return 'Unknown'

    if sra_id in cache:
        return cache[sra_id]

    try:
        # Step 1: resolve accession → internal numeric UID via esearch
        esearch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
        esearch_params = {'db': 'sra', 'term': sra_id, 'retmode': 'json'}
        esearch_query = f"{esearch_url}?{urllib.parse.urlencode(esearch_params)}"
        with urllib.request.urlopen(esearch_query, timeout=20) as resp:
            esearch_data = json.load(resp)
        uid_list = esearch_data.get('esearchresult', {}).get('idlist', [])
        if not uid_list:
            raise ValueError("No UID found for accession")
        uid = uid_list[0]

        # Step 2: esummary with numeric UID → expxml contains ScientificName
        esummary_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
        esummary_params = {'db': 'sra', 'id': uid, 'retmode': 'json'}
        esummary_query = f"{esummary_url}?{urllib.parse.urlencode(esummary_params)}"
        with urllib.request.urlopen(esummary_query, timeout=20) as resp:
            esummary_data = json.load(resp)
        doc = esummary_data.get('result', {}).get(uid, {})
        expxml = doc.get('expxml', '')

        # Organism is stored as <Organism taxid="..." ScientificName="..."/>
        m = re.search(r'<Organism[^>]+ScientificName="([^"]+)"', expxml)
        if m:
            organism = m.group(1).strip()
            cache[sra_id] = organism
            with open(cache_path, 'w', encoding='utf-8') as fh:
                json.dump(cache, fh, indent=2)
            return organism
    except Exception:
        pass

    # Second fallback: pysradb ile direkt SRA metadata çek
    try:
        from pysradb import SRAweb
        db = SRAweb()
        df = db.sra_metadata(sra_id)
        if not df.empty and 'organism_name' in df.columns:
            organism = str(df.loc[0, 'organism_name']).strip()
            if organism:
                cache[sra_id] = organism
                with open(cache_path, 'w', encoding='utf-8') as fh:
                    json.dump(cache, fh, indent=2)
                return organism
    except Exception:
        # pysradb olabilir yüklü değil veya internet/dönüş yok
        pass

    cache[sra_id] = 'Unknown'
    with open(cache_path, 'w', encoding='utf-8') as fh:
        json.dump(cache, fh, indent=2)

    return 'Unknown'


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
    'CT_Tag', 'Category', 'CycleID', 'SRA_ID', 'SRA_Organism',
    'CT_Length_bp', 'Score', 'Score2', 'NumIS', 'IS_Group', 'IS_Family', 'IS_Length_bp',
    'IS_Names', 'NumAMR', 'Antibiotic_Class', 'Resistance_Gene',
    'ARO_ID', 'ARO_Name', 'ARO_Drug_Class',
    'NumXeno', 'Xenobiotic_Functions',
    'NumCompTN', 'Known_CompTN',
]

# Clean summary columns matching the human-readable view (image format)
SUMMARY_COLS = [
    'Category', 'CT_Tag', 'IS_Group', 'IS_Family',
    'Antibiotic_Class', 'Resistance_Gene',
    'CT_Length_bp', 'IS_Length_bp', 'Score2',
    'SRA_Organism', 'SRA_ID', 'Known_CompTN',
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
            score2_str = record.get('score2', '0')
            try:
                score2 = round(float(score2_str), 1)
            except ValueError:
                score2 = 0.0

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
                'Score2':             score2,
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
                row['Resistance_Gene'] = gene

                aro_id, aro_name, aro_drug_class = get_aro_metadata(gene)
                if aro_id:
                    row['ARO_ID'] = aro_id
                    row['ARO_Name'] = aro_name or ''
                    row['ARO_Drug_Class'] = aro_drug_class or cls
                    # if the ARO class exists, prefer it
                    if aro_drug_class:
                        row['Antibiotic_Class'] = aro_drug_class
                else:
                    row['ARO_ID'] = ''
                    row['ARO_Name'] = ''
                    row['ARO_Drug_Class'] = ''

                if sra_id:
                    row['SRA_Organism'] = query_sra_organism(sra_id)
                else:
                    row['SRA_Organism'] = 'Unknown'

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


def write_summary_csv(rows: list, out_csv_path: str) -> int:
    """
    Write a clean human-readable summary CSV from enriched rows.

    Columns match the spreadsheet view (Category, CT_Tag, IS_Group, IS_Family,
    Antibiotic_Class, Resistance_Gene, CT_Length_bp, IS_Length_bp, Score2,
    SRA_Organism, SRA_ID, Known_CompTN).

    Parameters
    ----------
    rows : list[dict]
        Enriched rows (from build_enriched_rows or combined across samples).
    out_csv_path : str
        Destination path.

    Returns
    -------
    int  Number of rows written.
    """
    if not rows:
        return 0

    with open(out_csv_path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

    return len(rows)
