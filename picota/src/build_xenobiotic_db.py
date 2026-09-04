#!/usr/bin/env python3
"""Rebuild the xenobiotic-degradation database from pathway definitions.

Why rebuild rather than patch
-----------------------------
The shipped DBs/Xenobiotics/Xenobiotics_classified.fasta (82,164 sequences) was
assembled by keyword-matching NCBI protein descriptions, and measurement of the
published tables shows what that costs:

  * 464 mobile-element proteins (transposase, integrase, recombinase, resolvase)
    are in it. They were on clean_xeno.py's exclude list, but its final line
    writes `enzymatic + xenobiotic + excluded` -- the excluded set is written
    out. For PICOTA this is not noise: a composite transposon is IS elements
    flanking cargo, so counting a transposase as cargo lets a bare IS element
    look like a composite transposon.

  * They are 0.6% of the database but 29.3% of the xenobiotic annotations in the
    published tables -- a 49x enrichment, because PICOTA searches transposon
    cycles and those are exactly where transposon proteins are.

  * A further 28.3% of published annotations are antibiotic resistance genes
    (chloramphenicol acetyltransferase and friends), which belong to CARD. Scored
    from both databases, one gene raises the cargo term twice.

  * Only 5.4% of the file matches the script's own xenobiotic keyword list. The
    other 94.6% matched `enzymatic_core` -- "transporter", "regulator", "kinase",
    "oxygenase" -- which is checked first, so a real xenobiotic enzyme is filed
    under the generic term and never tested against the specific one. "oxygenase"
    alone matches 62,469 entries.

Neither the source (kegg_download.py, which despite its name queries NCBI with
those same generic terms) nor the classifier carries any functional evidence.

What this builds instead
------------------------
KEGG's own category "Xenobiotics biodegradation and metabolism" -- 22 pathways,
597 KEGG Orthology groups -- is a published, citable definition of the gene set.
Each KO carries an EC number and its pathway membership, so every sequence in
the output can say which degradation pathway it belongs to, which the current
headers cannot. Sequences come from UniProt: reviewed (Swiss-Prot) entries are
the manually curated core; unreviewed bacterial entries add breadth under the
same EC, capped so one over-studied enzyme cannot dominate.

Then the filters that failed before, applied at the end and verified:
mobile elements out, CARD overlap out, redundancy collapsed at 90% identity.

Usage
-----
    build_xenobiotic_db.py --stage kegg      # KO metadata  -> kegg_ko.tsv
    build_xenobiotic_db.py --stage seqs      # UniProt      -> raw.fasta
    build_xenobiotic_db.py --stage filter    # filters+cdhit-> Xenobiotics_kegg.fasta
    build_xenobiotic_db.py --stage report

Stages are separate because the first two are network-bound and slow; each
writes its output and the next reads it, so an interrupted run resumes.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request

# KEGG "Xenobiotics biodegradation and metabolism" (BRITE br08901), plus
# map01220, the aromatic-degradation overview that collects the ring-cleavage
# enzymes the individual pathways share.
PATHWAYS = {
    "map00362": "Benzoate degradation",
    "map00627": "Aminobenzoate degradation",
    "map00364": "Fluorobenzoate degradation",
    "map00625": "Chloroalkane and chloroalkene degradation",
    "map00361": "Chlorocyclohexane and chlorobenzene degradation",
    "map00623": "Toluene degradation",
    "map00622": "Xylene degradation",
    "map00633": "Nitrotoluene degradation",
    "map00642": "Ethylbenzene degradation",
    "map00643": "Styrene degradation",
    "map00791": "Atrazine degradation",
    "map00930": "Caprolactam degradation",
    "map00363": "Bisphenol degradation",
    "map00621": "Dioxin degradation",
    "map00626": "Naphthalene degradation",
    "map00624": "Polycyclic aromatic hydrocarbon degradation",
    "map00365": "Furfural degradation",
    "map00984": "Steroid degradation",
    "map01220": "Degradation of aromatic compounds",
}
# Drug metabolism by cytochrome P450 is mostly vertebrate liver biology. It is
# in KEGG's xenobiotics category but adds eukaryotic P450s to a database used to
# annotate bacterial cargo, so it is opt-in rather than default.
PATHWAYS_OPTIONAL = {
    "map00980": "Metabolism of xenobiotics by cytochrome P450",
    "map00982": "Drug metabolism - cytochrome P450",
    "map00983": "Drug metabolism - other enzymes",
}

MOBILE = [
    "transposase", "integrase", "recombinase", "resolvase", "transposon",
    "insertion sequence", "is element", "intein", "retron", "conjugal transfer",
    "relaxase", "mobilization protein",
]

# The database is built from EC numbers, so every entry is an enzyme by
# construction and a KO without an EC -- a transporter, a regulator -- never
# reaches stage 2. This list is the second net: UniProt occasionally carries an
# EC on an entry whose name is a regulator or a structural subunit, and one such
# entry in a cargo database is one gene that can be miscalled as cargo.
NON_ENZYMATIC = [
    "transcriptional regulator", "transcription regulator", "response regulator",
    "two-component", "sigma factor", "repressor", "activator protein",
    "abc transporter", "mfs transporter", "permease", "porin",
    "ribosomal protein", "elongation factor", "chaperone", "cell division",
    "hypothetical protein", "uncharacterized protein", "domain-containing protein",
]

# Progress that only appears when the job ends is not progress. These stages
# run for minutes against rate-limited APIs and are normally redirected to a
# log, where Python's block buffering would hold every line until exit.
try:
    sys.stdout.reconfigure(line_buffering=True)
except AttributeError:  # pragma: no cover - Python < 3.7
    pass

KEGG = "https://rest.kegg.jp"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


def get(url, tries=3, pause=0.4):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as exc:
            if attempt == tries - 1:
                print(f"    ! {url[:90]}: {exc}", file=sys.stderr)
                return ""
            time.sleep(2 ** attempt)
    return ""


# ─── stage 1: KEGG ───────────────────────────────────────────────────────────

def stage_kegg(out_tsv, pathways):
    ko_paths = {}
    for m, name in pathways.items():
        text = get(f"{KEGG}/link/ko/{m}")
        kos = [ln.split("\t")[1].replace("ko:", "")
               for ln in text.strip().splitlines() if "\t" in ln]
        for k in kos:
            ko_paths.setdefault(k, set()).add(m)
        print(f"  {m}  {name[:46]:<46} {len(kos):>4} KO")
        time.sleep(0.35)

    print(f"\n  {len(ko_paths)} unique KO; fetching EC and name for each")
    rows = []
    for i, (ko, maps) in enumerate(sorted(ko_paths.items()), 1):
        rec = get(f"{KEGG}/get/ko:{ko}")
        symbol = name = ""
        ecs = []
        for line in rec.splitlines():
            if line.startswith("SYMBOL"):
                symbol = line.split(None, 1)[1].strip()
            elif line.startswith("NAME"):
                name = line.split(None, 1)[1].strip()
                ecs = re.findall(r"EC:([0-9.\-]+)", name)
                # A KO's NAME can carry several ECs: "[EC:1.14.13.- 1.14.13.7]"
                tail = re.search(r"\[EC:([^\]]+)\]", name)
                if tail:
                    ecs = [e for e in tail.group(1).split() if e]
                    name = name[:tail.start()].strip()
        rows.append({
            "ko": ko, "symbol": symbol, "name": name,
            "ec": ";".join(ecs), "pathways": ";".join(sorted(maps)),
        })
        if i % 50 == 0:
            print(f"    {i}/{len(ko_paths)}")
        time.sleep(0.35)

    with open(out_tsv, "w") as fh:
        fh.write("ko\tsymbol\tname\tec\tpathways\n")
        for r in rows:
            fh.write(f"{r['ko']}\t{r['symbol']}\t{r['name']}\t{r['ec']}\t{r['pathways']}\n")
    n_ec = sum(1 for r in rows if r["ec"])
    print(f"\n  wrote {out_tsv}: {len(rows)} KO, {n_ec} with an EC number")


# ─── stage 1b: pathway specificity ───────────────────────────────────────────

def stage_specific(ko_tsv, out_tsv):
    """Keep only the KOs that belong to xenobiotic degradation and nothing else.

    Anchoring on KEGG's pathway list reduced the generic-enzyme problem but did
    not remove it, because a KEGG pathway legitimately includes enzymes that
    also serve core metabolism. Measured on the 36,497 sequences the EC queries
    returned, the largest contributors were urease (700), acylphosphatase (425),
    3-hydroxyacyl-CoA dehydrogenase (403) and nitrogenase (362) -- urease is in
    the atrazine pathway because it hydrolyses a breakdown product, and it is
    also an ordinary housekeeping enzyme in a great many bacteria. Counting one
    as xenobiotic cargo repeats the mistake the rebuild exists to correct.

    So specificity is decided by KEGG's own category tree (BRITE br08901): a KO
    is kept when it appears in the "Xenobiotics biodegradation and metabolism"
    category and in no other metabolic category. The global maps -- Metabolic
    pathways, Microbial metabolism in diverse environments, Degradation of
    aromatic compounds -- are overviews that almost everything joins, so they
    count as evidence of neither.
    """
    brite = get(f"{KEGG}/get/br:br08901")
    category, pathway_cat = None, {}
    for line in brite.splitlines():
        if line.startswith("B") and len(line) > 2:
            category = line[1:].strip()
        elif line.startswith("C") and category:
            parts = line[1:].strip().split(None, 1)
            if parts and parts[0].isdigit():
                pathway_cat["map" + parts[0]] = category
    print(f"  BRITE: {len(pathway_cat)} pathways in {len(set(pathway_cat.values()))} categories")

    XENO_CAT = "Xenobiotics biodegradation and metabolism"
    OVERVIEW = "Global and overview maps"

    links = get(f"{KEGG}/link/pathway/ko")
    ko_all = {}
    for line in links.strip().splitlines():
        if "\t" not in line:
            continue
        ko, path = line.split("\t")
        ko = ko.replace("ko:", "")
        path = path.replace("path:", "")
        if path.startswith("map"):
            ko_all.setdefault(ko, set()).add(path)
    print(f"  KEGG link: pathway membership for {len(ko_all)} KO")

    with open(ko_tsv) as fh:
        header = next(fh)
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]

    kept, dropped = [], []
    for r in rows:
        ko = r[0]
        cats = {pathway_cat.get(p) for p in ko_all.get(ko, set())}
        cats.discard(None)
        cats.discard(OVERVIEW)
        other = cats - {XENO_CAT}
        (kept if (XENO_CAT in cats and not other) else dropped).append((r, sorted(other)))

    with open(out_tsv, "w") as fh:
        fh.write(header)
        for r, _ in kept:
            fh.write("\t".join(r) + "\n")

    print(f"\n  kept    {len(kept):>4} KO specific to xenobiotic degradation")
    print(f"  dropped {len(dropped):>4} KO shared with other metabolic categories")
    from collections import Counter
    why = Counter(c for _, cats in dropped for c in cats)
    for c, n in why.most_common(8):
        print(f"      {n:>4}  also in: {c}")
    print(f"\n  wrote {out_tsv}")


# ─── stage 2: UniProt ────────────────────────────────────────────────────────

def uniprot_fasta(query, size=500):
    url = f"{UNIPROT}?query={urllib.parse.quote(query)}&format=fasta&size={size}"
    return get(url)


def stage_seqs(ko_tsv, out_fasta, max_unreviewed):
    with open(ko_tsv) as fh:
        next(fh)
        kos = [dict(zip(("ko", "symbol", "name", "ec", "pathways"),
                        ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]

    ec_map = {}
    for k in kos:
        for ec in filter(None, k["ec"].split(";")):
            if ec.endswith("-"):
                continue        # a partial EC would pull in a whole subclass
            ec_map.setdefault(ec, []).append(k)
    print(f"  {len(ec_map)} complete EC numbers from {len(kos)} KO")

    seen, written, stats = set(), 0, {"reviewed": 0, "unreviewed": 0}
    with open(out_fasta, "w") as out:
        for i, (ec, entries) in enumerate(sorted(ec_map.items()), 1):
            ko_ids = ";".join(sorted({e["ko"] for e in entries}))
            paths = ";".join(sorted({p for e in entries
                                     for p in e["pathways"].split(";") if p}))
            label = entries[0]["name"] or entries[0]["symbol"] or ec

            for reviewed, cap in (("true", 500), ("false", max_unreviewed)):
                if cap <= 0:
                    continue
                q = f"ec:{ec} AND reviewed:{reviewed} AND taxonomy_id:2"
                txt = uniprot_fasta(q, size=min(cap, 500))
                for block in txt.split("\n>"):
                    block = block.strip()
                    if not block:
                        continue
                    if not block.startswith(">"):
                        block = ">" + block
                    head, _, seq = block.partition("\n")
                    seq = "".join(seq.split())
                    if not seq:
                        continue
                    acc = head.split("|")[1] if "|" in head else head[1:].split()[0]
                    if acc in seen:
                        continue
                    seen.add(acc)
                    desc = head.split(" ", 1)[1] if " " in head else ""
                    org = ""
                    mo = re.search(r"OS=(.+?)\s+(OX|GN|PE|SV)=", desc)
                    if mo:
                        org = mo.group(1)
                    prot = desc.split(" OS=")[0].strip()
                    # Provenance in the header: which KO, which EC, which
                    # pathway. The old headers were raw NCBI descriptions, so a
                    # reader could not tell benzoate degradation from ubiquinone
                    # biosynthesis.
                    out.write(f">{acc}|KO:{ko_ids}|EC:{ec}|PATH:{paths}|"
                              f"{prot.replace('|', '/')}|{org}\n{seq}\n")
                    written += 1
                    stats["reviewed" if reviewed == "true" else "unreviewed"] += 1
                time.sleep(0.3)
            if i % 25 == 0:
                print(f"    {i}/{len(ec_map)} EC, {written} sequences")

    print(f"\n  wrote {out_fasta}: {written} sequences "
          f"({stats['reviewed']} reviewed, {stats['unreviewed']} unreviewed)")


# ─── stage 3: filters ────────────────────────────────────────────────────────

def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def stage_filter(raw_fasta, card_fasta, out_fasta, identity):
    records = list(read_fasta(raw_fasta))
    print(f"  input: {len(records)} sequences")

    kept = [(h, s) for h, s in records
            if not any(m in h.lower() for m in MOBILE)]
    print(f"  after mobile-element filter:  {len(kept):>7}  (-{len(records) - len(kept)})")

    before = len(kept)
    kept = [(h, s) for h, s in kept
            if not any(m in h.lower() for m in NON_ENZYMATIC)]
    print(f"  after non-enzymatic filter:   {len(kept):>7}  (-{before - len(kept)})")

    # CARD overlap. A resistance gene scored from both databases raises the
    # cargo term twice for one gene; CARD is the right home for it.
    tmp_in = out_fasta + ".tmp.faa"
    if card_fasta and os.path.exists(card_fasta):
        with open(tmp_in, "w") as fh:
            for h, s in kept:
                fh.write(f">{h}\n{s}\n")
        db = out_fasta + ".card"
        hits = out_fasta + ".card.tsv"
        try:
            subprocess.run(["diamond", "makedb", "--in", card_fasta, "-d", db],
                           check=True, capture_output=True)
            subprocess.run(["diamond", "blastp", "-q", tmp_in, "-d", db, "-o", hits,
                            "--id", "60", "--query-cover", "70", "--max-target-seqs", "1",
                            "--quiet", "--threads", "8"], check=True, capture_output=True)
            drop = {ln.split("\t")[0] for ln in open(hits)}
            before = len(kept)
            kept = [(h, s) for h, s in kept if h.split("|")[0] not in drop
                    and h.split()[0] not in drop]
            print(f"  after CARD-overlap filter: {len(kept)}  (-{before - len(kept)})")
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            print(f"  ! CARD filter skipped ({exc}); install diamond to enable")
        finally:
            for f in (tmp_in, hits, db + ".dmnd"):
                if os.path.exists(f):
                    os.remove(f)
    else:
        print("  ! CARD filter skipped (no CARD fasta given)")

    pre = out_fasta + ".pre"
    with open(pre, "w") as fh:
        for h, s in kept:
            fh.write(f">{h}\n{s}\n")

    # Redundancy. Without this one heavily sequenced enzyme family contributes
    # hundreds of near-identical entries and its hit count reflects how often it
    # was deposited, not how often it is present.
    try:
        subprocess.run(["cd-hit", "-i", pre, "-o", out_fasta, "-c", str(identity),
                        "-n", "5", "-M", "4000", "-T", "8", "-d", "0"],
                       check=True, capture_output=True)
        n_final = sum(1 for _ in read_fasta(out_fasta))
        print(f"  after cd-hit at {identity:.0%}: {n_final}  (-{len(kept) - n_final})")
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"  ! cd-hit skipped ({exc}); writing unclustered")
        os.replace(pre, out_fasta)
        n_final = len(kept)
    finally:
        for f in (pre, out_fasta + ".clstr"):
            if os.path.exists(f) and f.endswith(".pre"):
                os.remove(f)

    # Verify, rather than assume. This is precisely the check whose absence let
    # 464 mobile-element proteins ship inside the previous database.
    leaked = [h for h, _ in read_fasta(out_fasta)
              if any(m in h.lower() for m in MOBILE + NON_ENZYMATIC)]
    if leaked:
        print(f"\n  FAIL: {len(leaked)} excluded sequences survived the filters:")
        for h in leaked[:5]:
            print(f"    {h[:100]}")
        raise SystemExit(1)
    print(f"  verified: 0 mobile-element and 0 non-enzymatic sequences in {out_fasta}")


def stage_report(fasta, ko_tsv):
    from collections import Counter
    paths = Counter()
    n = 0
    for h, _ in read_fasta(fasta):
        n += 1
        mo = re.search(r"PATH:([^|]*)", h)
        if mo:
            for p in filter(None, mo.group(1).split(";")):
                paths[p] += 1
    names = dict(PATHWAYS, **PATHWAYS_OPTIONAL)
    print(f"{fasta}: {n} sequences\n")
    print(f"  {'pathway':<12} {'sequences':>9}   name")
    for p, c in paths.most_common():
        print(f"  {p:<12} {c:>9}   {names.get(p, '')}")


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True,
                    choices=["kegg", "specific", "seqs", "filter", "report"])
    ap.add_argument("--workdir", default=os.path.join(here, "xenodb"))
    ap.add_argument("--include-cyp450", action="store_true",
                    help="add the drug-metabolism P450 pathways (mostly vertebrate)")
    ap.add_argument("--max-unreviewed", type=int, default=200,
                    help="per EC; 0 keeps only Swiss-Prot")
    ap.add_argument("--identity", type=float, default=0.90,
                    help="cd-hit clustering identity")
    ap.add_argument("--card", default="/home/lin-bio/00_GithubCodes/picotaDev/picota/"
                                      "picota/DBs/Antibiotics/protein_fasta_protein_homolog_model.fasta")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    ko_tsv = os.path.join(args.workdir, "kegg_ko.tsv")
    ko_specific = os.path.join(args.workdir, "kegg_ko_specific.tsv")
    raw = os.path.join(args.workdir, "raw.fasta")
    final = os.path.join(args.workdir, "Xenobiotics_kegg.fasta")

    if args.stage == "kegg":
        paths = dict(PATHWAYS)
        if args.include_cyp450:
            paths.update(PATHWAYS_OPTIONAL)
        stage_kegg(ko_tsv, paths)
    elif args.stage == "specific":
        stage_specific(ko_tsv, ko_specific)
    elif args.stage == "seqs":
        stage_seqs(ko_specific if os.path.exists(ko_specific) else ko_tsv,
                   raw, args.max_unreviewed)
    elif args.stage == "filter":
        stage_filter(raw, args.card, final, args.identity)
    else:
        stage_report(final, ko_tsv)


if __name__ == "__main__":
    main()
