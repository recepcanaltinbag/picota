import os
import subprocess
from Bio import SeqIO
import pandas as pd

def make_diamond_db(fasta_path, db_path):
    if os.path.exists(db_path):
        print(f"Diamond DB already exists: {db_path}, skipping makedb.")
        return
    print(f"Creating Diamond DB for {fasta_path} ...")
    cmd = ["diamond", "makedb", "--in", fasta_path, "-d", db_path]
    subprocess.run(cmd, check=True)

def run_prodigal_with_gff(input_fasta, output_proteins, output_genes, output_gff):
    cmd = [
        "prodigal",
        "-i", input_fasta,
        "-a", output_proteins,
        "-d", output_genes,
        "-f", "gff",
        "-o", output_gff,
        "-p", "meta"
    ]
    subprocess.run(cmd, check=True)

def parse_gff_cds(gff_file):
    cds_list = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if parts[2] != "CDS":
                continue
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]
            # ID extraction (prodigal genelde ID=xxx şeklinde verir)
            id_field = ""
            for attr in attributes.split(';'):
                if attr.startswith("ID="):
                    id_field = attr[3:]
                    break
            cds_list.append({"id": id_field, "start": start, "end": end})
    return cds_list

def overlap_ratio(start1, end1, start2, end2):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_len = max(0, overlap_end - overlap_start + 1)
    length1 = end1 - start1 + 1
    length2 = end2 - start2 + 1
    return overlap_len / min(length1, length2)

def filter_overlapping_cds(cds_list, overlap_threshold=0.8):
    # cds_list sıralı start'a göre
    cds_list_sorted = sorted(cds_list, key=lambda x: x["start"])
    filtered = []
    for cds in cds_list_sorted:
        overlap_found = False
        for fcds in filtered:
            ov = overlap_ratio(cds["start"], cds["end"], fcds["start"], fcds["end"])
            if ov > overlap_threshold:
                overlap_found = True
                break
        if not overlap_found:
            filtered.append(cds)
    return filtered

def filter_protein_fasta(proteins_faa, filtered_cds, output_filtered_faa):
    # filtered_cds id'lerini set yap
    keep_ids = set(cds["id"] for cds in filtered_cds)
    records = (rec for rec in SeqIO.parse(proteins_faa, "fasta") if rec.id in keep_ids)
    SeqIO.write(records, output_filtered_faa, "fasta")


def run_diamond(query_fasta, db_path, out_file):
    cmd = [
        "diamond", "blastp",
        "--query", query_fasta,
        "--db", db_path,
        "--out", out_file,
        "--outfmt", "6",
        "--max-target-seqs", "1",
        "--quiet"
    ]
    subprocess.run(cmd, check=True)




def scoring_main(cycle_folder, picota_out_folder, path_to_antibiotics, path_to_xenobiotics, path_to_ises):
    os.makedirs(picota_out_folder, exist_ok=True)
    temp_folder = os.path.join(picota_out_folder, "temp")
    os.makedirs(temp_folder, exist_ok=True)
    
    # Diamond DB dosyalarının isimlerini fasta isimlerinden türetelim:
    path_to_antibiotics_db = path_to_antibiotics + ".dmnd"
    path_to_xenobiotics_db = path_to_xenobiotics + ".dmnd"
    path_to_ises_db = path_to_ises + ".dmnd"

    # DB yoksa oluştur:
    make_diamond_db(path_to_antibiotics, path_to_antibiotics_db)
    make_diamond_db(path_to_xenobiotics, path_to_xenobiotics_db)
    make_diamond_db(path_to_ises, path_to_ises_db)

    summary = []

    for fasta_file in os.listdir(cycle_folder):
        if not fasta_file.endswith(".fasta") and not fasta_file.endswith(".fa"):
            continue

        fasta_path = os.path.join(cycle_folder, fasta_file)
        base_name = os.path.splitext(fasta_file)[0]
        out_prefix = os.path.join(temp_folder, base_name)

        proteins_faa = out_prefix + ".proteins.faa"
        genes_fna = out_prefix + ".genes.fna"
        gff_file = out_prefix + ".gff"
        filtered_proteins_faa = out_prefix + ".proteins.filtered.faa"

        print(f"Running Prodigal on {fasta_file} with GFF output...")
        run_prodigal_with_gff(fasta_path, proteins_faa, genes_fna, gff_file)


        cds_list = parse_gff_cds(gff_file)
        filtered_cds = filter_overlapping_cds(cds_list, overlap_threshold=0.5)

        print(f"Filtering proteins by overlap threshold, kept {len(filtered_cds)} CDS from {len(cds_list)}")
        filter_protein_fasta(proteins_faa, filtered_cds, filtered_proteins_faa)

        diamond_out_antibiotics = out_prefix + ".antibiotics.tsv"
        diamond_out_xenobiotics = out_prefix + ".xenobiotics.tsv"
        diamond_out_ises = out_prefix + ".ises.tsv"

        print(f"Running Diamond on {base_name} proteins against Antibiotics DB...")
        run_diamond(proteins_faa, path_to_antibiotics, diamond_out_antibiotics)
        print(f"Running Diamond on {base_name} proteins against Xenobiotics DB...")
        run_diamond(proteins_faa, path_to_xenobiotics, diamond_out_xenobiotics)
        print(f"Running Diamond on {base_name} proteins against ISes DB...")
        run_diamond(proteins_faa, path_to_ises, diamond_out_ises)

        protein_to_contig = {}
        query_lengths = {}
        for record in SeqIO.parse(proteins_faa, "fasta"):
            header = record.id  # tüm header, örn: Cycle_1-len18725-_1_23
            contig_id = header.rsplit('_', 1)[0]  # sondan ilk '_' e kadar al
            protein_to_contig[header] = contig_id
            query_lengths[header] = len(record.seq)

        def count_hits_and_list(diamond_tsv, query_lengths, coverage_threshold=0.5):
            contig_counts = {}
            contig_hits = {}
            if not os.path.exists(diamond_tsv):
                return contig_counts, contig_hits
            with open(diamond_tsv) as f:
                for line in f:
                    if not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    # Diamond outfmt 6 kolonları sabit sırada
                    qseqid = parts[0]
                    sseqid = parts[1]
                    align_length = int(parts[3])
                    qstart = int(parts[6])
                    qend = int(parts[7])
                    qlen = query_lengths.get(qseqid)
                    if qlen is None:
                        continue
                    coverage = align_length / qlen
                    if coverage < coverage_threshold:
                        continue
                    contig = protein_to_contig.get(qseqid)
                    if contig is None:
                        continue
                    contig_counts[contig] = contig_counts.get(contig, 0) + 1
                    contig_hits.setdefault(contig, set()).add(sseqid)
            return contig_counts, contig_hits


        antibiotic_counts, antibiotic_hits = count_hits_and_list(diamond_out_antibiotics, query_lengths, coverage_threshold=0.5)
        xenobiotic_counts, xenobiotic_hits = count_hits_and_list(diamond_out_xenobiotics, query_lengths, coverage_threshold=0.5)
        ises_counts, ises_hits = count_hits_and_list(diamond_out_ises, query_lengths, coverage_threshold=0.5)


        all_contigs = set(list(antibiotic_counts.keys()) + list(xenobiotic_counts.keys()) + list(ises_counts.keys()))

        for contig in all_contigs:
            summary.append({
                "cycle_file": fasta_file,
                "contig": contig,
                "antibiotics_hits": antibiotic_counts.get(contig, 0),
                "antibiotics_list": ",".join(sorted(antibiotic_hits.get(contig, []))),
                "xenobiotics_hits": xenobiotic_counts.get(contig, 0),
                "xenobiotics_list": ",".join(sorted(xenobiotic_hits.get(contig, []))),
                "ises_hits": ises_counts.get(contig, 0),
                "ises_list": ",".join(sorted(ises_hits.get(contig, [])))
            })

    df_summary = pd.DataFrame(summary)
    summary_csv = os.path.join(picota_out_folder, "cycles_contig_hits_summary.csv")
    df_summary.to_csv(summary_csv, index=False)
    print(f"Summary saved to {summary_csv}")
