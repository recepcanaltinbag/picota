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

def run_prodigal(input_fasta, output_proteins, output_genes):
    cmd = [
        "prodigal",
        "-i", input_fasta,
        "-a", output_proteins,
        "-d", output_genes,
        "-p", "meta"
    ]
    subprocess.run(cmd, check=True)

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

        print(f"Running Prodigal on {fasta_file}...")
        run_prodigal(fasta_path, proteins_faa, genes_fna)

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
        for record in SeqIO.parse(proteins_faa, "fasta"):
            header = record.id  # tüm header, örn: Cycle_1-len18725-_1_23
            contig_id = header.rsplit('_', 1)[0]  # sondan ilk '_' e kadar al
            protein_to_contig[header] = contig_id

        def count_hits(diamond_tsv):
            contig_counts = {}
            if not os.path.exists(diamond_tsv):
                return contig_counts
            with open(diamond_tsv) as f:
                for line in f:
                    if not line.strip():
                        continue
                    qseqid = line.split('\t')[0]
                    contig = protein_to_contig.get(qseqid)
                    if contig is None:
                        continue
                    contig_counts[contig] = contig_counts.get(contig, 0) + 1
            return contig_counts

        antibiotic_counts = count_hits(diamond_out_antibiotics)
        xenobiotic_counts = count_hits(diamond_out_xenobiotics)
        ises_counts = count_hits(diamond_out_ises)

        all_contigs = set(list(antibiotic_counts.keys()) + list(xenobiotic_counts.keys()) + list(ises_counts.keys()))

        for contig in all_contigs:
            summary.append({
                "cycle_file": fasta_file,
                "contig": contig,
                "antibiotics_hits": antibiotic_counts.get(contig, 0),
                "xenobiotics_hits": xenobiotic_counts.get(contig, 0),
                "ises_hits": ises_counts.get(contig, 0)
            })

    df_summary = pd.DataFrame(summary)
    summary_csv = os.path.join(picota_out_folder, "cycles_contig_hits_summary.csv")
    df_summary.to_csv(summary_csv, index=False)
    print(f"Summary saved to {summary_csv}")
