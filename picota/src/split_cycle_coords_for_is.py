from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pandas as pd
import re

def split_cycles_from_picota(picota_tab, cycle_fasta_file, output_dir, split_min_score=50):
    """
    PICOTA final tab ve cycle FASTA dosyasını alır,
    her cycle için sadece 1 transposon ve 1 cargo bölgesi ayırır.
    Cargo sınırları Antibiotic/Xenobiotic genlere göre belirlenir.
    Transposon başta veya sonda olabilir.
    Eğer cargo baştaysa bitiş transposonun başlangıcı olur.
    Eğer cargo sondaysa başlangıç transposonun bitişinden hemen sonra olur.
    Output FASTA ID formatı: CycleID_start_end_transposon / cargo
    """
    os.makedirs(output_dir, exist_ok=True)
    fasta_records = {rec.id: rec for rec in SeqIO.parse(cycle_fasta_file, "fasta")}
    df = pd.read_csv(picota_tab, sep="\t", header=0)
    created_fasta_files = []
    for idx, row in df.iterrows():
        cycle_id = row['CycleID']
        # Skor kontrolü
        if pd.isna(row['score2']) or float(row['score2']) < split_min_score:
            continue
        if pd.isna(row['IScoords']) or row['IScoords'] == "":
            continue

        match = re.search(r'-len(\d+)-comp(\d+)-', cycle_id)
        if match:
            seq_len = int(match.group(1))
            comp_number = int(match.group(2))
        else:
            continue  # format uymuyorsa atla
        if cycle_id not in fasta_records:
            continue

        seq_record = fasta_records[cycle_id]
        seq_str = str(seq_record.seq)

        # === Transposon koordinatları ===
        is_coords_list = row['IScoords'].split(';')
        is_ranges = [tuple(map(int, c.split('-'))) for c in is_coords_list if "-" in c]
        if not is_ranges:
            continue
        is_ranges_sorted = sorted(is_ranges)
        trans_start, trans_end = is_ranges_sorted[0][0], is_ranges_sorted[-1][1]

        # === Antibiotic + Xenobiotic koordinatları birleşimi ===
        ant_xeno_coords = []
        if not pd.isna(row.get("Antcoords")) and row["Antcoords"] != "":
            for c in str(row["Antcoords"]).split(";"):
                if "-" in c:
                    s, e = map(int, c.split("-"))
                    ant_xeno_coords.append((s, e))
        if not pd.isna(row.get("Xenocoords")) and row["Xenocoords"] != "":
            for c in str(row["Xenocoords"]).split(";"):
                if "-" in c:
                    s, e = map(int, c.split("-"))
                    ant_xeno_coords.append((s, e))

        cargo_start, cargo_end = None, None

        # === Cargo sınırları ===
        if ant_xeno_coords:
            min_gene = min(s for s, e in ant_xeno_coords)
            max_gene = max(e for s, e in ant_xeno_coords)

            if min_gene > trans_end:
                # Cargo sondadır → başı transposon bitişi
                cargo_start = trans_end + 1
                cargo_end = max_gene
            elif max_gene < trans_start:
                # Cargo baştadır → sonu transposon başlangıcı
                cargo_start = min_gene
                cargo_end = trans_start - 1
        else:
            # gen yoksa → transposondan önce veya sonra en büyük parçayı seç
            left_len = trans_start - 1
            right_len = seq_len - trans_end
            if right_len >= left_len and right_len > 0:
                cargo_start, cargo_end = trans_end + 1, seq_len
            elif left_len > right_len:
                cargo_start, cargo_end = 1, trans_start - 1

        # === FASTA dosyası yaz ===
        out_file = os.path.join(output_dir, f"{cycle_id}_split.fasta")
        with open(out_file, "w") as f:
            # Transposon
            trans_seq = seq_str[trans_start-1:trans_end]
            rec_trans = SeqRecord(
                Seq(trans_seq),
                id=f"{cycle_id}_{trans_start}_{trans_end}_transposon",
                description=""
            )
            SeqIO.write(rec_trans, f, "fasta")

            # Cargo
            if cargo_start and cargo_end and cargo_end > cargo_start:
                cargo_seq = seq_str[cargo_start-1:cargo_end]
                rec_cargo = SeqRecord(
                    Seq(cargo_seq),
                    id=f"{cycle_id}_{cargo_start}_{cargo_end}_cargo",
                    description=""
                )
                SeqIO.write(rec_cargo, f, "fasta")

        created_fasta_files.append(out_file)

    return created_fasta_files
