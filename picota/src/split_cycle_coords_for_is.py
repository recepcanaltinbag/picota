from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pandas as pd

def split_cycles_from_picota(picota_tab, cycle_fasta_file, output_dir, min_score=80, min_cargo_gap=5):
    """
    PICOTA final tab ve cycle FASTA dosyasını alır,
    her cycle için transposon ve cargo bölgelerini ayırır (tek cargo segmenti: en büyük olan),
    output_dir altında FASTA dosyaları oluşturur ve liste olarak döner.

    Args:
        picota_tab (str): PICOTA final tab dosya yolu
        cycle_fasta_file (str): Cycle FASTA dosya yolu
        output_dir (str): Split FASTA dosyalarının kaydedileceği klasör
        min_score (float): Üçüncü score için minimum değer; altındakiler atlanır
        min_cargo_gap (int): Çok küçük cargo segmentleri göz ardı etmek için threshold

    Returns:
        list: Oluşturulan FASTA dosyalarının tam yolları
    """
    os.makedirs(output_dir, exist_ok=True)
    fasta_records = {rec.id: rec for rec in SeqIO.parse(cycle_fasta_file, "fasta")}

    df = pd.read_csv(picota_tab, sep="\t", header=0)
    created_fasta_files = []

    for idx, row in df.iterrows():
        cycle_id = row['CycleID']

        # Üçüncü score kontrolü
        if pd.isna(row['score2']) or float(row['score2']) < min_score:
            continue

        # IS koordinatları yoksa atla
        if pd.isna(row['IScoords']) or row['IScoords'] == "":
            continue

        # Cycle uzunluğunu çek
        try:
            seq_len = int(cycle_id.split('-len')[-1].replace('-', ''))
        except:
            continue

        if cycle_id not in fasta_records:
            continue

        seq_record = fasta_records[cycle_id]
        seq_str = str(seq_record.seq)

        # IS koordinatlarını tuple listesine çevir ve sırala
        is_coords_list = row['IScoords'].split(';')
        is_ranges = [tuple(map(int, c.split('-'))) for c in is_coords_list]
        is_ranges_sorted = sorted(is_ranges)

        # Transposon segmentleri
        transposon_segments = [seq_str[start-1:end] for start, end in is_ranges_sorted]

        # Cargo segmenti: transposon dışındaki tüm segmentler
        cargo_candidates = []
        prev_end = 0
        for start, end in is_ranges_sorted:
            if start-1 - prev_end > min_cargo_gap:
                cargo_candidates.append(seq_str[prev_end:start-1])
            prev_end = end
        # son transposon sonrası
        if seq_len - prev_end > min_cargo_gap:
            cargo_candidates.append(seq_str[prev_end:seq_len])

        # Tek cargo segmenti: en büyük olan
        cargo_seq = max(cargo_candidates, key=len) if cargo_candidates else ""

        # FASTA dosyası oluştur
        out_file = os.path.join(output_dir, f"{cycle_id}_split.fasta")
        with open(out_file, "w") as f:
            # Transposon yaz
            for i, seg in enumerate(transposon_segments):
                if len(seg) == 0:
                    continue
                rec = SeqRecord(Seq(seg), id=f"{cycle_id}_transposon_{i+1}", description="")
                SeqIO.write(rec, f, "fasta")
            # Cargo yaz (tek segment: en büyük)
            if len(cargo_seq) > 0:
                rec = SeqRecord(Seq(cargo_seq), id=f"{cycle_id}_cargo", description="")
                SeqIO.write(rec, f, "fasta")

        created_fasta_files.append(out_file)

    return created_fasta_files
