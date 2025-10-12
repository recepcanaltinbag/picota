import os
import re
import pandas as pd
from ast import literal_eval

mapping_dir = "mapping"
scoring_dir = "scoring"
output_csv = "merged_summary.csv"


def load_sra_map():
    """selected_sra.csv dosyasını yükler ve long→short dict döndürür"""
    sra_map = {}
    try:
        df = pd.read_csv("selected_sra.csv")
        if "sra_short_id" in df.columns and "sra_long_id" in df.columns:
            # long → short mapping
            sra_map = dict(zip(df["sra_long_id"], df["sra_short_id"]))
            print(f"✅ {len(sra_map)} SRA eşlemesi yüklendi (long→short).")
        else:
            print("⚠️ selected_sra.csv sütunları eksik (sra_short_id, sra_long_id gerekli)")
    except FileNotFoundError:
        print("⚠️ selected_sra.csv bulunamadı, mapping ID'ler doğrudan kullanılacak.")
    return sra_map


def parse_mapping_file(file_path):
    """mapping dosyasından sadece pattern bilgilerini çıkarır"""
    patterns = {}
    with open(file_path, "r") as f:
        lines = f.read().splitlines()

    for line in lines:
        line = line.strip()
        if not line or line.startswith("==="):
            continue

        if line.startswith(("SRR", "DRR", "ERR")):
            if ":" in line and "=" in line:
                right_part = line.split(":", 1)[1]
                pairs = [p.strip() for p in right_part.split(",") if "=" in p]
                for pair in pairs:
                    key, val = pair.split("=")
                    try:
                        patterns[key.strip()] = int(val.strip())
                    except ValueError:
                        pass
        elif line.startswith("{"):
            try:
                dict_data = literal_eval(line)
                if isinstance(dict_data, dict):
                    patterns.update(dict_data)
            except Exception:
                pass
    return patterns


def extract_cycle_id(filename):
    """dosya adından CycleID'yi çıkartır"""
    match = re.search(r"srr_summary_?(Cycle_.+?)_split\.fasta\.txt", filename)
    return match.group(1) if match else None


def collect_mapping_data(sra_map):
    mapping_data = []
    unknown_counter = 1

    for root, dirs, files in os.walk(mapping_dir):
        for file in files:
            if file.endswith(".txt") and "srr_summary" in file:
                path = os.path.join(root, file)
                patterns = parse_mapping_file(path)
                cycle_id = extract_cycle_id(file)

                # mapping klasör adı long id (ör: SRR5884790)
                long_id = os.path.basename(root)

                # eşleşen short id'yi bul
                srr_id = sra_map.get(long_id)
                if not srr_id:
                    srr_id = f"Unknown_SRR_{unknown_counter}"
                    unknown_counter += 1
                    print(f"⚠️ Eşleşme bulunamadı, long id: {long_id} -> {srr_id}")

                if cycle_id:
                    mapping_data.append({
                        "SRR_long": long_id,
                        "SRR": srr_id,  # short ID (scoring klasörüne denk)
                        "CycleID": cycle_id,
                        **patterns
                    })
                else:
                    print(f"⚠️ CycleID çıkarılamadı: {file}")

    df = pd.DataFrame(mapping_data)
    if df.empty:
        print("⚠️ mapping_df boş, SRR veya CycleID eksik olabilir.")
    else:
        print(f"✅ {len(df)} mapping kaydı toplandı.")
    return df


def merge_with_scoring(mapping_df):
    merged = []
    for srr in mapping_df["SRR"].unique():
        mapping_subset = mapping_df[mapping_df["SRR"] == srr]
        scoring_path = os.path.join(scoring_dir, srr, "picota_final_tab")

        if not os.path.exists(scoring_path):
            print(f"⚠️ Skipped (scoring bulunamadı): {srr}")
            continue

        try:
            scoring_df = pd.read_csv(scoring_path, sep="\t")
        except Exception as e:
            print(f"⚠️ Skipped (okuma hatası): {srr} ({e})")
            continue

        merged_df = pd.merge(mapping_subset, scoring_df,
                             on="CycleID", how="left",
                             suffixes=("_pattern", "_score"))
        merged.append(merged_df)

    if merged:
        return pd.concat(merged, ignore_index=True)
    else:
        return pd.DataFrame()


def main():
    print("📥 SRA ID eşlemesi yükleniyor...")
    sra_map = load_sra_map()

    print("🧩 Mapping dosyaları toplanıyor...")
    mapping_df = collect_mapping_data(sra_map)
    print(f"{len(mapping_df)} mapping kaydı bulundu.")
    print(mapping_df.head())

    print("🔗 Scoring dosyalarıyla birleştiriliyor...")
    merged_df = merge_with_scoring(mapping_df)
    print(f"{len(merged_df)} toplam birleşik kayıt elde edildi.")

    merged_df.to_csv(output_csv, index=False)
    print(f"✅ Birleşik tablo '{output_csv}' olarak kaydedildi.")
    return merged_df


if __name__ == "__main__":
    df = main()

