

import os
import glob
import logging
import argparse
import yaml
import subprocess
import csv

from src.cycle_finderv2 import cycle_analysis
from src.sra_download import run_sra_down
from src.assembly import assembly_main
from src.scoringv3_blast import scoring_main




def load_sra_pairs(sra_id_file):
    pairs = []
    with open(sra_id_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            short_id = row["sra_short_id"].strip()
            long_id = row["sra_long_id"].strip() if row["sra_long_id"].strip() not in ("", "-", "null") else None
            pairs.append((short_id, long_id))
    return pairs



def run_minimap2(ref_fasta, fastq_file, bam_out="mapping.bam", threads=4, run_dir=None):
    """Minimap2 ile mapping yapar ve sorted BAM döner.
    Temp dosyalar ayrı klasörde tutulur. Var olanlar tekrar kullanılmaz."""
    
    if run_dir is None:
        run_dir = "tmp_mapping"
    os.makedirs(run_dir, exist_ok=True)

    sam_out = os.path.join(run_dir, os.path.basename(bam_out).replace(".bam", ".sam"))
    bam_out = os.path.join(run_dir, os.path.basename(bam_out))
    sorted_bam = bam_out.replace(".bam", "_sorted.bam")
    bai_out = sorted_bam + ".bai"

    # Eğer sorted bam + index varsa hiç çalıştırma
    if os.path.exists(sorted_bam) and os.path.exists(bai_out):
        print(f"[OK] {sorted_bam} ve {bai_out} bulundu, mapping atlandı.")
        return sorted_bam

    # Minimap2
    if not os.path.exists(sam_out):
        cmd = [
            "minimap2", "-ax", "map-ont", "-t", str(threads),
            ref_fasta, fastq_file
        ]
        print(f"[RUN] Minimap2: {' '.join(cmd)}")
        with open(sam_out, "w") as f:
            subprocess.run(cmd, stdout=f, check=True)

    # SAM -> BAM
    if not os.path.exists(bam_out):
        subprocess.run(["samtools", "view", "-bS", sam_out, "-o", bam_out], check=True)

    # BAM -> sorted BAM
    if not os.path.exists(sorted_bam):
        subprocess.run(["samtools", "sort", bam_out, "-o", sorted_bam], check=True)

    # BAM index
    if not os.path.exists(bai_out):
        subprocess.run(["samtools", "index", sorted_bam], check=True)

    if os.path.exists(sam_out):
        os.remove(sam_out)

    return sorted_bam



# === Logger ===
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.StreamHandler()]
)

# === Helper Functions ===
def run_sra_download(acc, out_dir, sra_folder, fastq_dump_path):
    raw_files = [os.path.join(sra_folder, f"{acc}_{i}.fastq") for i in (1, 2)]
    if all(os.path.exists(f) for f in raw_files):
        logging.info(f"[{acc}] FASTQ zaten mevcut, indirilmiyor.")
    else:
        logging.info(f"[{acc}] FASTQ indiriliyor...")
        run_sra_down(acc, out_dir, sra_folder, fastq_dump_path, keep_sra_file=False, the_force=True)
    return raw_files

def run_longread_download(acc, out_dir, sra_folder, fastq_dump_path):
    fastq_file = os.path.join(sra_folder, f"{acc}.fastq")
    if os.path.exists(fastq_file):
        logging.info(f"[{acc}] Long-read FASTQ zaten mevcut, indirilmiyor.")
    else:
        logging.info(f"[{acc}] Long-read indiriliyor...")
        run_sra_down(acc, out_dir, sra_folder, fastq_dump_path, keep_sra_file=False, the_force=True)
    return fastq_file


def run_assembly(acc, raw_files, out_folder, cfg):
    gfa_files = glob.glob(os.path.join(out_folder, "*.gfa"))
    if gfa_files:
        logging.info(f"[{acc}] Assembly zaten mevcut, atlanıyor.")
        return gfa_files

    logging.info(f"[{acc}] Assembly başlatılıyor...")
    assembly_main(
        acc, raw_files, out_folder,
        cfg["assembly_threads"], cfg["assembly_k_mer_list"], cfg["assembly_quiet"],
        cfg["assembly_keep_temp_files"], cfg["assembly_path_of_spades"],
        cfg["assembly_path_of_fastp"], cfg["assembly_skip_filtering"],
        cfg["assembler_type"], cfg["assembly_path_of_megahit"],
        cfg["gfa_tools_path"], cfg["path_of_bandage"]
    )
    gfa_files = glob.glob(os.path.join(out_folder, "*.gfa"))
    return gfa_files


def run_cycle_analysis(acc, gfa_file, out_file, cfg):
    if os.path.exists(out_file):
        logging.info(f"[{acc}] Cycle analizi zaten yapılmış.")
        return

    logging.info(f"[{acc}] Cycle analizi başlatılıyor...")
    cycle_analysis(
        gfa_file, out_file, cfg["find_all_path"], cfg["path_limit"],
        cfg["min_size_of_cycle"], cfg["max_size_of_cycle"], cfg["name_prefix_cycle"],
        cfg["min_component_number"], cfg["max_component_number"],
        cfg["k_mer_sim"], cfg["threshold_sim"]
    )


def run_scoring(acc, cycle_file, out_folder, cfg):
    logging.info(f"[{acc}] Scoring başlatılıyor...")
    scoring_main(
        cycle_file, out_folder,
        cfg["path_to_antibiotics"], cfg["path_to_xenobiotics"], cfg["path_to_ises"]
    )


# === Pipeline ===
def process_accession(short_acc, long_acc, cfg, project_root):
    logging.info(f"=== Başlatıldı: {short_acc} ===")

    # klasör yapısı
    sra_folder = os.path.join(project_root, "raw", short_acc)
    asm_folder = os.path.join(project_root, "assembly", short_acc)
    cyc_folder = os.path.join(project_root, "cycles")
    scr_folder = os.path.join(project_root, "scoring", short_acc)
    os.makedirs(sra_folder, exist_ok=True)
    os.makedirs(asm_folder, exist_ok=True)
    os.makedirs(cyc_folder, exist_ok=True)
    os.makedirs(scr_folder, exist_ok=True)

    # 1) SRA Download
    raw_files = run_sra_download(short_acc, asm_folder, sra_folder, cfg["path_of_fastq_dump"])

    # 2) Assembly
    gfa_files = run_assembly(short_acc, raw_files, asm_folder, cfg)
    if not gfa_files:
        logging.warning(f"[{short_acc}] Assembly başarısız, GFA bulunamadı.")
        return

    # 3) Cycle Analysis
    gfa_file = gfa_files[0]  # ilk gfa dosyası
    out_cycle_file = os.path.join(cyc_folder, f"{short_acc}_{os.path.basename(gfa_file)}.fasta")
    run_cycle_analysis(short_acc, gfa_file, out_cycle_file, cfg)

    # 4) Scoring
    run_scoring(short_acc, out_cycle_file, scr_folder, cfg)

    # 5) Temizlik
    if cfg.get("delete_fastq_files", False):
        for f in raw_files:
            if os.path.exists(f):
                os.remove(f)
                logging.info(f"[{acc}] Silindi: {f}")

    # long-read varsa mapping
    if long_acc:
        logging.info(f"[{short_acc}] için long-read bulundu: {long_acc}")
        map_folder = os.path.join(project_root, "mapping", short_acc)
        os.makedirs(map_folder, exist_ok=True)
        long_sra_folder = os.path.join(project_root, "raw_long", long_acc)
        os.makedirs(long_sra_folder, exist_ok=True)

        try:
            long_fastq = run_longread_download(
                long_acc, map_folder, long_sra_folder, cfg["path_of_fastq_dump"]
            )
            sorted_bam = run_minimap2(
                out_cycle_file, long_fastq,
                bam_out=f"{long_acc}_mapping.bam",
                threads=cfg.get("mapping_threads", 4),
                run_dir=map_folder
            )
            logging.info(f"[{short_acc}] Long-read mapping tamamlandı: {sorted_bam}")
        except Exception as e:
            logging.error(f"[{short_acc}] Mapping sırasında hata oluştu: {e}")
    else:
        logging.info(f"[{short_acc}] için long-read bulunamadı, mapping atlandı.")


    logging.info(f"✅ Tamamlandı: {short_acc}")






    




logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.StreamHandler()]
)

def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)


def main():
    # --- MAİN İÇİ DEFAULTLAR ---
    default_config = "picota/config.yaml"   # buraya default yaml dosya ismi

    # --- ARGPARSE ---
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="YAML config dosyası", default=default_config)
    args = parser.parse_args()

    cfg = load_config(args.config)
    project_root = cfg.get("outdir")
    os.makedirs(project_root, exist_ok=True)

    # short–long eşleşmeleri oku
    sra_pairs = load_sra_pairs(cfg["sra_id_file"])

    for short_id, long_id in sra_pairs:
        process_accession(short_id, long_id, cfg, project_root)

    logging.info("🎉 Tüm işlemler tamamlandı.")


if __name__ == "__main__":
    main()
