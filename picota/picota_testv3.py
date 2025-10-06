import os
import glob
import logging
import argparse
import subprocess
import csv

from src.cycle_finderv2 import cycle_analysis
from src.sra_download import run_sra_down
from src.assembly import assembly_main
from src.scoringv4ProtBlast import scoring_main
from src.split_cycle_coords_for_is import split_cycles_from_picota
from src.bam_analyse import bam_file_analyze
from src.analyze_blocksv3 import analyze_blocks
from src.logger_setup import setup_logger_from_config
from src.config_loader import load_config, Config, LoggingConfig, ToleranceConfig




logger: logging.Logger = None

'''
# --- Global ---
cfg: Config = load_config("picota/config.yaml")
logger: logging.Logger = setup_logger_from_config(cfg)
'''

# --- SRA pairs can be short, long sra ids ---
def load_sra_pairs(sra_id_file: str):
    pairs = []
    with open(sra_id_file) as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            short_id = row["sra_short_id"].strip()
            val = row.get("sra_long_id")  # Can be none
            if val is not None:
                val = val.strip()
            long_id = val if val not in (None, "", "-", "null") else None
            pairs.append((short_id, long_id))
    return pairs


# --- Minimap2 mapping ---
def run_minimap2(ref_fasta, fastq_file, bam_out="mapping.bam", threads=4, run_dir=None):
    if run_dir is None:
        run_dir = "tmp_mapping"
    os.makedirs(run_dir, exist_ok=True)

    sam_out = os.path.join(run_dir, os.path.basename(bam_out).replace(".bam", ".sam"))
    bam_out = os.path.join(run_dir, os.path.basename(bam_out))
    sorted_bam = bam_out.replace(".bam", "_sorted.bam")
    bai_out = sorted_bam + ".bai"

    if os.path.exists(sorted_bam) and os.path.exists(bai_out):
        logger.info(f"[OK] {sorted_bam} ve {bai_out} bulundu, mapping atlandı.")
        return sorted_bam

    if not os.path.exists(sam_out):
        cmd = ["minimap2", "-ax", "map-ont", "-t", str(threads), ref_fasta, fastq_file]
        logger.info(f"[RUN] Minimap2: {' '.join(cmd)}")
        with open(sam_out, "w") as f:
            subprocess.run(cmd, stdout=f, check=True)

    if not os.path.exists(bam_out):
        subprocess.run(["samtools", "view", "-bS", sam_out, "-o", bam_out], check=True)

    if not os.path.exists(sorted_bam):
        subprocess.run(["samtools", "sort", bam_out, "-o", sorted_bam], check=True)

    if not os.path.exists(bai_out):
        subprocess.run(["samtools", "index", sorted_bam], check=True)

    if os.path.exists(sam_out):
        os.remove(sam_out)

    return sorted_bam


# --- Helper Functions ---
def run_sra_download(acc, out_dir, sra_folder, fastq_dump_path, logger_name):
    os.makedirs(sra_folder, exist_ok=True)
    expected_files = [os.path.join(out_dir, f"{acc}_{i}.fastq") for i in (1, 2)]
    missing = [f for f in expected_files if not os.path.exists(f)]
    if missing:
        logger.info(f"[{acc}] FASTQ eksik, indiriliyor...")
        run_sra_down(acc, out_dir, sra_folder, fastq_dump_path, keep_sra_file=True, the_force=False, logger_name=logger_name)
    else:
        logger.info(f"[{acc}] FASTQ mevcut, atlandı.")
    return [f for f in expected_files if os.path.exists(f)]


def run_longread_download(acc, out_dir, sra_folder, fastq_dump_path, logger_name):
    os.makedirs(sra_folder, exist_ok=True)
    fastq_file = os.path.join(out_dir, f"{acc}_1.fastq")
    if not os.path.exists(fastq_file):
        logger.info(f"[{acc}] Long-read FASTQ indiriliyor...")
        run_sra_down(acc, out_dir, sra_folder, fastq_dump_path, keep_sra_file=True, the_force=False, logger_name=logger_name)
    else:
        logger.info(f"[{acc}] Long-read FASTQ mevcut, atlandı.")
    return fastq_file


# --- Assembly ---
def run_assembly(acc, raw_files, out_folder, cfg: Config):
    gfa_files = glob.glob(os.path.join(out_folder, "*.gfa"))
    if gfa_files:
        logger.info(f"[{acc}] Assembly zaten var, atlandı.")
        return gfa_files

    logger.info(f"[{acc}] Assembly başlatılıyor...")
    assembly_main(
        acc, raw_files, out_folder,
        cfg.paths.assembly_threads, cfg.paths.assembly_k_mer_list, cfg.paths.quiet,
        cfg.paths.keep_temp_files, cfg.paths.path_of_spades,
        cfg.paths.path_of_fastp, cfg.paths.skip_filtering,
        cfg.paths.assembler_type, cfg.paths.path_of_megahit,
        cfg.paths.gfa_tools_path, cfg.paths.path_of_bandage
    )
    return glob.glob(os.path.join(out_folder, "*.gfa"))


# --- Cycle Analysis ---
def run_cycle_analysis(acc, gfa_file, out_file, cfg: Config):
    if os.path.exists(out_file):
        logger.info(f"[{acc}] Cycle analysis zaten yapılmış.")
        return
    logger.info(f"[{acc}] Cycle analysis başlatılıyor...")
    cycle_analysis(
        gfa_file, out_file, cfg.paths.find_all_path, cfg.paths.path_limit,
        cfg.options.min_size_of_cycle, cfg.options.max_size_of_cycle, cfg.options.name_prefix_cycle,
        cfg.options.min_component_number, cfg.options.max_component_number,
        cfg.options.k_mer_sim, cfg.options.threshold_sim
    )


# --- Scoring ---
def run_scoring(acc, cycle_file, out_folder, cfg: Config):
    picota_final_tab = os.path.join(out_folder, 'picota_final_tab')
    if os.path.exists(picota_final_tab):
        logger.info(f"[{acc}] Scoring zaten yapılmış.")
        return picota_final_tab

    logger.info(f"[{acc}] Scoring başlatılıyor...")
    scoring_main(
        cycle_file, out_folder,
        cfg.paths.path_to_antibiotics, cfg.paths.path_to_xenobiotics, cfg.paths.path_to_ises, 
        cfg.options.mean_of_CompTns,
        cfg.options.std_of_CompTns,
        cfg.options.total_score_type,
        cfg.options.threshold_final_score,
        cfg.options.max_z,
        cfg.options.dist_type,
        cfg.options.path_of_prodigal,
        cfg.options.path_of_blastn,
        cfg.options.path_of_makeblastdb,
        cfg.options.path_of_blastx,
        cfg.options.path_of_blastp,
        cfg.logging.logger_name
    )
    return picota_final_tab


# --- Pipeline per accession ---
def process_accession(short_acc, long_acc, cfg: Config):
    logger.info(f"=== Başlatıldı: {short_acc} ===")
    project_root = cfg.paths.outdir

    sra_folder = os.path.join(project_root, "raw", short_acc)
    asm_folder = os.path.join(project_root, "assembly", short_acc)
    cyc_folder = os.path.join(project_root, "cycles")
    scr_folder = os.path.join(project_root, "scoring", short_acc)
    annot_folder = os.path.join(project_root, "annot", short_acc)
    os.makedirs(sra_folder, exist_ok=True)
    os.makedirs(asm_folder, exist_ok=True)
    os.makedirs(cyc_folder, exist_ok=True)
    os.makedirs(scr_folder, exist_ok=True)
    os.makedirs(annot_folder, exist_ok=True)

    # 1) SRA download
    
    raw_files = run_sra_download(short_acc, asm_folder, sra_folder, cfg.paths.fastq_dump, cfg.logging.logger_name)

    # 2) Assembly
    gfa_files = run_assembly(short_acc, raw_files, asm_folder, cfg)
    if not gfa_files:
        logger.warning(f"[{short_acc}] Assembly başarısız, GFA bulunamadı.")
        return

    # 3) Cycle analysis
    gfa_file = gfa_files[0]
    out_cycle_file = os.path.join(cyc_folder, f"{short_acc}_{os.path.basename(gfa_file)}.fasta")
    run_cycle_analysis(short_acc, gfa_file, out_cycle_file, cfg)

    # 4) Scoring
    picota_final_tab = run_scoring(short_acc, out_cycle_file, scr_folder, cfg)
    annotated_fastas = split_cycles_from_picota(picota_final_tab, out_cycle_file, annot_folder)

    # 5) Long-read download + mapping
    long_fastq = None
    if long_acc:
        map_folder = os.path.join(project_root, "mapping", long_acc)
        os.makedirs(map_folder, exist_ok=True)
        long_sra_folder = os.path.join(project_root, "raw_long", long_acc)
        os.makedirs(long_sra_folder, exist_ok=True)

        try:
            long_fastq = run_longread_download(long_acc, map_folder, long_sra_folder, cfg.paths.fastq_dump, cfg.logging.logger_name)
        except Exception as e:
            logger.error(f"[{long_acc}] Long-read download sırasında hata: {e}")

    for fasta_record in annotated_fastas:
        if long_fastq:
            try:
                sorted_bam = run_minimap2(
                    fasta_record, long_fastq,
                    bam_out=f"{long_acc}_{os.path.basename(fasta_record)}_mapping.bam",
                    threads=cfg.options.mapping_threads,
                    run_dir=os.path.join(project_root, "mapping", long_acc)
                )
                out_bam_analysis = os.path.join(map_folder, f"{long_acc}_{os.path.basename(fasta_record)}_full")
                out_bam_analysis_partial = os.path.join(map_folder, f"{long_acc}_{os.path.basename(fasta_record)}_partial")
                bam_file_analyze(sorted_bam, out_bam_analysis, out_bam_analysis_partial)
                analyze_blocks(out_bam_analysis_partial, out_bam_analysis, map_folder, os.path.basename(fasta_record), cfg)
            except Exception as e:
                logger.error(f"[{long_acc}] Mapping veya BAM analiz sırasında hata: {e}")

    logger.info(f"✅ Tamamlandı: {short_acc}")


# --- Main ---
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="YAML config dosyası", default="picota/config.yaml")
    args = parser.parse_args()


    print(args.config)
    cfg: Config = load_config(args.config)

    os.makedirs(cfg.paths.outdir, exist_ok=True)
    sra_pairs = load_sra_pairs(cfg.paths.sra_id_file)
    global logger
    logger = setup_logger_from_config(cfg)
    logger.debug("DEBUG control")

    for short_id, long_id in sra_pairs:
        process_accession(short_id, long_id, cfg)

    logger.info("🎉 Tüm işlemler tamamlandı.")


if __name__ == "__main__":
    main()
