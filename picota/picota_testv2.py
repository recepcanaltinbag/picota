from src.cycle_finderv2 import cycle_analysis

from src.sra_download import run_sra_down
from src.assembly import assembly_main
from src.scoringv2 import scoring_main


import os
import glob


# === Global Ayarlar ===
base_dir = "/media/lin-bio/back2/picota_test"
picota_out_folder = "/media/lin-bio/back2/picota_out"
os.makedirs(picota_out_folder, exist_ok=True)

delete_fastq_files = True

# Eğer birden fazla accession varsa buraya yaz
sra_acc_files = [
    "SRR32600069",
    "SRR32600072",
    "SRR29656294",
    "SRR29656296",
    "SRR29656298",
    "SRR18925654",
    "SRR12917028",
    "SRR12917029",
    "SRR12917032",
    "SRR12917033",
    "SRR12917034",
    "SRR12917049",
    "SRR12917050",
    "SRR28272156"
]



# Ortak parametreler
path_of_fastq_dump = "parallel-fastq-dump"
assembly_threads = 24
assembly_k_mer_list = "79,99"
assembly_quiet = False
assembly_keep_temp_files = True
assembly_path_of_spades = "spades.py"
assembly_path_of_fastp = "fastp"
assembly_skip_filtering = True
assembler_type = "megahit"
assembly_path_of_megahit = "megahit"
gfa_tools_path = "./picota/tools/gfaview/misc/fastg2gfa"
path_of_bandage = "./picota/tools/Bandage_Ubuntu-x86-64_v0.9.0_AppDir/Bandage_Ubuntu-x86-64_v0.9.0/usr/bin/bandage"

# Scoring veri tabanları
path_to_antibiotics = "picota/DBs/Antibiotics/protein_fasta_protein_homolog_model.fasta"
path_to_xenobiotics = "picota/DBs/Xenobiotics/Xenobiotics.fasta"
path_to_ises = "picota/DBs/ISes/clusters.single.faa"

# === Pipeline Döngüsü ===
for sra_acc_file in sra_acc_files:
    print(f"\n=== İşlem başlatıldı: {sra_acc_file} ===")
    
    # SRA download

    main_out_folder = os.path.join(base_dir, sra_acc_file)

    sra_folder = os.path.join(main_out_folder, "sra_files")
    os.makedirs(sra_folder, exist_ok=True)

    out_dir = os.path.join(sra_folder)
    raw_file_list = [os.path.join(sra_folder, f"{sra_acc_file}_{i}.fastq") for i in (1, 2)]

    for file_path2 in raw_file_list:
        if os.path.exists(file_path2):
            print('fastq exist skip sra')
        else:
            run_sra_down(sra_acc_file, out_dir, sra_folder, path_of_fastq_dump, keep_sra_file=False, the_force=True)
    
    # Assembly
    
    
    assembly_main(
        sra_acc_file, raw_file_list, main_out_folder,
        assembly_threads, assembly_k_mer_list, assembly_quiet,
        assembly_keep_temp_files, assembly_path_of_spades,
        assembly_path_of_fastp, assembly_skip_filtering,
        assembler_type, assembly_path_of_megahit,
        gfa_tools_path, path_of_bandage
    )
    
    # Cycle Analysis (örnek olarak tek bir gfa dosyası kullanılıyor)
    #gfa_file_name = "k119.gfa"  # Not: dinamik hale getirmek istiyorsan os.listdir ile taranabilir
    #path_to_data = os.path.join(base_dir, "gfa_files", "best_gfa", gfa_file_name)
    
    gfa_search_path = os.path.join(base_dir, sra_acc_file, "*.gfa")
    gfa_files = glob.glob(gfa_search_path)
    if not gfa_files:
        print(f"⚠️ GFA dosyası bulunamadı: {gfa_search_path}")
        continue  # Bu SRA için işlemi atla
    else:
        if delete_fastq_files:
            for file_path2 in raw_file_list:
                if os.path.exists(file_path2):
                    os.remove(file_path2)
                    print(f"🗑️ Silindi: {file_path2}")
                else:
                    print(f"⚠️ Dosya bulunamadı (zaten silinmiş olabilir): {file_path2}")


    path_to_data = gfa_files[0]  # İlk bulunan GFA dosyasını kullan

    out_cycle_folder = os.path.join(base_dir, "test_cycle")
    os.makedirs(out_cycle_folder, exist_ok=True)
    out_cycle_file = os.path.join(out_cycle_folder, f"{sra_acc_file}_{gfa_file_name}.fasta")
    
    find_all_path = True
    path_limit = 25
    min_size_of_cycle = 2000
    max_size_of_cycle = 50000
    name_prefix_cycle = 'Cycle'
    min_component_number = 1
    max_component_number = 50
    k_mer_sim = 50
    threshold_sim = 50
    
    cycle_analysis(
        path_to_data, out_cycle_file, find_all_path, path_limit,
        min_size_of_cycle, max_size_of_cycle, name_prefix_cycle,
        min_component_number, max_component_number, k_mer_sim, threshold_sim
    )
    
    # Scoring
    scoring_out_folder = os.path.join(picota_out_folder, sra_acc_file)

    scoring_main(
        cycle_file, scoring_out_folder,
        path_to_antibiotics, path_to_xenobiotics, path_to_ises
    )

    print(f"✅ Tamamlandı: {sra_acc_file}")

print("🎉 Tüm işlemler tamamlandı.")
