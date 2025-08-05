from src.cycle_finderv2 import cycle_analysis

from src.sra_download import run_sra_down
from src.assembly import assembly_main
from src.scoringv2 import scoring_main



path_to_data = "/media/lin-bio/back2/picota_assembly/gfa_files/best_gfa/k119.gfa"
path_to_data = "picota/test_data/testNitro.gfa"
find_all_path = True
path_limit = 25
out_cycle_file = "/media/lin-bio/back2/picota_assembly/test_cycle/testNitro.fasta"

min_size_of_cycle = 2000
max_size_of_cycle = 50000
name_prefix_cycle = 'Cycle'
min_component_number = 1
max_component_number = 50

k_mer_sim = 50
threshold_sim = 50

cycle_analysis(path_to_data, out_cycle_file, find_all_path, path_limit, min_size_of_cycle, max_size_of_cycle,\
                name_prefix_cycle, min_component_number, max_component_number, k_mer_sim, threshold_sim)



cycle_folder = "/media/lin-bio/back2/picota_assembly/test_cycle"
picota_out_folder = "/media/lin-bio/back2/picota_assembly/picota_out"
path_to_antibiotics = "picota/DBs/Antibiotics/protein_fasta_protein_homolog_model.fasta"
path_to_xenobiotics = "picota/DBs/Xenobiotics/Xenobiotics.fasta"
path_to_ises = "picota/DBs/ISes/clusters.single.faa"


scoring_main(cycle_folder, picota_out_folder, path_to_antibiotics, path_to_xenobiotics, path_to_ises)



print('NEd')
input()













sra_acc_file = "SRR12917050"
out_dir = "picota/sra_files" + "/" + sra_acc_file
sra_folder = "picota/sra_files"
path_of_fastq_dump = "parallel-fastq-dump"

'''
run_sra_down(sra_acc_file, out_dir, sra_folder, path_of_fastq_dump, True, True)
'''
print('End')


name_for_assembly = sra_acc_file
raw_file_list = ["picota/sra_files" + "/" + sra_acc_file + "/" + sra_acc_file + "_1.fastq", "picota/sra_files" + "/" + sra_acc_file + "/" + sra_acc_file + "_2.fastq"]
main_out_folder = "/media/lin-bio/back2/picota_assembly"
assembly_threads = 24
assembly_k_mer_list = "79,99"
assembly_quiet = False
assembly_keep_temp_files = True
assembly_path_of_spades = "spades.py"
assembly_path_of_fastp = "fastp"
assembly_skip_filtering = True
assembler_type="megahit"
assembly_path_of_megahit="megahit"
gfa_tools_path = "./picota/tools/gfaview/misc/fastg2gfa"

path_of_bandage = "./picota/tools/Bandage_Ubuntu-x86-64_v0.9.0_AppDir/Bandage_Ubuntu-x86-64_v0.9.0/usr/bin/bandage"


'''
assembly_main(name_for_assembly, raw_file_list, main_out_folder, assembly_threads, assembly_k_mer_list, \
        assembly_quiet, assembly_keep_temp_files, \
        assembly_path_of_spades, assembly_path_of_fastp, assembly_skip_filtering, assembler_type, assembly_path_of_megahit, gfa_tools_path, path_of_bandage)

'''


path_to_data = "/media/lin-bio/back2/picota_assembly/gfa_files/best_gfa/k119.gfa"
path_to_data = "picota/test_data/testNitro.gfa"
find_all_path = True
path_limit = 25
out_cycle_file = "/media/lin-bio/back2/picota_assembly/testNitro.fasta"

min_size_of_cycle = 2000
max_size_of_cycle = 50000
name_prefix_cycle = 'Cycle'
min_component_number = 1
max_component_number = 50

k_mer_sim = 50
threshold_sim = 50

cycle_analysis(path_to_data, out_cycle_file, find_all_path, path_limit, min_size_of_cycle, max_size_of_cycle,\
                name_prefix_cycle, min_component_number, max_component_number, k_mer_sim, threshold_sim)







