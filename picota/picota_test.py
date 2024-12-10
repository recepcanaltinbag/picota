

from src.cycle_finder import cycle_anaylsis


#zst_file()
#input()


path_to_data = "assembly_graph_big.gfa"

find_all_path = True
path_limit = 25

out_cycle_file = "picota/test_data/cyclesOutBig.fasta"

min_size_of_cycle = 3000
max_size_of_cycle = 100000
name_prefix_cycle = 'Cycle'
min_component_number = 1
max_component_number = 25

k_mer_sim = 200
threshold_sim = 99


cycle_anaylsis(path_to_data, out_cycle_file, find_all_path, path_limit, min_size_of_cycle, max_size_of_cycle,\
                name_prefix_cycle, min_component_number, max_component_number, k_mer_sim, threshold_sim)

    
    