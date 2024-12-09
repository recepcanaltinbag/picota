import subprocess
import glob
import os
import shutil


# ----------------------------------------------------------------
# 1- FILTERING THE RAW FILE WITH FASTP
# ----------------------------------------------------------------
'''
fastp -i in.fq -o out.fq
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

'''
# file path is list
def raw_read_filtering(raw_file, out_folder, fastp_path, q, reads_to_process, quiet_mode):

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    if len(raw_file) == 1:
        raw_file_1_path = raw_file[0]
        filtered_raw_file_1_name = "filtered_" + raw_file_1_path.split('/')[-1].split('.')[0] 
        f1_path = out_folder + "/" + filtered_raw_file_1_name +".fastq"
        args = f"{fastp_path} -i {raw_file_1_path} -o {f1_path} -h {f1_path}.html -q {q} --reads_to_process {reads_to_process}"
        print('Command will be run:')
        print(args)
        print('-------')
        my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True, capture_output=quiet_mode)


    elif len(raw_file) == 2:
        raw_file_1_path = raw_file[0]
        raw_file_2_path = raw_file[1]
        filtered_raw_file_1_name = "filtered_" + raw_file_1_path.split('/')[-1].split('.')[0]
        filtered_raw_file_2_name = "filtered_" + raw_file_2_path.split('/')[-1].split('.')[0]
        f1_path = out_folder + "/" + filtered_raw_file_1_name +".fastq"
        f2_path = out_folder + "/" + filtered_raw_file_2_name +".fastq"
        args = f"{fastp_path} -i {raw_file_1_path} -I {raw_file_2_path} \
            -o {f1_path} -O {f2_path} -h {f1_path}.html -q {q} --reads_to_process {reads_to_process}" 
        print('Command will be run:')
        print(args)
        print('-------')
        my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True, capture_output=quiet_mode)


    else:
        print('More than 2 file or no file exist, error!')



#k_mer_list = '39,59,79,99'
def assembly_driver_spades(spades_path, file_path, out_folder, gfa_folder, gfa_name, threads, k_mer, quiet_mode, assembly_keep_temp_files):
    if len(file_path) == 1:
        args = f"{spades_path} -1 {file_path[0]} -o {out_folder} -t {str(threads)}  -k {k_mer}" 
        print('Command will be run:')
        print(args)
        print('-------')
        my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True, capture_output=quiet_mode)
        gfa_files = glob.glob(out_folder + '/*.gfa')
        if len(gfa_files) == 0:
            raise Exception('There is no gfa files, there can be error in assembly process')
        else:
            args = f"cp {gfa_files[0]} {os.path.join(gfa_folder, gfa_name)}"
            my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True)

        if assembly_keep_temp_files == False:
            shutil.rmtree(out_folder)
            for file_pt in file_path:
                if os.path.exists(file_pt):
                    os.remove(file_pt)
            print('Temp Files deleted., if you want to keep them use --keep_temp_files')
                

    elif len(file_path) == 2:
        args = f"{spades_path} -1 {file_path[0]} -2 {file_path[1]} -o {out_folder} -t {str(threads)} -k {k_mer}" 
        print('Command will be run:')
        print(args)
        print('-------')
        my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True, capture_output=quiet_mode)
        gfa_files = glob.glob(out_folder + '/*.gfa')

        if len(gfa_files) == 0:
            raise Exception('There is no gfa files, there can be error in assembly process')
        else:
            args = f"cp {gfa_files[0]} {os.path.join(gfa_folder, gfa_name)}"
            my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True)
        
        

    else:
        print('Error, there is no fastq file or more than two fastq file!')




def assembly_main(name_for_assembly, raw_file_list, main_out_folder, assembly_threads, assembly_k_mer_list, \
        assembly_quiet, assembly_reads_to_process, assembly_fastp_q, assembly_keep_temp_files, \
        assembly_path_of_spades, assembly_path_of_fastp, assembly_skip_filtering):
    
    if not os.path.exists(main_out_folder):
        os.makedir(main_out_folder)
    
    the_final_file_list = []
    if assembly_skip_filtering == False:
        out_folder_for_filtering = os.path.join(main_out_folder, 'filtering')
        if not os.path.exists(out_folder_for_filtering):
            os.mkdir(out_folder_for_filtering)
        
        out_filtering = os.path.join(out_folder_for_filtering, name_for_assembly)
        raw_read_filtering(raw_file_list, out_filtering, assembly_path_of_fastp, assembly_fastp_q, assembly_reads_to_process, assembly_quiet)
        filtered_file_list = glob.glob(out_filtering + '/*.fastq')
        the_final_file_list = filtered_file_list
    else:
        the_final_file_list = raw_file_list

    for k_mer_l in assembly_k_mer_list.split(','):
        gfa_folder = os.path.join(main_out_folder, 'gfa_files')
        if not os.path.exists(gfa_folder):
            os.mkdir(gfa_folder)

        gfa_name = name_for_assembly + '_' + k_mer_l + '.gfa'
        out_assembly_main = os.path.join(main_out_folder, 'assembly')
        if not os.path.exists(out_assembly_main):
            os.mkdir(out_assembly_main)
        
        out_assembly = os.path.join(out_assembly_main, name_for_assembly + '_' + k_mer_l)
        if not os.path.exists(out_assembly):
            os.mkdir(out_assembly)

        if os.path.exists(gfa_folder + '/' + gfa_name):
            print('GFA File exist, skipping,', gfa_folder + '/' + gfa_name)
        else:
            assembly_driver_spades(assembly_path_of_spades, filtered_file_list, out_assembly, gfa_folder, \
            gfa_name, assembly_threads, k_mer_l, assembly_quiet, assembly_keep_temp_files)