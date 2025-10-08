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
def raw_read_filtering(raw_file, out_folder, fastp_path, quiet_mode):

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    if len(raw_file) == 1:
        raw_file_1_path = raw_file[0]
        filtered_raw_file_1_name = "filtered_" + raw_file_1_path.split('/')[-1].split('.')[0] 
        f1_path = out_folder + "/" + filtered_raw_file_1_name +".fastq"
        args = f"{fastp_path} -i {raw_file_1_path} -o {f1_path} -h {f1_path}.html"
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
            -o {f1_path} -O {f2_path} -h {f1_path}.html" 
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




def assembly_main(name_for_assembly, raw_file_list, main_out_folder, assembly_threads, assembly_k_mer_list,
        assembly_quiet, assembly_keep_temp_files,
        assembly_path_of_spades, assembly_path_of_fastp, assembly_skip_filtering, 
        assembler_type="spades", assembly_path_of_megahit=None, gfa_tools_path="gfatools", path_of_bandage="bandage"):

    
    if not os.path.exists(main_out_folder):
        os.mkdir(main_out_folder)
    
    the_final_file_list = []
    out_folder_for_filtering = os.path.join(main_out_folder, 'filtering')
    if not os.path.exists(out_folder_for_filtering):
        os.mkdir(out_folder_for_filtering)
    
    out_filtering = os.path.join(out_folder_for_filtering, name_for_assembly)
    if assembly_skip_filtering == False:
        raw_read_filtering(raw_file_list, out_filtering, assembly_path_of_fastp, assembly_quiet)
        filtered_file_list = glob.glob(out_filtering + '/*.fastq')
        if len(filtered_file_list) == 0: 
            the_final_file_list = raw_file_list
        else:
            the_final_file_list = filtered_file_list
    else:
        filtered_file_list = glob.glob(out_filtering + '/*.fastq')
        if len(filtered_file_list) == 0: 
            the_final_file_list = raw_file_list
        else:
            the_final_file_list = filtered_file_list
    
    gfa_folder = os.path.join(main_out_folder, 'gfa_files')
    if not os.path.exists(gfa_folder):
        os.mkdir(gfa_folder)
    out_assembly_main = os.path.join(main_out_folder, 'assembly')
    if not os.path.exists(out_assembly_main):
        os.mkdir(out_assembly_main)
    
    
    if assembler_type.lower() == "megahit":

        gfa_name = name_for_assembly + '_' + 'megahit' + '.gfa'
        out_assembly = os.path.join(out_assembly_main, name_for_assembly + '_' + 'megahit')
        if os.path.exists(out_assembly):
            shutil.rmtree(out_assembly)
            print(f"{out_assembly} deleted, overwrite will be truw.")

        if os.path.exists(gfa_folder + '/' + gfa_name):
            print('GFA File exist, skipping,', gfa_folder + '/' + gfa_name)
        else:
            assembly_driver_megahit(assembly_path_of_megahit, the_final_file_list, out_assembly, main_out_folder, gfa_folder, gfa_name, assembly_threads, assembly_quiet, assembly_keep_temp_files, gfa_tools_path, path_of_bandage)
    else:
        gfa_files = []
        for k_mer_l in assembly_k_mer_list.split(','):
            gfa_name = name_for_assembly + '_' + k_mer_l + '.gfa'
            out_assembly = os.path.join(out_assembly_main, name_for_assembly + '_' + k_mer_l)
            if not os.path.exists(out_assembly):
                os.mkdir(out_assembly)
            gfa_files.append(gfa_folder + '/' + gfa_name)
            if os.path.exists(gfa_folder + '/' + gfa_name):
                print('GFA File exist, skipping,', gfa_folder + '/' + gfa_name)
            else:
                assembly_driver_spades(assembly_path_of_spades, the_final_file_list, out_assembly, gfa_folder, \
                gfa_name, assembly_threads, k_mer_l, assembly_quiet, assembly_keep_temp_files)

        best_gfa = process_gfa_files(gfa_files, path_of_bandage)
        if best_gfa:
            # Dosyanın bulunduğu dizini al
            destination_path = os.path.join(main_out_folder, os.path.basename(best_gfa))
            shutil.copy(best_gfa, destination_path)

            print(f"Best GFA copied to: {destination_path}")
        else:
            print("No best GFA file found to copy.")





#MEGAHIT


def assembly_driver_megahit(megahit_path, file_path, out_folder, main_out_folder, gfa_folder, gfa_name, threads, quiet_mode, assembly_keep_temp_files, gfa_tools_path, path_of_bandage):
    
    
    
    
    if len(file_path) == 1:
        args = f"{megahit_path} -r {file_path[0]} -o {out_folder} -t {str(threads)} --k-min 55"
    elif len(file_path) == 2:
        args = f"{megahit_path} -1 {file_path[0]} -2 {file_path[1]} -o {out_folder} -t {str(threads)} --k-min 55"
    else:
        print('Error: there is no fastq file or more than two fastq file!')
        return

    print('Command will be run:')
    print(args)
    print('-------')
    my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True, capture_output=quiet_mode)

    intermediate_dir = os.path.join(out_folder, "intermediate_contigs")
    contig_files = glob.glob(os.path.join(intermediate_dir, "k*.contigs.fa"))

    gfa_files = []
    for contig_file in contig_files:
        if os.path.getsize(contig_file) == 0:
            print(f"Skipping {contig_file}: file is empty.")
            continue
        
        kmer = os.path.basename(contig_file).split('.')[0]
        fastg_path = os.path.join(gfa_folder, f"{kmer}.fastg")
        gfa_path = os.path.join(gfa_folder, f"{kmer}.gfa")

        if not os.path.isfile(gfa_path):
            cmd = f"megahit_core contig2fastg {kmer[1:]} {contig_file} > {fastg_path}"
            print(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)

            cmd = f"{gfa_tools_path} {fastg_path} > {gfa_path}"
            print(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)

        # Skoru GFA'dan hesapla
        gfa_files.append(gfa_path)


    best_gfa = process_gfa_files(gfa_files, path_of_bandage)

    if best_gfa:
        # Dosyanın bulunduğu dizini al
        destination_path = os.path.join(main_out_folder, os.path.basename(best_gfa))
        shutil.copy(best_gfa, destination_path)

        print(f"Best GFA copied to: {destination_path}")
    else:
        print("No best GFA file found to copy.")


    if assembly_keep_temp_files == False:
        shutil.rmtree(out_folder)
        for file_pt in file_path:
            if os.path.exists(file_pt):
                os.remove(file_pt)
        print('Temp Files deleted., if you want to keep them use --keep_temp_files')





def process_gfa_files(gfa_files, path_of_bandage):
    scores = []

    for gfa_file in gfa_files:
        if os.path.getsize(gfa_file) == 0:
            print(f"Skipping empty GFA: {gfa_file}")
            continue

        gfa_path, contigs, dead_ends, score = compute_score_from_gfa(gfa_file, path_of_bandage)
        scores.append((gfa_path, contigs, dead_ends, score))

    if not scores:
        print("No valid GFA files found.")
        return None

    # Normalize skorlar
    score_values = [s[3] for s in scores]
    min_score = min(score_values)
    max_score = max(score_values)
    range_score = max_score - min_score if max_score != min_score else 1

    normalized_scores = [
        (path, contigs, dead, (score - min_score) / range_score)
        for path, contigs, dead, score in scores
    ]

    # Tüm normalize skorları yazdır
    print("\nGFA Scores (Normalized):")
    for path, contigs, dead_ends, norm_score in normalized_scores:
        print(f"{os.path.basename(path)} | contigs: {contigs}, dead ends: {dead_ends}, normalized score: {norm_score:.4f}")

    # En yüksek normalize skora sahip GFA dosyasını bul
    best_gfa = max(normalized_scores, key=lambda x: x[3])
    print(f"\nBest GFA: {os.path.basename(best_gfa[0])} with score {best_gfa[3]:.4f}")
    return best_gfa[0]





def compute_score_from_gfa(gfa_file, path_of_bandage):
    # Bandage ile dead ends sayısını al
    print('gfa_process:', gfa_file)
    cmd = f'{path_of_bandage} info {gfa_file} | grep "Dead ends" | grep -oP "\\d+"'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable='/bin/bash')
    try:
        dead_ends = int(result.stdout.strip())
    except ValueError:
        dead_ends = 0

    # Contig ve edge (bağlantı) sayılarını hesapla
    contigs = 0
    edges = 0
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                contigs += 1
            elif line.startswith('L'):
                edges += 1

    # Skoru hesapla: bağlantı sayısı / (contig sayısı * (dead_ends + 1)^2)
    if contigs == 0:
        score = 0  # bölme hatasını önlemek için
    else:
        score = edges / (contigs * (dead_ends + 1)**2)

    return gfa_file, contigs, dead_ends, score