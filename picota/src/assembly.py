import subprocess
import glob
import os
import shutil
import logging
import shlex
from typing import List

logger: logging.Logger = None

# ----------------------------------------------------------------
# 1- FILTERING THE RAW FILE WITH FASTP
# ----------------------------------------------------------------
'''
fastp -i in.fq -o out.fq
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

'''
# file path is list
def raw_read_filtering(raw_file: List[str], out_folder: str, fastp_path: str, quiet_mode: bool) -> None:
    """Filter raw reads with fastp (supports single or paired-end)"""
    from pathlib import Path
    
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    if len(raw_file) == 1:
        raw_file_1_path = raw_file[0]
        filtered_raw_file_1_name = "filtered_" + Path(raw_file_1_path).stem
        f1_path = os.path.join(out_folder, filtered_raw_file_1_name + ".fastq")
        cmd = [fastp_path, '-i', raw_file_1_path, '-o', f1_path, '-h', f1_path + '.html']
        
    elif len(raw_file) == 2:
        raw_file_1_path = raw_file[0]
        raw_file_2_path = raw_file[1]
        filtered_raw_file_1_name = "filtered_" + Path(raw_file_1_path).stem
        filtered_raw_file_2_name = "filtered_" + Path(raw_file_2_path).stem
        f1_path = os.path.join(out_folder, filtered_raw_file_1_name + ".fastq")
        f2_path = os.path.join(out_folder, filtered_raw_file_2_name + ".fastq")
        cmd = [fastp_path, '-i', raw_file_1_path, '-I', raw_file_2_path, 
               '-o', f1_path, '-O', f2_path, '-h', f1_path + '.html']
    else:
        raise ValueError(f'Invalid number of input files: {len(raw_file)}. Expected 1 or 2.')
    
    logger.info(f'Running fastp: {" ".join(cmd)}')
    stdout_pipe = subprocess.DEVNULL if quiet_mode else None
    subprocess.run(cmd, check=True, stdout=stdout_pipe, stderr=subprocess.DEVNULL, text=True)



#k_mer_list = '39,59,79,99'
def _run_spades_assembly(spades_path: str, file_path: List[str], out_folder: str, 
                        gfa_folder: str, gfa_name: str, threads: int, k_mer: str, 
                        quiet_mode: bool, assembly_keep_temp_files: bool,
                        memory_gb: int = None, timeout_s: int = None,
                        delete_inputs: bool = True) -> None:
    """Helper function to run SPAdes assembly (reduces code duplication)

    memory_gb caps SPAdes' own allocator via -m. Without it SPAdes assumes
    250 GB, so on a smaller machine a high-coverage run does not fail -- it
    swaps, and the sample stops making progress instead of stopping. Capping
    below physical RAM turns that hang into an error the caller can record.

    timeout_s bounds the wall clock the same way: a sample that would run for
    days is killed and reported rather than blocking the queue behind it.
    """
    from pathlib import Path
    
    # Input validation
    for fpath in file_path:
        if not Path(fpath).exists():
            raise FileNotFoundError(f"Input file not found: {fpath}")
    
    # Build SPAdes command using list (prevents shell injection)
    cmd = [spades_path, '-o', out_folder, '-t', str(threads), '-k', k_mer]
    if memory_gb:
        cmd.extend(['-m', str(int(memory_gb))])
    
    if len(file_path) == 1:
        cmd.extend(['-1', file_path[0]])
    elif len(file_path) == 2:
        cmd.extend(['-1', file_path[0], '-2', file_path[1]])
    else:
        raise ValueError(f"Invalid number of input files: {len(file_path)}. Expected 1 or 2.")
    
    # Run assembly
    logger.info(f"Running SPAdes: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=quiet_mode, text=True, timeout=timeout_s)
    except subprocess.CalledProcessError as exc:
        # quiet_mode sends SPAdes' own diagnosis into the exception rather than
        # the terminal, and the caller sees only an exit status. Put the tail of
        # what SPAdes said into the message and the log, or a batch run records
        # a thousand identical "non-zero exit status" lines and no reason.
        tail = ""
        for stream in (exc.stderr, exc.stdout):
            if stream:
                tail = "\n".join(str(stream).strip().splitlines()[-15:])
                if tail:
                    break
        spades_log = os.path.join(out_folder, "spades.log")
        if not tail and os.path.exists(spades_log):
            with open(spades_log) as fh:
                tail = "\n".join(fh.read().splitlines()[-15:])
        logger.error(f"SPAdes exit {exc.returncode}. Last output:\n{tail}")
        raise RuntimeError(
            f"SPAdes exit {exc.returncode}: {tail.strip().splitlines()[-1] if tail.strip() else 'no output captured'}"
        ) from exc
    
    # Find and copy GFA file
    gfa_files = glob.glob(os.path.join(out_folder, '*.gfa'))
    if not gfa_files:
        raise RuntimeError('No GFA files generated during assembly. Check SPAdes output.')
    
    # Copy best GFA to destination
    shutil.copy(gfa_files[0], os.path.join(gfa_folder, gfa_name))
    logger.info(f'GFA file copied to {os.path.join(gfa_folder, gfa_name)}')
    
    # Cleanup temporary files if requested. The SPAdes work directory is this
    # k's alone and always goes; the input reads are shared by every k in
    # assembly_k_mer_list, so only the last pass may remove them.
    if not assembly_keep_temp_files:
        shutil.rmtree(out_folder)
        if delete_inputs:
            for file_pt in file_path:
                if Path(file_pt).exists():
                    Path(file_pt).unlink()
        logger.info('Temporary files deleted. Use --keep_temp_files to preserve.')


def assembly_driver_spades(spades_path, file_path, out_folder, gfa_folder, gfa_name, threads, k_mer, quiet_mode, assembly_keep_temp_files,
                           memory_gb=None, timeout_s=None, delete_inputs=True):
    """Main SPAdes driver - calls helper function"""
    _run_spades_assembly(spades_path, file_path, out_folder, gfa_folder, gfa_name, threads, k_mer, quiet_mode, assembly_keep_temp_files,
                         memory_gb=memory_gb, timeout_s=timeout_s, delete_inputs=delete_inputs)




def assembly_main(name_for_assembly, raw_file_list, main_out_folder, assembly_threads, assembly_k_mer_list,
        assembly_quiet, assembly_keep_temp_files,
        assembly_path_of_spades, assembly_path_of_fastp, assembly_skip_filtering, 
        assembler_type="spades", assembly_path_of_megahit=None, gfa_tools_path="gfatools", path_of_bandage="bandage", logger_name="picota_analysis",
        spades_memory_gb=None, spades_timeout_s=None):

    global logger
    logger = logging.getLogger(logger_name)

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
            print(f'GFA File exists, removing and recomputing: {gfa_folder}/{gfa_name}')
            os.remove(gfa_folder + '/' + gfa_name)

        best_gfa = assembly_driver_megahit(assembly_path_of_megahit, the_final_file_list, out_assembly, main_out_folder, gfa_folder, gfa_name, assembly_threads, assembly_quiet, assembly_keep_temp_files, gfa_tools_path, path_of_bandage)
        return best_gfa
    else:
        gfa_files = []
        k_list = [k.strip() for k in assembly_k_mer_list.split(',') if k.strip()]
        for k_index, k_mer_l in enumerate(k_list):
            gfa_name = k_mer_l + '.gfa'
            out_assembly = os.path.join(out_assembly_main, name_for_assembly + '_' + k_mer_l)
            if not os.path.exists(out_assembly):
                os.mkdir(out_assembly)
            gfa_files.append(gfa_folder + '/' + gfa_name)
            if os.path.exists(gfa_folder + '/' + gfa_name):
                print(f'GFA File exists, removing and recomputing: {gfa_folder}/{gfa_name}')
                os.remove(gfa_folder + '/' + gfa_name)

            assembly_driver_spades(assembly_path_of_spades, the_final_file_list, out_assembly, gfa_folder, \
                gfa_name, assembly_threads, k_mer_l, assembly_quiet, assembly_keep_temp_files,
                memory_gb=spades_memory_gb, timeout_s=spades_timeout_s,
                delete_inputs=(k_index == len(k_list) - 1))

        best_gfa = process_gfa_files(gfa_files, path_of_bandage)
        if best_gfa:
            # Dosyanın bulunduğu dizini al
            destination_path = os.path.join(main_out_folder, os.path.basename(best_gfa))
            shutil.copy(best_gfa, destination_path)

            print(f"Best GFA copied to: {destination_path}")
            return destination_path
        else:
            print("No best GFA file found to copy.")
            return None





#MEGAHIT


def assembly_driver_megahit(megahit_path, file_path, out_folder, main_out_folder, gfa_folder, gfa_name, threads, quiet_mode, assembly_keep_temp_files, gfa_tools_path, path_of_bandage):
    
    
    
    
    if len(file_path) == 1:
        args = f"{megahit_path} -r {file_path[0]} -o {out_folder} -t {str(threads)} --k-min 55"
    elif len(file_path) == 2:
        args = f"{megahit_path} -1 {file_path[0]} -2 {file_path[1]} -o {out_folder} -t {str(threads)} --k-min 55"
    else:
        logger.info('Error: there is no fastq file or more than two fastq file!')
        return

    logger.info('Command will be run:')
    logger.info(f"{args}")
    logger.info('-------')
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
            logger.info(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)

            cmd = f"{gfa_tools_path} {fastg_path} > {gfa_path}"
            logger.info(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)

        # Skoru GFA'dan hesapla
        gfa_files.append(gfa_path)


    best_gfa = process_gfa_files(gfa_files, path_of_bandage)

    if best_gfa:
        # Dosyanın bulunduğu dizini al
        destination_path = os.path.join(main_out_folder, os.path.basename(best_gfa))
        shutil.copy(best_gfa, destination_path)

        logger.info(f"Best GFA copied to: {destination_path}")
        result_gfa = destination_path
    else:
        logger.info("No best GFA file found to copy.")
        result_gfa = None


    if assembly_keep_temp_files == False:
        shutil.rmtree(out_folder)
        for file_pt in file_path:
            if os.path.exists(file_pt):
                os.remove(file_pt)
        logger.info('Temp Files deleted., if you want to keep them use --keep_temp_files')

    return result_gfa





def process_gfa_files(gfa_files, path_of_bandage):
    scores = []

    # Cache path (aynı dizinde gfa_scores.json / metadata sakla)
    cache_path = None
    if gfa_files:
        cache_path = os.path.join(os.path.dirname(os.path.abspath(gfa_files[0])), 'gfa_scores.json')

    # Eğer cache aynı gfa setine aitse tekrar hesaplama yapma
    if cache_path and os.path.exists(cache_path):
        try:
            with open(cache_path, 'r') as f:
                cached = json.load(f)
            cached_paths = [os.path.abspath(p) for p in cached.keys()]
            input_paths = [os.path.abspath(p) for p in gfa_files]
            if set(cached_paths) == set(input_paths):
                for path in gfa_files:
                    data = cached.get(os.path.abspath(path)) or cached.get(path)
                    if data is None:
                        raise KeyError
                    scores.append((path, data['contigs'], data['dead_ends'], data['score']))
                if logger is not None:
                    logger.info('Reusing cached GFA scores from %s', cache_path)
                else:
                    print('Reusing cached GFA scores from', cache_path)
        except Exception:
            scores = []  # bozulmuş cache durumunda yeniden hesaplama yapacağız

    for gfa_file in gfa_files:
        if os.path.getsize(gfa_file) == 0:
            if logger is not None:
                logger.info(f"Skipping empty GFA: {gfa_file}")
            else:
                print(f"Skipping empty GFA: {gfa_file}")
            continue

        # Eğer zaten cache kullanıldıysa ve gfa dosyaları bulunduysa skor eklemeyebiliriz.
        if scores and any(os.path.abspath(t[0]) == os.path.abspath(gfa_file) for t in scores):
            continue

        gfa_path, contigs, dead_ends, score = compute_score_from_gfa(gfa_file, path_of_bandage)
        scores.append((gfa_path, contigs, dead_ends, score))

    if not scores:
        if logger is not None:
            logger.info("No valid GFA files found.")
        else:
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
    if logger is not None:
        logger.info("\nGFA Scores (Normalized):")
    print("\nGFA Scores (Normalized):")
    for path, contigs, dead_ends, norm_score in normalized_scores:
        msg = f"{os.path.basename(path)} | contigs: {contigs}, dead ends: {dead_ends}, normalized score: {norm_score:.4f}"
        if logger is not None:
            logger.info(msg)
        print(msg)

    # Cache yaz
    if cache_path:
        cached_out = {}
        for path, contigs, dead_ends, score in scores:
            cached_out[os.path.abspath(path)] = {
                'contigs': contigs,
                'dead_ends': dead_ends,
                'score': score
            }
        try:
            with open(cache_path, 'w') as f:
                json.dump(cached_out, f, indent=2)
            if logger is not None:
                logger.info(f"Saved GFA score cache to {cache_path}")
            else:
                print(f"Saved GFA score cache to {cache_path}")
        except Exception:
            pass

    # En yüksek normalize skora sahip GFA dosyasını bul
    best_gfa = max(normalized_scores, key=lambda x: x[3])
    best_msg = f"\nBest GFA: {os.path.basename(best_gfa[0])} with score {best_gfa[3]:.4f}"
    if logger is not None:
        logger.info(best_msg)
    print(best_msg)
    return best_gfa[0]





def compute_score_from_gfa(gfa_file, path_of_bandage):
    # Bandage ile dead ends sayısını al
    if logger is not None:
        logger.info(f'gfa_process: {gfa_file}')
    else:
        print(f'gfa_process: {gfa_file}')
    cmd = f'{shlex.quote(path_of_bandage)} info {shlex.quote(gfa_file)} | grep "Dead ends" | grep -oP "\\d+"'
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