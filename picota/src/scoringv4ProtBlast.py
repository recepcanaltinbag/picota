import subprocess
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import date
import pandas as pd
import shutil
import logging
# --------------------------- Classes ---------------------------


logger: logging.Logger = None


class CodingRegion:
    def __init__(self, start, end, strand, fullname, r_type, score):
        self.start = start
        self.end = end
        self.strand = strand
        self.fullname = fullname
        self.r_type = r_type
        self.score = score
        self.product = ''
        self.gene = ''

class GeneticInfo:
    def __init__(self, seqacc, qseqid, seq_description, feature_list, nuc_seq, score0, score1, score2):
        self.seq_acc = seqacc
        self.seq_id = qseqid
        self.seq_description = seq_description
        self.feature_list = feature_list
        self.nuc_seq = nuc_seq
        self.score0 = score0
        self.score1 = score1
        self.score2 = score2

# --------------------------- Helper Functions ---------------------------

def genbak_create(nuc_seq, seq_acc, seq_id, seq_description, feature_list, out_file_path):
    feature_list.sort(key=lambda x: x.start)
    sequence_object = Seq(nuc_seq)
    today = date.today()
    record = SeqRecord(
        sequence_object,
        id=seq_id.split('-')[0],
        annotations={"molecule_type": "DNA", "date": f"{today.day}-{today.month}-{today.year}"},
        name=seq_id.split('-')[0],
        description=seq_description
    )

    for feature in feature_list:
        record.features.append(
            SeqFeature(
                FeatureLocation(start=feature.start, end=feature.end, strand=feature.strand),
                type="CDS",
                qualifiers={
                    'gene': feature.gene,
                    'product': feature.product,
                    'res_type': feature.r_type,
                    'info': feature.fullname,
                    'score': feature.score
                }
            )
        )

    with open(out_file_path, 'w') as output_file:
        SeqIO.write(record, output_file, 'genbank')


def calculate_total_score(total_score_type, dist_type, max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe):
    min_z = 0
    z = (abs(len_of_cycle - mean_of_CompTns))/std_of_CompTns
    if z > max_z:
        z = max_z
    if dist_type == 1:
        if len_of_cycle < mean_of_CompTns:
            z = 0


    z_c_l = 1 - (z - min_z)/(max_z - min_z)
    total_score = 0
    antc = 0
    isc = 0
    xc = 0
    if total_score_type == 0:
        for ant in lst_ant:
            antc += ant
        for ist in lst_is:
            isc += ist
        for xet in lst_xe:
            xc += xet
        total_score = (antc + isc + xc)**z_c_l

    elif total_score_type == 1:
        if len(lst_ant) > 0:
            antc = 100
        if len(lst_is) > 0:
            isc = 100
        if len(lst_xe) > 0:
            xc = 100
        total_score = (antc + isc + xc)**z_c_l

    elif total_score_type == 2:
        for ant in lst_ant:
            antc += ant
        if len(lst_is) > 0:
            isc = 1
        for xet in lst_xe:
            xc += xet
        if xc + antc == 0:
            antc = 50
        total_score = ((antc + xc)*isc)**z_c_l
    else:
        raise Exception('Error, total_score_type is no valid, it can one of these: 0, 1, 2')

    return total_score


def control_cycle_file(cycle_file_path):
    with open(cycle_file_path) as f:
        the_lines = f.readlines()
    if not the_lines or 'Error' in the_lines[0]:
        return False
    return True


def split_fasta(in_file, out_folder):
    os.makedirs(out_folder, exist_ok=True)
    with open(in_file, "r") as f_open:
        for rec in SeqIO.parse(f_open, "fasta"):
            id_file_path = os.path.join(out_folder, rec.id)
            with open(id_file_path, "w") as id_file:
                id_file.write(f">{rec.id}\n{rec.seq}")


def delete_blast_db(db_dir):
    try:
        shutil.rmtree(db_dir)
    except OSError as e:
        logger.error(f"Error: {e.filename} - {e.strerror}.")


# --------------------------- External Tools ---------------------------

def prodigal_driver(path_of_prodigal, cycle_file, out_file_gbk, out_file_nuc, out_file_prot):
    args = f"{path_of_prodigal} -i {cycle_file} -p meta -o {out_file_gbk} -d {out_file_nuc} -a {out_file_prot} -q"
    subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True)


def make_blast_db(path_of_makeblastdb, db_input, db_output, db_type="nucl"):
    ext = ".nin" if db_type == "nucl" else ".pin"
    check_file = f"{db_output}{ext}"
    if os.path.exists(check_file):
        return
    args = f"{path_of_makeblastdb} -in {db_input} -dbtype {db_type} -out {db_output}"
    subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    logger.info(f"[+] BLAST DB created: {db_output} ({db_type})")


def run_blast(path_of_blast, query, database, output):
    if not os.path.exists(output):
        extras = '"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"'
        args = f'{path_of_blast} -db {database} -query {query} -out {output} -outfmt {extras}'
        subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True,
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parsing_blast_file(blast_result_file, r_type, threshold_blast, info_prod_dict):
    list_of_cds = []

    try:
        blast_result = pd.read_table(blast_result_file, header=None)
    except pd.errors.EmptyDataError:
        return list_of_cds

    cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen'.split()
    blast_result.columns = cols
    df_filtered = blast_result[(blast_result['pident'] >= 80.0) & (blast_result['evalue'] < 1e10)]

    for qseqid, frame in df_filtered.groupby('qseqid'):
        the_best_list = []
        
        for idx in range(len(frame)):
            match_len = int(frame.iloc[idx]['length'])
            slen = int(frame.iloc[idx]['slen'])
            score = (match_len / slen) * float(frame.iloc[idx]['pident'])
            if score > threshold_blast:
                
                the_best_list.append((score, idx))
        if the_best_list:
            best_score, best_idx = sorted(the_best_list, key=lambda x: x[0], reverse=True)[0]
            start = int(frame.iloc[best_idx]['qstart'])
            end   = int(frame.iloc[best_idx]['qend'])
            fullname = frame.iloc[best_idx]['sseqid']

            # Protein hits -> nucleotide koordinatına çevir
            if r_type != 'InsertionSequences':
                offset = int(info_prod_dict[qseqid][0])
                start = (start - 1) * 3 + 1 + offset
                end   = end * 3 + offset

            strand = 1
            if start > end:
                strand = -1
                start, end = end, start

            list_of_cds.append(CodingRegion(start, end, strand, fullname, r_type, best_score))

    return list_of_cds


def blast_driver(path_of_makeblastdb, path_of_blast, out_blast_folder, db_path, blast_query,
                 r_type, info_prod_dict, threshold_blast, db_type="nucl"):

    cycle_file_name = os.path.basename(db_path)
    db_dir = os.path.join(out_blast_folder, "blast_temp")
    db_output = os.path.join(db_dir, cycle_file_name)
    db_input = db_path
    result_path = os.path.join(out_blast_folder, "blast_files", f"{os.path.basename(blast_query)}_{r_type}.out")

    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(os.path.join(out_blast_folder, "blast_files"), exist_ok=True)

    if not os.path.exists(blast_query):
        logger.warning('No available Blast Query file')
        return []

    try:
        make_blast_db(path_of_makeblastdb, db_input, db_output, db_type=db_type)
        run_blast(path_of_blast, blast_query, db_output, result_path)
    except Exception as e:
        raise UserWarning('Blast Error') from e

    return parsing_blast_file(result_path, r_type, threshold_blast, info_prod_dict)


# --------------------------- Main Scoring ---------------------------

def scoring_main(cycle_folder, picota_out_folder,
                 path_to_antibiotics, path_to_xenobiotics, path_to_ises,
                 mean_of_CompTns=5850, std_of_CompTns=2586,
                 total_score_type=0, threshold_final_score=50,
                 max_z=20, dist_type=1,
                 path_of_prodigal="prodigal",
                 path_of_blastn="blastn",
                 path_of_makeblastdb="makeblastdb",
                 path_of_blastx="blastx",
                 path_of_blastp="blastp", logger_name="picota_analysis"):

    global logger
    logger = logging.getLogger(logger_name)

    picota_temp_folder = os.path.join(picota_out_folder, "Pico_Temp")
    splitted_cycle_folder = os.path.join(picota_temp_folder, "Splitted_Cycles")
    prodigal_out_folder = os.path.join(picota_temp_folder, "Prodigal_Out")
    out_blast_folder = os.path.join(picota_temp_folder, "Blast_Out")
    picota_final_tab = os.path.join(picota_out_folder, 'picota_final_tab')
    for folder in [picota_out_folder, picota_temp_folder, splitted_cycle_folder, prodigal_out_folder, out_blast_folder]:
        os.makedirs(folder, exist_ok=True)

    # Cycle dosyaları
    if os.path.isfile(cycle_folder):
        cycle_files = [cycle_folder]
    elif os.path.isdir(cycle_folder):
        cycle_files = glob.glob(os.path.join(cycle_folder, "*.fasta"))
    else:
        cycle_files = []

    logger.info(f"Found Cycle Files: {cycle_files}")
    final_list_comps = []

    for cycle_file in cycle_files:
        if not control_cycle_file(cycle_file):
            continue

        splitted_cycle_single_folder = os.path.join(splitted_cycle_folder, os.path.basename(cycle_file).split('.')[0])
        prodigal_out_for_cycle = os.path.join(prodigal_out_folder, os.path.basename(cycle_file).split('.')[0])
        picota_out_for_cycle = os.path.join(picota_out_folder, os.path.basename(cycle_file).split('.')[0])
        for folder in [splitted_cycle_single_folder, prodigal_out_for_cycle, picota_out_for_cycle]:
            os.makedirs(folder, exist_ok=True)

        genetic_info_list = []

        # Cycle split
        split_fasta(cycle_file, splitted_cycle_single_folder)
        splitted_cycles = glob.glob(os.path.join(splitted_cycle_single_folder, "*"))

        for splitted_cycle in splitted_cycles:
            #print('.', end='', flush=True)
            cds_list = []

            # Prodigal outputs
            out_file_gbk  = os.path.join(prodigal_out_for_cycle, f'{os.path.basename(splitted_cycle)}.gbk')
            out_file_nuc  = os.path.join(prodigal_out_for_cycle, f'{os.path.basename(splitted_cycle)}.fasta')
            out_file_prot = os.path.join(prodigal_out_for_cycle, f'{os.path.basename(splitted_cycle)}.faa')

            prodigal_driver(path_of_prodigal, splitted_cycle, out_file_gbk, out_file_nuc, out_file_prot)

            # Prodigal gen başlangıçları
            info_prod_dict = {}
            with open(out_file_nuc, 'r') as f_nuc:
                for rec in SeqIO.parse(f_nuc, 'fasta'):
                    parts = rec.description.split(' # ')
                    info_prod_dict[rec.id] = (int(parts[1]), int(parts[2]))

            # BLAST
            if os.path.exists(path_to_antibiotics):
                cds_list.extend(blast_driver(path_of_makeblastdb, path_of_blastp, out_blast_folder,
                                             path_to_antibiotics, out_file_prot, 'Antibiotics', info_prod_dict,
                                             threshold_blast=50, db_type="prot"))
            if os.path.exists(path_to_xenobiotics):
                cds_list.extend(blast_driver(path_of_makeblastdb, path_of_blastp, out_blast_folder,
                                             path_to_xenobiotics, out_file_prot, 'Xenobiotics', info_prod_dict,
                                             threshold_blast=50, db_type="prot"))
            if os.path.exists(path_to_ises):
                cds_list.extend(blast_driver(path_of_makeblastdb, path_of_blastn, out_blast_folder,
                                             path_to_ises, splitted_cycle, 'InsertionSequences', info_prod_dict,
                                             threshold_blast=50, db_type="nucl"))
            
            # CDS score listeleri
            lst_ant, lst_is, lst_xe = [], [], []
            for the_cds in cds_list:
                if the_cds.r_type == 'Antibiotics':
                    lst_ant.append(the_cds.score)
                    try:
                        the_cds.product = the_cds.fullname.split('|')[-2]
                        the_cds.gene = the_cds.fullname.split('|')[-1]
                    except:
                        the_cds.product = the_cds.fullname
                        the_cds.gene = the_cds.fullname
                elif the_cds.r_type == 'Xenobiotics':
                    lst_xe.append(the_cds.score)
                    try:
                        the_cds.product = the_cds.fullname.split('|')[0].split(':')[0]
                        the_cds.gene = the_cds.fullname.split('|')[0].split(':')[1]
                    except:
                        the_cds.product = the_cds.fullname
                        the_cds.gene = the_cds.fullname
                elif the_cds.r_type == 'InsertionSequences':
                    lst_is.append(the_cds.score)
                    the_cds.product = the_cds.fullname

            # Len of cycle
            with open(splitted_cycle, 'r') as sc_f:
                the_cy_lines = sc_f.readlines()
            len_of_cycle = sum(len(line.strip()) for line in the_cy_lines if not line.startswith('>'))
            nuc_of_cycle = ''.join(line.strip() for line in the_cy_lines if not line.startswith('>'))


            # Score hesaplama
            score0 = calculate_total_score(0, dist_type, max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe)
            score1 = calculate_total_score(1, dist_type, max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe)
            score2 = calculate_total_score(2, dist_type, max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe)

            t_score = [score0, score1, score2][total_score_type]

            if t_score > threshold_final_score:
                logger.info(f'Analyzing: {os.path.basename(splitted_cycle)}')
                logger.info(
                    f"{os.path.basename(splitted_cycle)} (score0): {score0}, (score1): {score1}, (score2): {score2}"
                )
                logger.info(
                    f"Antibiotics: {len(lst_ant)}, Xenobiotics: {len(lst_xe)}, ISes: {len(lst_is)}"
                )
                logger.info('--------------------------------------------------------------------')

            # GeneticInfo objesi
            qseqid = os.path.basename(splitted_cycle)
            seqacc = 'No Accession'
            seq_description = f'this annotation made by PICOTA pipeline, score0:{score0}, score1:{score1}, score2:{score2}'
            the_gen_info = GeneticInfo(seqacc, qseqid, seq_description, cds_list, nuc_of_cycle, score0, score1, score2)

            if t_score > threshold_final_score:
                genetic_info_list.append((the_gen_info, t_score))

        # GeneticInfo sıralama ve final çıktılar
        sorted_genetic_info_list = sorted(genetic_info_list, key=lambda l: l[1], reverse=True)

        for gen_info in sorted_genetic_info_list:
            IS_str, IS_coords = [], []
            Ant_str, Ant_coords = [], []
            Xeno_str, Xeno_coords = [], []

            # önce IS'leri topla
            for the_g_cds in gen_info[0].feature_list:
                if the_g_cds.r_type == 'InsertionSequences':
                    IS_str.append(the_g_cds.product)
                    IS_coords.append(f"{the_g_cds.start}-{the_g_cds.end}")

            # sonra diğerlerini ekle
            for the_g_cds in gen_info[0].feature_list:
                start, end = the_g_cds.start, the_g_cds.end

                if the_g_cds.r_type == 'Antibiotics':
                    Ant_str.append(the_g_cds.product)
                    Ant_coords.append(f"{start}-{end}")

                elif the_g_cds.r_type == 'Xenobiotics':
                    # IS koordinatlarıyla çakışma kontrolü
                    overlap = any(
                        not (end < int(is_start) or start > int(is_end))
                        for is_start, is_end in (coord.split('-') for coord in IS_coords)
                    )
                    if not overlap:  # çakışmıyorsa ekle
                        Xeno_str.append(the_g_cds.product)
                        Xeno_coords.append(f"{start}-{end}")


            final_list_comps.append('\t'.join((
                gen_info[0].seq_id,
                os.path.basename(cycle_file).split('.')[0].split('_')[0],
                os.path.basename(cycle_file).split('.')[0].split('_')[1],
                str(gen_info[0].score0),
                str(gen_info[0].score1),
                str(gen_info[0].score2),
                str(len(IS_str)), ';'.join(IS_str), ';'.join(IS_coords),
                str(len(Ant_str)), ';'.join(Ant_str), ';'.join(Ant_coords),
                str(len(Xeno_str)), ';'.join(Xeno_str), ';'.join(Xeno_coords)
            )))

            # GenBank kaydı
            out_genbank_path = os.path.join(picota_out_for_cycle, f"{gen_info[0].seq_id}.gbk")
            genbak_create(gen_info[0].nuc_seq, gen_info[0].seq_acc, gen_info[0].seq_id,
                          gen_info[0].seq_description, gen_info[0].feature_list, out_genbank_path)

    # Final tab yazımı
    with open(picota_final_tab, 'w') as f_out:
        f_out.write('\t'.join(['CycleID', 'SRAID', 'kmer', 'score0', 'score1', 'score2',
                               'NumIS', 'ISproducts', 'IScoords',
                               'NumAnt', 'Antproducts', 'Antcoords',
                               'NumXeno', 'Xenoproducts', 'Xenocoords']) + '\n')
        for line in final_list_comps:
            f_out.write(line + '\n')
