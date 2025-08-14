import subprocess
import os
import glob
from Bio import SeqIO
import pandas as pd
import shutil

'''
Usage:  prodigal 
         -d:  Write nucleotide sequences of genes to the selected file.
         -i:  Specify FASTA/Genbank input file (default reads from stdin).
         -o:  Specify output file (default writes to stdout).
         -p:  Select procedure (single or meta).  Default is single.
         -q:  Run quietly (suppress normal stderr output).
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from datetime import date

def genbak_create(nuc_seq, seq_acc, seq_id, seq_description, feature_list, out_file_path): 

    feature_list.sort(key=lambda x: x.start)

    sequence_object = Seq(nuc_seq)


    record = SeqRecord(sequence_object, id=seq_id.split('-')[0],  annotations={"molecule_type": "DNA", "date":f"{date.day}-{date.month}-{date.year}"}, name=seq_id.split('-')[0], description=seq_description)
    
    for feature in feature_list:
        record.features.append(SeqFeature(FeatureLocation(start=feature.start, end=feature.end, strand=feature.strand), type="CDS",
            qualifiers={
                'gene' : feature.gene,
                'product' : feature.product,
                'res_type' : feature.r_type,
                'info' : feature.fullname,
                'score' : feature.score
            }
         ))
    
    output_file = open(out_file_path, 'w')
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

# cds is coding region list
class GeneticInfo:
    def __init__(self, seqacc, qseqid, seq_description, feature_list,  nuc_seq, score0, score1, score2):
        self.seq_acc = seqacc
        self.seq_id = qseqid
        self.seq_description = seq_description
        self.feature_list = feature_list
        self.nuc_seq = nuc_seq
        self.score0 = score0
        self.score1 = score1
        self.score2 = score2


def control_cycle_file(cycle_file_path):
    with open(cycle_file_path) as f:
        the_lines = f.readlines()
        len_of_lines = len(the_lines)
        if len_of_lines == 0:
            return False
        if 'Error' in the_lines[0]:
            return False
        if len_of_lines > 1:
            return True

def prodigal_driver(path_of_prodigal, cycle_file,  out_file_gbk, out_file_nuc):
    args = f"{path_of_prodigal} -i {cycle_file} -p meta -o {out_file_gbk} -d {out_file_nuc} -q" 
    my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True)
   

def split_fasta(in_file, out_folder):
    f_open = open(in_file, "r")
    for rec in SeqIO.parse(f_open, "fasta"):
        the_id = rec.id
        seq = rec.seq
        id_file = open(out_folder + '/' + the_id , "w")
        id_file.write(">"+str(the_id)+"\n"+str(seq))
        id_file.close()



def delete_blast_db(db_dir):
    try:
        shutil.rmtree(db_dir)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))

def make_blast_db(path_of_makeblastdb, db_input, db_output):
    nin_file = f"{db_output}.nin"  # nucleotide DB için
    if os.path.exists(nin_file):
        #print(f"[OK] BLAST DB zaten var: {db_output}")
        return
    
    args = f"{path_of_makeblastdb} -in {db_input} -dbtype nucl -out {db_output}"
    subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"[+] BLAST DB oluşturuldu: {db_output}")
def run_blast(path_of_blastn, query, database, output):
    extras = '"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"'
    args = f'{path_of_blastn} -db {database} -query {query} -out {output} -outfmt {extras}'
    my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True,  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def blast_driver(path_of_makeblastdb, path_of_blastn, out_blast_folder, db_path, blast_query, r_type, info_prod_dict, threshold_blast):
    cycle_file_name = os.path.basename(db_path)
    db_output = os.path.join(out_blast_folder,"blast_temp", cycle_file_name)
    db_dir = out_blast_folder + "/blast_temp"
    db_input = db_path
    result_path = os.path.join(out_blast_folder,"blast_files", os.path.basename(blast_query) + "_" + r_type + ".out")

    if not os.path.exists(db_dir):
        os.mkdir(db_dir)

    if not os.path.exists(out_blast_folder + "/blast_files"):
        os.mkdir(out_blast_folder + "/blast_files")

    if not os.path.exists(blast_query):
        print('No available Blast Query file')
        return False

    try:
        make_blast_db(path_of_makeblastdb, db_input, db_output)
        run_blast(path_of_blastn, blast_query, db_output, result_path)
    except Exception as e:
        print('Blast Error')
        raise UserWarning('Blast Error')
    
    cds_list = parsing_blast_file(result_path, r_type, threshold_blast, info_prod_dict)
    #if len(cds_list) == 0:
    #    print('no good result')
    #else:
    #    print('Good result')
   
    #delete_blast_db(db_dir)
    return cds_list


def parsing_blast_file(blast_result_file, r_type, threshold_blast, info_prod_dict, name_of_query=''):
    list_of_cds = []

    try:
        blast_result = pd.read_table(blast_result_file, header=None)
    except pd.errors.EmptyDataError:
        #print('empty blast file')
        return list_of_cds

    default_outfmt6_cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen'.strip().split(' ')
    blast_result.columns = default_outfmt6_cols
    df_filtered = blast_result[(blast_result['pident'] >= 80.0) & (blast_result['evalue'] < 1e10)]

    

    if len(df_filtered) > 0:
        if name_of_query == '':
            name_of_query = (df_filtered.iloc[0]['qseqid'])
        
        groups = df_filtered.groupby('qseqid')
        #print(groups['sseqid'])

        for qseqid, frame in groups:
            #print(frame.iloc[0]['sstart'])
            #print(len(frame))
            the_best_list = []
            for trv_in_group in range(0,len(frame)):
                
                the_matching_len = int(frame.iloc[trv_in_group]['length'])
                the_slen = int(frame.iloc[trv_in_group]['slen'])
                the_score_of_seq = (the_matching_len/the_slen)*float(frame.iloc[trv_in_group]['pident'])
                
                if threshold_blast < the_score_of_seq:
                    the_best_list.append((the_score_of_seq, trv_in_group))
                    
                    #higher than threshold

            the_sorted_list = sorted(the_best_list, key=lambda l:l[0], reverse=True)
            if len(the_sorted_list) != 0:
                the_score_of_seq = the_sorted_list[0][0]
                the_best_trv = the_sorted_list[0][1]
                the_start = int(frame.iloc[the_best_trv]['qstart'])
                the_end = int(frame.iloc[the_best_trv]['qend'])
                the_fullname = frame.iloc[the_best_trv]['sseqid']
                the_r_type = r_type

                if r_type != 'InsertionSequences':
                    the_start = the_start + int(info_prod_dict[frame.iloc[the_best_trv]['qseqid']][0])
                    the_end = the_end + int(info_prod_dict[frame.iloc[the_best_trv]['qseqid']][0])
                the_strand = '1'
                if int(the_start) > int(the_end):
                    the_strand = '-1'
                    temp_t = the_end
                    the_end = the_start
                    the_start = temp_t

                list_of_cds.append(CodingRegion(int(the_start), int(the_end), int(the_strand), the_fullname, the_r_type, the_score_of_seq))

    return list_of_cds




#Variables
#mean_of_CompTns = 5850
#std_of_CompTns = 2586
#total_score_type = 0
#threshold_final_score = 50
# Total score type can be 0, 1, 2 
#max_z = 20
#dist_type = 1
# dist type 1 related with negative z scores is 0
# dist type 0 is normal dist of z scores
# 1 is more accurate because lower lenghts means low number of genes, so score will be lower anyway
def scoring_main(cycle_folder, picota_out_folder, \
    path_to_antibiotics, path_to_xenobiotics, path_to_ises, \
    mean_of_CompTns = 5850, std_of_CompTns = 2586, total_score_type = 0, threshold_final_score = 50, \
    max_z = 20, dist_type = 1, path_of_prodigal = "prodigal", path_of_blastn = "blastn", path_of_makeblastdb = "makeblastdb", path_of_blastx = "blastx"):

    #FOLDERS
    picota_temp_folder = os.path.join(picota_out_folder, "Pico_Temp")
    splitted_cycle_folder = os.path.join(picota_temp_folder, "Splitted_Cycles")
    prodigal_out_folder = os.path.join(picota_temp_folder, "Prodigal_Out")
    out_blast_folder = os.path.join(picota_temp_folder, "Blast_Out")
    picota_final_tab = os.path.join(picota_out_folder, 'picota_final_tab')
    
    if not os.path.exists(picota_out_folder):
        os.mkdir(picota_out_folder)
    if not os.path.exists(picota_temp_folder):
        os.mkdir(picota_temp_folder)
    if not os.path.exists(splitted_cycle_folder):
        os.mkdir(splitted_cycle_folder)
    if not os.path.exists(prodigal_out_folder):
        os.mkdir(prodigal_out_folder)
    if not os.path.exists(out_blast_folder):
        os.mkdir(out_blast_folder)

    #INPUT CYCLE FILES
    final_list_comps = []

    cycle_files = glob.glob(cycle_folder+"/*")
    print('Number of cycle files: ',len(cycle_files))

    for cycle_file in cycle_files:


        if control_cycle_file(cycle_file):

            

            splitted_cycle_single_folder = splitted_cycle_folder + '/' + os.path.basename(cycle_file).split('.')[0]
            if not os.path.exists(splitted_cycle_single_folder):
                os.mkdir(splitted_cycle_single_folder)

            prodigal_out_for_cycle = prodigal_out_folder + '/' + os.path.basename(cycle_file).split('.')[0]
            if not os.path.exists(prodigal_out_for_cycle):
                os.mkdir(prodigal_out_for_cycle)

            picota_out_for_cycle = picota_out_folder +  '/' + os.path.basename(cycle_file).split('.')[0]
            if not os.path.exists(picota_out_for_cycle):
                os.mkdir(picota_out_for_cycle)

            genetic_info_list = []

            split_fasta(cycle_file, splitted_cycle_single_folder)
            splitted_cycles = glob.glob(splitted_cycle_single_folder+"/*")



            # For every cycle, blast for Antibiotics, Insertion Seqs, and Xenobiotics
            for splitted_cycle in splitted_cycles:

                print('.',end='',flush=True)
                cds_list = []

                out_file_gbk = prodigal_out_for_cycle + '/' + os.path.basename(splitted_cycle) + '.gbk'
                out_file_nuc = prodigal_out_for_cycle + '/' + os.path.basename(splitted_cycle) + '.fasta'
                try:
                    prodigal_driver(path_of_prodigal, splitted_cycle,  out_file_gbk, out_file_nuc)
                except:
                    raise Exception('Prodigal Error, control the program')
                else:
                    #print('Prodigal was finished, filen in:', out_file_nuc)
                    
                    info_prod_dict = {}
                    with open(out_file_nuc, 'r') as prod_fil_nuc_a:
                        prod_lines = prod_fil_nuc_a.readlines()
                    for prod_lin in prod_lines:
                        if '>' in prod_lin:
                            the_pr_id = prod_lin.replace('>','').split(' # ')[0]
                            the_pr_start = int(prod_lin.replace('>','').split(' # ')[1])
                            the_pr_end = int(prod_lin.replace('>','').split(' # ')[2])
                            info_prod_dict[the_pr_id] = (the_pr_start, the_pr_end)
                            
                            

                    
                    #Antibiotics
                    if os.path.exists(path_to_antibiotics):
                        r_type = 'Antibiotics'
                        cds_list.extend(blast_driver(path_of_makeblastdb, path_of_blastn, out_blast_folder, path_to_antibiotics, out_file_nuc, r_type, info_prod_dict, threshold_blast=80))
                    else:
                        print('No available Antibiotics DB, control the path')
                    #Xenobiotics
                    if os.path.exists(path_to_xenobiotics):
                        r_type = 'Xenobiotics'
                        cds_list.extend(blast_driver(path_of_makeblastdb, path_of_blastn, out_blast_folder, path_to_xenobiotics, out_file_nuc, r_type, info_prod_dict, threshold_blast=80))
                    else:
                        print('No available Xenobiotics DB, control the path')
                    #Insertion Sequences
                    if os.path.exists(path_to_ises):
                        r_type = 'InsertionSequences'
                        cds_list.extend(blast_driver(path_of_makeblastdb, path_of_blastn, out_blast_folder, path_to_ises, splitted_cycle, r_type, info_prod_dict, threshold_blast=80))
                    else:
                        print('No available InsertionSequences DB, control the path')

                    lst_ant = []
                    lst_is = []
                    lst_xe = []
                    for the_cds in cds_list:
                        
                    
                        if the_cds.r_type == 'Antibiotics':
                            lst_ant.append(the_cds.score)
                            try:
                                the_cds.product = the_cds.fullname.split('|')[-2]
                                the_cds.gene = the_cds.fullname.split('|')[-1]
                            except:
                                the_cds.product = the_cds.fullname
                                the_cds.gene = the_cds.fullname
                        if the_cds.r_type == 'Xenobiotics':
                            lst_xe.append(the_cds.score)
                            try:
                                the_cds.product = the_cds.fullname.split('|')[0].split(':')[0]
                                the_cds.gene = the_cds.fullname.split('|')[0].split(':')[1]
                            except:
                                the_cds.product = the_cds.fullname
                                the_cds.gene = the_cds.fullname
            
                        if the_cds.r_type == 'InsertionSequences':
                            lst_is.append(the_cds.score)
                            the_cds.product = the_cds.fullname
                    

                    with open(splitted_cycle, 'r') as sc_f:
                        the_cy_lines = sc_f.readlines()
                    len_of_cycle = len(the_cy_lines[1])
                    nuc_of_cycle = the_cy_lines[1]


                    score0 = calculate_total_score(0, int(dist_type), max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe)
                    score1 = calculate_total_score(1, int(dist_type), max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe)
                    score2 = calculate_total_score(2, int(dist_type), max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe)

                    t_score = score0
                    if int(total_score_type) == 0:
                        t_score = score0
                    elif int(total_score_type) == 1:
                        t_score = score1
                    elif int(total_score_type) == 2:
                        t_score = score2
                    else:
                        raise Exception('Error, total_score_type is no valid, it can oen of these: 0, 1, 2')

                    if t_score > threshold_final_score:
                        print(f'\nAnalyzing: {os.path.basename(splitted_cycle)}\n')
                        print(os.path.basename(splitted_cycle), ' (score0):', score0, ' (score1):', score1, ' (score2):', score2)
                        print('Antibiotics: ', len(lst_ant), 'Xenobiotics: ', len(lst_xe), 'ISes: ', len(lst_is))
                        print('--------------------------------------------------------------------')

                    qseqid = os.path.basename(splitted_cycle)
                    seqacc = 'No Accession'
                    seq_description = f'this annotation made by PICOTA pipeline, score0:{score0},  score1:{score1},  score2:{score2}'

                    the_gen_info = GeneticInfo(seqacc, qseqid, seq_description, cds_list,  nuc_of_cycle, score0, score1, score2)
                    
                   
                    
                    if threshold_final_score < t_score:
                        genetic_info_list.append((the_gen_info, t_score))
            
            sorted_genetic_info_list = sorted(genetic_info_list, key=lambda l:l[1], reverse=True)
            for gen_info in sorted_genetic_info_list:
                IS_str = []
                Ant_str = []
                Xeno_str = []
                for the_g_cds in gen_info[0].feature_list:
                    if the_g_cds.r_type == 'Antibiotics':
                        Ant_str.append(the_g_cds.product)
                    if the_g_cds.r_type == 'Xenobiotics':
                        Xeno_str.append(the_g_cds.product)
                    if the_g_cds.r_type == 'InsertionSequences':
                        IS_str.append(the_g_cds.product)

                #Write this functions
                # seq_id, SRR_acc, k-mer, score0, score1, score2
                final_list_comps.append('\t'.join((gen_info[0].seq_id, os.path.basename(cycle_file).split('.')[0].split('_')[0], \
                    os.path.basename(cycle_file).split('.')[0].split('_')[1], str(gen_info[0].score0), str(gen_info[0].score1), \
                        str(gen_info[0].score2), str(len(IS_str)), ';'.join(IS_str), str(len(Ant_str)),';'.join(Ant_str),\
                            str(len(Xeno_str)),';'.join(Xeno_str)))+'\n')
                picota_out_gbk_path = picota_out_for_cycle + '/' + str(int(gen_info[1])) + '_' + gen_info[0].seq_id + '.gbk'
                picota_out_fasta_path = picota_out_for_cycle + '/' + str(int(gen_info[1])) + '_' + gen_info[0].seq_id + '.fasta'
                genbak_create(gen_info[0].nuc_seq, gen_info[0].seq_acc, gen_info[0].seq_id, gen_info[0].seq_description, gen_info[0].feature_list, picota_out_gbk_path)

                pic_fasta = open(picota_out_fasta_path, 'w')
                pic_fasta.write(f'>{gen_info[0].seq_id}\n{gen_info[0].nuc_seq}')
                pic_fasta.close()

    pic_fin1 = open(picota_final_tab, 'w')
    pic_fin1.writelines(final_list_comps)
    pic_fin1.close()

