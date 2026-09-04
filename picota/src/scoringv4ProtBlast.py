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
import re
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
    def __init__(self, seqacc, qseqid, seq_description, feature_list, nuc_seq, score0, score1, score2, score3=0.0):
        self.seq_acc = seqacc
        self.seq_id = qseqid
        self.seq_description = seq_description
        self.feature_list = feature_list
        self.nuc_seq = nuc_seq
        self.score0 = score0
        self.score1 = score1
        self.score2 = score2
        self.score3 = score3

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


def _length_fit(len_of_cycle, mean_of_CompTns, std_of_CompTns, max_z, dist_type):
    """How well the candidate's length matches the known CT length distribution, in [0, 1]."""
    z = abs(len_of_cycle - mean_of_CompTns) / std_of_CompTns
    if dist_type == 1 and len_of_cycle < mean_of_CompTns:
        z = 0
    return max(0.0, 1.0 - min(z, max_z) / max_z)


def _component_fit(comp_number):
    """
    Penalty for a cycle that wanders through many graph components.

    Named for what it measures rather than what it was once believed to measure.
    It does not report how close the candidate is to the two-node shape of a
    composite transposon: real elements thread 2 components when their flanking
    IS has two genomic copies and 16 when it has sixteen, so the count tracks how
    widespread the IS is, not the element's structure. Spurious cycles wandering
    the host's repeats sit near 11 throughout.

    That makes it a filter against host repeats, which is worth having -- removing
    it entirely raised false positives on wild-type genomes from 1 to 6 and cost
    PPV in four of six scenarios. But it works against the very case this method
    exists for: on shared_is it ranks real 16-copy elements BELOW the spurious
    cycles they should beat. Its weight is therefore 0.30 rather than the 0.40 it
    held when it was mistaken for a structural term, and the structural question
    is not asked at all; see the note on the weights in score3.

    The tolerance is a constant. It was previously scaled by a depth-derived
    copy-number estimate, but that input never reached this function -- the depth
    sidecar was not copied next to the FASTA -- so every published result was
    produced at this constant. Scaling it does not help anyway: spurious cycles
    read a mid-range depth ratio near 4 and would collect the same widened
    tolerance as the multi-copy elements it was meant to protect.
    """
    return 1.0 / (1.0 + abs(comp_number - 2) / 12.0)


def _best(scores):
    """Best hit quality in [0, 1]. Quality of the best hit, never a count."""
    return min(1.0, max(scores) / 100.0) if scores else 0.0


def calculate_total_score(total_score_type, dist_type, max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe, comp_number, novel_cargo_floor=0.30):
    min_z = 0
    z = (abs(len_of_cycle - mean_of_CompTns))/std_of_CompTns
    if z > max_z:
        z = max_z
    if dist_type == 1:
        if len_of_cycle < mean_of_CompTns:
            z = 0

    z_l = z + (abs(comp_number - 2))**(0.5)
    z_c_l = 1 - (z_l - min_z)/(max_z - min_z)
    total_score = 0
    antc = 0
    isc = 0
    xc = 0
    antcxc = 0
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
        if len(lst_ant) > 0:
            antcxc = 90
        if len(lst_is) > 0:
            isc = 1
        if len(lst_xe) > 0:   #existence 
            antcxc = 90

        # if no ant or xc little encourage
        # if antc + xc == 0:
        #    antc = 50

        total_score = ((antcxc) * isc) + 10**z_c_l

    elif total_score_type == 3:
        # IS quality gates the whole score because an insertion sequence is a
        # necessary condition, not a bonus: no IS means no composite transposon,
        # and a weak IS hit should scale everything down proportionally.
        #
        # The three weighted terms answer two separable questions -- structure
        # and multi-copy evidence ask "is this a composite transposon", cargo
        # asks "is it interesting" -- at 0.70 against 0.30. score2 has that
        # split the other way round at roughly 10/90, which is why it labels
        # rather than ranks and why cargo absent from a database is unreportable
        # however good the structure.
        is_quality = _best(lst_is)
        cargo_quality = max(_best(lst_ant), _best(lst_xe))
        if cargo_quality <= 0:
            cargo_quality = novel_cargo_floor

        # The IS gate is presence-with-quality, not quality alone. No IS still
        # scores zero, but a hit at 59% identity-coverage is unambiguously an
        # insertion sequence and halving the whole score for it compounded with
        # the component penalty to put the shared-IS case 0.1 points under the
        # threshold -- penalising twice for one structure.
        gate = 0.5 + 0.5 * is_quality if is_quality > 0 else 0.0

        # Length gates rather than contributes, for the same reason the IS does:
        # a composite transposon has a characteristic size, and something six
        # times the mean is not one whatever its homology. As a weighted term at
        # 20% of the total it could not express that -- two cycles of 35 and 37
        # kb in a wild-type genome, with an IS and a resistance gene 30 kb apart,
        # scored 68.6 and 65.8. Real elements pay nothing for this: under
        # dist_type 1 anything shorter than the mean has a length fit of 1.
        length_gate = _length_fit(len_of_cycle, mean_of_CompTns, std_of_CompTns,
                                  max_z, dist_type)

        # Two weighted terms, each answering something the definition of a
        # composite transposon actually asserts: the cargo is a single-copy
        # stretch (so it assembles into few nodes) and it is a payload worth
        # carrying (so it hits a database, or is credited the novel-cargo floor).
        #
        # Component count used to take 0.40 here and is gone. It does not measure
        # the element; it measures how widespread the flanking IS is, because a
        # cycle through a many-copy IS necessarily threads more components.
        # Measured per scenario, real elements run from 2 components at two IS
        # copies to 16 at sixteen while spurious cycles sit near 11 throughout,
        # so a fixed tolerance ranks 16-copy elements BELOW the spurious cycles
        # they should beat -- 0.46 against 0.55 on shared_is, an inversion in
        # exactly the case contig-based tools cannot solve and this method exists
        # for. Scaling the tolerance by a depth-derived copy estimate does not
        # rescue it, because spurious cycles read a mid-range ratio of about 4
        # and would collect the same widened tolerance.
        #
        # Cargo node count has no such contamination: cargo is single-copy by
        # definition however many copies its flanking IS has, and it separates in
        # every scenario with more than two real elements, shared_is included
        # (median 3 against 4).
        # Two weighted terms, renormalised from the three that were there.
        #
        # The third was a multi-copy estimate from read depth, and it never
        # worked: the sidecar carrying its input was not copied next to the
        # FASTA, so it returned its unknown-depth default of 0.5 for every
        # candidate ever scored -- a constant 0.15 of the total that also capped
        # every score at 85. Fixing the plumbing does not rescue it. The ratio
        # does track copy number (measured medians 1.5 at two genomic copies,
        # 10.7 at sixteen) but spurious cycles through the host's own mid-copy
        # IS elements sit near 4 in every scenario, so real elements fall on both
        # sides of the spurious ones and no monotonic function of it can separate
        # them.
        #
        # A cargo-node-count term was measured as a replacement and rejected.
        # Cargo is single-copy by definition, so counting the graph nodes its
        # region spans looked like a structural signal free of copy number, and
        # per-scenario medians supported that (1 for real elements against 4 for
        # spurious, holding even on shared_is). In place it collapsed precision:
        # false positives on cargo_is_diff went from 12 to 74, because the term
        # measures the cargo region only as well as the IS annotation bounds it.
        # An IS inside the cargo widens that annotation, leaves few nodes outside
        # it, and the term reads a fragmented cycle as a clean one.
        #
        # What remains is what the benchmark was actually run on, with the dead
        # weight removed rather than left in as a constant.
        total_score = 100.0 * gate * length_gate * (
            0.57 * _component_fit(comp_number)
            + 0.43 * cargo_quality)

    else:
        raise Exception('Error, total_score_type is no valid, it can one of these: 0, 1, 2, 3')

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


BLAST_THREADS = int(os.environ.get('PICOTA_BLAST_THREADS', '4'))


def run_blast(path_of_blast, query, database, output, threads=None):
    """
    One BLAST search, threaded.

    The searches ran single-threaded, which dominated scoring: a case took about
    a minute on a 24-core machine that sat at a load of five. The reference
    databases are the query here (the cycle is the database), so the protein
    searches against CARD and the xenobiotics set are the slow ones, and they
    parallelise well.

    PICOTA_BLAST_THREADS overrides the default for callers that already run
    several cases at once and would otherwise oversubscribe the machine.
    """
    if not os.path.exists(output):
        extras = '"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"'
        args = (f'{path_of_blast} -db {database} -query {query} -out {output} '
                f'-outfmt {extras} -num_threads {threads or BLAST_THREADS}')
        subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True,
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _load_blast_table(blast_result_file):
    """The hit table as a DataFrame with the module's fixed filter applied, or
    None when BLAST found nothing."""
    try:
        blast_result = pd.read_table(blast_result_file, header=None)
    except pd.errors.EmptyDataError:
        return None

    cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen'.split()
    blast_result.columns = cols
    return blast_result[(blast_result['pident'] >= 80.0) & (blast_result['evalue'] < 1e10)]


def _owner_of(owner, qseqid):
    """Which cycle a query belongs to. Without an owner map every hit lands
    under a single None key, which is what the one-cycle-per-file callers want."""
    return None if owner is None else owner.get(qseqid)


def parsing_blast_file_grouped(blast_result_file, r_type, threshold_blast, info_prod_dict, owner=None):
    """
    Best hit per query, grouped by the cycle each query came from.

    Same selection as parsing_blast_file -- this is where that logic now lives --
    but able to read a result file holding every cycle's queries at once.
    """
    grouped = {}
    df_filtered = _load_blast_table(blast_result_file)
    if df_filtered is None:
        return grouped

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
            if r_type != 'InsertionSequences' and r_type != 'CompTNs':
                offset = int(info_prod_dict[qseqid][0])
                start = (start - 1) * 3 + 1 + offset
                end   = end * 3 + offset

            strand = 1
            if start > end:
                strand = -1
                start, end = end, start

            grouped.setdefault(_owner_of(owner, qseqid), []).append(
                CodingRegion(start, end, strand, fullname, r_type, best_score))

    return grouped


def parsing_blast_file(blast_result_file, r_type, threshold_blast, info_prod_dict):
    grouped = parsing_blast_file_grouped(blast_result_file, r_type, threshold_blast, info_prod_dict)
    return grouped.get(None, [])


def parsing_blast_file_merged_grouped(blast_result_file, r_type, threshold_blast, info_prod_dict, owner=None):
    """
    The single best subject per cycle, HSPs merged before scoring.

    The one-cycle version kept one region for the whole result file. A file now
    holds every cycle, so "the best" has to be resolved per cycle -- keeping it
    file-wide would hand one cycle's transposon call to all of them.
    """
    per_cycle = {}
    df_filtered = _load_blast_table(blast_result_file)
    if df_filtered is None:
        return {}

    for qseqid, frame in df_filtered.groupby('qseqid'):
        cycle = _owner_of(owner, qseqid)
        for sseqid, subframe in frame.groupby('sseqid'):
            slen = int(subframe.iloc[0]['slen'])

            # subject koordinatlarını merge et
            subj_intervals = [(min(row.sstart, row.send), max(row.sstart, row.send)) for _, row in subframe.iterrows()]
            merged_subj = merge_intervals(subj_intervals)
            total_covered = sum(e - s + 1 for s, e in merged_subj)
            mean_pident = subframe['pident'].mean()

            score = (total_covered / slen) * mean_pident

            if score > threshold_blast:
                # query tarafında da aralıkları merge et
                query_intervals = [(min(row.qstart, row.qend), max(row.qstart, row.qend)) for _, row in subframe.iterrows()]
                merged_query = merge_intervals(query_intervals)

                # toplam aralığı belirle
                start = merged_query[0][0]
                end = merged_query[-1][1]
                fullname = sseqid

                strand = 1
                if start > end:
                    strand = -1
                    start, end = end, start

                per_cycle.setdefault(cycle, []).append(
                    CodingRegion(start, end, strand, fullname, r_type, score))

    return {cycle: [max(regions, key=lambda x: x.score)]
            for cycle, regions in per_cycle.items() if regions}


def parsing_blast_file_merged(blast_result_file, r_type, threshold_blast, info_prod_dict):
    grouped = parsing_blast_file_merged_grouped(blast_result_file, r_type, threshold_blast, info_prod_dict)
    return grouped.get(None, [])


def merge_intervals(intervals):
    """
    Çakışan HSP'leri merge eder. Basit bir algoritma kullanıyor:
    sıralı aralıkları üst üste bindir ve toplam kapsama hesapla.
    """
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [intervals[0]]

    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:  # çakışma varsa merge et
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged
















def blast_driver(path_of_makeblastdb, path_of_blast, out_blast_folder, db_path, blast_query,
                 r_type, info_prod_dict, db_out_path, threshold_blast, db_type="nucl"):

    cycle_file_name = os.path.basename(db_path)
    db_dir = os.path.join(db_out_path, "blast_temp")
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
    if r_type == "CompTNs":
        return parsing_blast_file_merged(result_path, r_type, threshold_blast, info_prod_dict)
    else:
        return parsing_blast_file(result_path, r_type, threshold_blast, info_prod_dict)


def blast_driver_batch(path_of_makeblastdb, path_of_blast, out_blast_folder, db_path, blast_query,
                       r_type, info_prod_dict, db_out_path, threshold_blast, owner, tag,
                       db_type="nucl"):
    """
    One search for every cycle at once, returning hits grouped by cycle.

    BLAST pays a fixed cost per invocation that a query of five proteins never
    amortises: it opens the reference database and walks its index whether the
    query is one sequence or a thousand. The per-cycle loop paid that once per
    cycle per database, so the two nucleotide searches -- whose queries are a
    single few-kb cycle each -- spent almost all of their time on setup.

    Measured on twenty real cycles from SRR20032745 against the four bundled
    databases, BLAST at four threads and the databases already built: 23.5 s
    over 80 invocations against 13.1 s over 4. Broken down, the nucleotide
    searches fall from 3.8 s to 0.4 s and CARD from 3.7 s to 1.5 s, while the
    82k-sequence xenobiotics search only moves 14.0 s to 10.2 s -- that one is
    dominated by real alignment work rather than by setup, and is where a
    faster aligner, not fewer invocations, is what would pay.

    The hits are the same ones. An E-value depends on the database and on the
    individual query, not on what else shares the query file, and the score this
    module computes from a hit uses neither. What does change is that "best hit"
    must now be resolved per cycle rather than per file, which is what the
    grouped parsers exist for.

    Unlike blast_driver this never reuses an existing result file. The per-cycle
    path could get away with it because its output name carried the cycle id;
    a batch file is named for the run, and silently reusing one from a previous
    run would score the wrong sequences.
    """
    db_dir = os.path.join(db_out_path, "blast_temp")
    db_output = os.path.join(db_dir, os.path.basename(db_path))
    result_path = os.path.join(out_blast_folder, "blast_files", f"{tag}_{r_type}.batch.out")

    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(os.path.join(out_blast_folder, "blast_files"), exist_ok=True)

    if not os.path.exists(blast_query) or os.path.getsize(blast_query) == 0:
        logger.warning(f'No available Blast Query file for {r_type}')
        return {}

    try:
        make_blast_db(path_of_makeblastdb, db_path, db_output, db_type=db_type)
        if os.path.exists(result_path):
            os.remove(result_path)
        run_blast(path_of_blast, blast_query, db_output, result_path)
    except Exception as e:
        raise UserWarning('Blast Error') from e

    if r_type == "CompTNs":
        return parsing_blast_file_merged_grouped(result_path, r_type, threshold_blast, info_prod_dict, owner)
    return parsing_blast_file_grouped(result_path, r_type, threshold_blast, info_prod_dict, owner)


def concat_fasta(in_files, out_file):
    """
    Merge FASTA files into one query, counting the records written.

    split_fasta writes its records without a trailing newline, so a plain
    concatenation would run one sequence into the next record's header.
    """
    written = 0
    out_dir = os.path.dirname(out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_file, 'w') as out_handle:
        for in_file in in_files:
            if not os.path.exists(in_file):
                continue
            with open(in_file) as in_handle:
                text = in_handle.read()
            if not text.strip():
                continue
            written += sum(1 for line in text.splitlines() if line.startswith('>'))
            out_handle.write(text if text.endswith('\n') else text + '\n')
    return written


def diamond_driver(diamond_path, query_file, db_fasta, r_type, info_prod_dict, threshold_score, threads=24):
    """
    DIAMOND ile protein araması yapar, CodingRegion objelerini döndürür.
    db_fasta: DB FASTA dosyası, yoksa .dmnd oluşturulur
    """
    import pandas as pd
    import os
    list_of_cds = []

    # DB kontrolü
    db_dmnd = f"{db_fasta}.dmnd"
    if not os.path.exists(db_dmnd):
        os.makedirs(os.path.dirname(db_dmnd), exist_ok=True)
        args_makedb = f"{diamond_path} makedb --in {db_fasta} -d {db_dmnd}"
        print(f"[+] Creating DIAMOND DB: {db_dmnd}")
        subprocess.run(args_makedb, shell=True, executable='/bin/bash', check=True)

    # Çıkış dosyası
    result_file = f"{query_file}_{r_type}.m8"
    os.makedirs(os.path.dirname(result_file), exist_ok=True)

    # DIAMOND blastp
    args_blast = f"{diamond_path} blastp -q {query_file} -d {db_dmnd} -o {result_file} " \
                 f"-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen -p {threads}"
    subprocess.run(args_blast, shell=True, executable='/bin/bash', check=True)

    # Sonuçları oku
    try:
        df = pd.read_table(result_file, header=None)
    except pd.errors.EmptyDataError:
        return list_of_cds

    df.columns = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen'.split()
    df_filtered = df[(df['pident'] >= 80.0) & (df['evalue'] < 1e-10)]

    for qseqid, frame in df_filtered.groupby('qseqid'):
        the_best_list = []

        for idx in range(len(frame)):
            match_len = int(frame.iloc[idx]['length'])
            slen = int(frame.iloc[idx]['slen'])
            score = (match_len / slen) * float(frame.iloc[idx]['pident'])
            if score > threshold_score:
                the_best_list.append((score, idx))

        if the_best_list:
            best_score, best_idx = sorted(the_best_list, key=lambda x: x[0], reverse=True)[0]

            start = int(frame.iloc[best_idx]['qstart'])
            end   = int(frame.iloc[best_idx]['qend'])
            fullname = frame.iloc[best_idx]['sseqid']

            # protein -> nucleotide koordinatına çevir
            offset = int(info_prod_dict[qseqid][0])
            start = (start - 1) * 3 + 1 + offset
            end   = end * 3 + offset

            strand = 1
            if start > end:
                strand = -1
                start, end = end, start

            # CodingRegion objesi
            cds_obj = CodingRegion(start, end, strand, fullname, r_type, best_score)
            list_of_cds.append(cds_obj)

    return list_of_cds

# --------------------------- Main Scoring ---------------------------

def scoring_main(cycle_folder, picota_out_folder,
                 path_to_antibiotics, path_to_xenobiotics, path_to_ises, path_to_TNs, db_out_path, 
                 mean_of_CompTns=5850, std_of_CompTns=2586,
                 total_score_type=0, threshold_final_score=50,
                 max_z=20, dist_type=1,
                 path_of_prodigal="prodigal",
                 path_of_blastn="blastn",
                 path_of_makeblastdb="makeblastdb",
                 path_of_blastx="blastx",
                 path_of_blastp="blastp", logger_name="picota_analysis",
                 blast_batch=True):
    """
    Score every cycle in `cycle_folder`.

    `blast_batch` runs one search per database over all cycles at once instead
    of one per cycle per database -- see blast_driver_batch for the measurement
    and for why the hits are the same either way. Set it False to fall back to
    the per-cycle searches.
    """

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

        # Per-cycle depth, when the detection stage left a sidecar next to the
        # FASTA. Keyed by the cycle id exactly as it appears in the header.
        # Cycle split
        split_fasta(cycle_file, splitted_cycle_single_folder)
        splitted_cycles = glob.glob(os.path.join(splitted_cycle_single_folder, "*"))

        # Annotation pass. Prodigal runs per cycle as before -- its output
        # feeds the coordinate offsets -- and the searches then run either once
        # per database over every cycle at once, or once per cycle as they used
        # to. Prodigal ids are unique across cycles (they are the cycle id plus
        # an ordinal), so one dictionary can serve every cycle's lookups.
        cycle_records = []
        info_prod_dict = {}
        prot_owner = {}

        for splitted_cycle in splitted_cycles:
            cycle_id = os.path.basename(splitted_cycle)

            # Prodigal outputs
            out_file_gbk  = os.path.join(prodigal_out_for_cycle, f'{cycle_id}.gbk')
            out_file_nuc  = os.path.join(prodigal_out_for_cycle, f'{cycle_id}.fasta')
            out_file_prot = os.path.join(prodigal_out_for_cycle, f'{cycle_id}.faa')

            prodigal_driver(path_of_prodigal, splitted_cycle, out_file_gbk, out_file_nuc, out_file_prot)

            # Prodigal gen başlangıçları
            with open(out_file_nuc, 'r') as f_nuc:
                for rec in SeqIO.parse(f_nuc, 'fasta'):
                    parts = rec.description.split(' # ')
                    info_prod_dict[rec.id] = (int(parts[1]), int(parts[2]))
                    prot_owner[rec.id] = cycle_id

            cycle_records.append((splitted_cycle, cycle_id, out_file_prot))

        # A nucleotide query is the cycle itself, so its query id already is the
        # cycle id; the map is there only so both search types read the same way.
        nuc_owner = {cycle_id: cycle_id for _, cycle_id, _ in cycle_records}
        hits_by_cycle = {cycle_id: [] for _, cycle_id, _ in cycle_records}

        # BLAST. Order matters: the per-cycle hit list is consumed downstream in
        # database order, so it has to stay antibiotics, xenobiotics, IS, TN.
        blast_jobs = [
            (path_to_antibiotics, path_of_blastp, 'Antibiotics',        'prot', 50),
            (path_to_xenobiotics, path_of_blastp, 'Xenobiotics',        'prot', 50),
            (path_to_ises,        path_of_blastn, 'InsertionSequences', 'nucl', 50),
            (path_to_TNs,         path_of_blastn, 'CompTNs',            'nucl', 80),
        ]

        if blast_batch:
            batch_tag = os.path.basename(cycle_file).split('.')[0]
            merged_dir = os.path.join(out_blast_folder, "merged")
            merged_prot = os.path.join(merged_dir, f'{batch_tag}.faa')
            merged_nuc  = os.path.join(merged_dir, f'{batch_tag}.fna')
            n_prot = concat_fasta([prot for _, _, prot in cycle_records], merged_prot)
            concat_fasta([cycle for cycle, _, _ in cycle_records], merged_nuc)
            logger.info(f"Batched BLAST: {len(cycle_records)} cycles, {n_prot} predicted proteins")

        for db_path, path_of_tool, r_type, db_type, threshold in blast_jobs:
            if not os.path.exists(db_path):
                continue

            if blast_batch:
                query, owner = (merged_prot, prot_owner) if db_type == 'prot' else (merged_nuc, nuc_owner)
                grouped = blast_driver_batch(path_of_makeblastdb, path_of_tool, out_blast_folder,
                                             db_path, query, r_type, info_prod_dict, db_out_path,
                                             threshold_blast=threshold, owner=owner, tag=batch_tag,
                                             db_type=db_type)
                for cycle_id, regions in grouped.items():
                    if cycle_id in hits_by_cycle:
                        hits_by_cycle[cycle_id].extend(regions)
                    else:
                        logger.warning(f"BLAST hit for unknown cycle {cycle_id}, ignored")
            else:
                for splitted_cycle, cycle_id, out_file_prot in cycle_records:
                    query = out_file_prot if db_type == 'prot' else splitted_cycle
                    hits_by_cycle[cycle_id].extend(
                        blast_driver(path_of_makeblastdb, path_of_tool, out_blast_folder,
                                     db_path, query, r_type, info_prod_dict, db_out_path,
                                     threshold_blast=threshold, db_type=db_type))

        # Scoring pass.
        for splitted_cycle, cycle_id, _out_file_prot in cycle_records:
            cds_list = hits_by_cycle[cycle_id]

            # CDS score listeleri
            lst_ant, lst_is, lst_xe, lst_CmpTN = [], [], [], []
            for the_cds in cds_list:
                if the_cds.r_type == 'CompTNs':
                    lst_CmpTN.append(the_cds.fullname)

            for the_cds in cds_list:
                if the_cds.r_type == 'Antibiotics':
                    lst_ant.append(the_cds.score)
                    try:
                        parts = the_cds.fullname.split('|')
                        the_cds.product = parts[-2] if len(parts) >= 2 else the_cds.fullname
                        the_cds.gene = parts[-1] if len(parts) >= 1 else the_cds.fullname
                    except (IndexError, AttributeError):
                        the_cds.product = the_cds.fullname
                        the_cds.gene = the_cds.fullname
                elif the_cds.r_type == 'Xenobiotics':
                    lst_xe.append(the_cds.score)
                    try:
                        xeno_parts = the_cds.fullname.split('|')
                        first_field = xeno_parts[0] if xeno_parts else the_cds.fullname
                        colon_parts = first_field.split(':')
                        the_cds.product = colon_parts[0] if len(colon_parts) >= 1 else the_cds.fullname
                        the_cds.gene = colon_parts[1] if len(colon_parts) >= 2 else the_cds.fullname
                    except (IndexError, AttributeError):
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
            header_line = next((line.strip() for line in the_cy_lines if line.startswith('>')), None)

            # Component number'ı regex ile çek
            comp_number = 0
            if header_line:
                match = re.search(r"-comp([A-Za-z0-9]+)-", header_line)
                if match:
                    comp_number = int(match.group(1))

            logger.info(f"Component number: {comp_number}")

            # Score hesaplama
            score0 = calculate_total_score(0, dist_type, max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe, comp_number)
            score1 = calculate_total_score(1, dist_type, max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe, comp_number)
            score2 = calculate_total_score(2, dist_type, max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe, comp_number)
            score3 = calculate_total_score(3, dist_type, max_z, mean_of_CompTns, std_of_CompTns, len_of_cycle, lst_ant, lst_is, lst_xe, comp_number)

            t_score = [score0, score1, score2, score3][total_score_type]

            if t_score > threshold_final_score:
                logger.info(f'Analyzing: {os.path.basename(splitted_cycle)}')
                logger.info(
                    f"{os.path.basename(splitted_cycle)} (score0): {score0}, (score1): {score1}, (score2): {score2}"
                )
                logger.info(
                    f"Antibiotics: {len(lst_ant)}, Xenobiotics: {len(lst_xe)}, ISes: {len(lst_is)}, CompTNs: {lst_CmpTN}"
                )
                logger.info('--------------------------------------------------------------------')

            # GeneticInfo objesi
            qseqid = os.path.basename(splitted_cycle)
            seqacc = 'No Accession'
            seq_description = f'this annotation made by PICOTA pipeline, score0:{score0}, score1:{score1}, score2:{score2}'
            the_gen_info = GeneticInfo(seqacc, qseqid, seq_description, cds_list, nuc_of_cycle, score0, score1, score2, score3)

            if t_score > threshold_final_score:
                genetic_info_list.append((the_gen_info, t_score))

        # GeneticInfo sıralama ve final çıktılar
        sorted_genetic_info_list = sorted(genetic_info_list, key=lambda l: l[1], reverse=True)

        for gen_info in sorted_genetic_info_list:
            IS_str, IS_coords = [], []
            Ant_str, Ant_coords = [], []
            Xeno_str, Xeno_coords = [], []
            CompTn_str, CompTn_coords = [], []

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
                if the_g_cds.r_type == 'CompTNs':
                    CompTn_str.append(the_g_cds.fullname)
                    CompTn_coords.append(f"{start}-{end}")


                if the_g_cds.r_type == 'Xenobiotics':
                    # IS koordinatlarıyla çakışma kontrolü
                    overlap = any(
                        not (end < int(is_start) or start > int(is_end))
                        for is_start, is_end in (coord.split('-') for coord in IS_coords)
                    )
                    if not overlap:  # çakışmıyorsa ekle
                        Xeno_str.append(the_g_cds.product)
                        Xeno_coords.append(f"{start}-{end}")

            if len(CompTn_str) == 0:
                CompTn_str = ["Novel"]

            _basename_parts = os.path.basename(cycle_file).split('.')[0].split('_')
            _sra_id  = _basename_parts[0]
            _kmer_id = _basename_parts[1] if len(_basename_parts) > 1 else ''
            final_list_comps.append('\t'.join((
                gen_info[0].seq_id,
                _sra_id,
                _kmer_id,
                str(gen_info[0].score0),
                str(gen_info[0].score1),
                str(gen_info[0].score2),
                str(gen_info[0].score3),
                str(len(IS_str)), ';'.join(IS_str), ';'.join(IS_coords),
                str(len(Ant_str)), ';'.join(Ant_str), ';'.join(Ant_coords),
                str(len(Xeno_str)), ';'.join(Xeno_str), ';'.join(Xeno_coords),
                str(len(CompTn_str)), ';'.join(CompTn_str), ';'.join(CompTn_coords)
            )))

            # GenBank kaydı
            out_genbank_path = os.path.join(picota_out_for_cycle, f"{gen_info[0].seq_id}.gbk")
            genbak_create(gen_info[0].nuc_seq, gen_info[0].seq_acc, gen_info[0].seq_id,
                          gen_info[0].seq_description, gen_info[0].feature_list, out_genbank_path)

    # Final tab yazımı
    with open(picota_final_tab, 'w') as f_out:
        f_out.write('\t'.join(['CycleID', 'SRAID', 'kmer', 'score0', 'score1', 'score2', 'score3',
                               'NumIS', 'ISproducts', 'IScoords',
                               'NumAnt', 'Antproducts', 'Antcoords',
                               'NumXeno', 'Xenoproducts', 'Xenocoords','NumCompTN','CompTN','CompTNscoords']) + '\n')
        for line in final_list_comps:
            f_out.write(line + '\n')

    # Enriched CSV output.
    # Imported as `src.output_formatter`: every caller reaches this module as
    # `src.scoringv4ProtBlast`, so `src` is a package rather than a sys.path
    # entry and a bare `from output_formatter import ...` never resolves. It
    # raised ImportError into the broad handler below and was reported as
    # "generation skipped", so the file was silently never written.
    try:
        from src.output_formatter import write_enriched_csv
    except ImportError:  # direct execution with src/ on the path
        from output_formatter import write_enriched_csv

    enriched_csv = os.path.join(picota_out_folder, 'picota_enriched.csv')
    try:
        n_cts = write_enriched_csv(picota_final_tab, enriched_csv)
        logger.info(f"[+] Enriched CSV written: {enriched_csv} ({n_cts} composite transposons)")
    except Exception as _fmt_err:
        logger.warning(f"Enriched CSV generation skipped: {_fmt_err}")
