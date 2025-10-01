import pysam
from collections import defaultdict
import os
import re


'''
def collect_annotated_blocks_single_readlen(bam_file, boundary=820, gap_tolerance=5, min_block_len=300, min_len_tol=0.1, len_tol=0.1):
    
    
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except OSError as e:
        print(f"\033[91m[HATA] {bam_file} açılırken sorun: {e}. Dosya atlanıyor...\033[0m")
        return {}
    
    read_blocks = defaultdict(lambda: defaultdict(list))
    ref_lengths = {ref: bam.get_reference_length(ref) for ref in bam.references}

    # 1) Tüm alignment bloklarını read bazında topla
    temp_blocks = defaultdict(lambda: defaultdict(list))
    read_lengths = defaultdict(int)  # her read için tek read_len

    for aln in bam:
        if aln.is_unmapped:
            continue
        qname = aln.query_name
        ref = aln.reference_name

        # read uzunluğu: query_length + hard clip baş ve son
        aln_read_len = aln.query_length
        if aln.cigartuples:
            # baştaki hard clip
            if aln.cigartuples[0][0] == 5:
                aln_read_len += aln.cigartuples[0][1]
            # sondaki hard clip
            if aln.cigartuples[-1][0] == 5:
                aln_read_len += aln.cigartuples[-1][1]

        # read_len dictionary’sinde max olanı al
        read_lengths[qname] = max(read_lengths[qname], aln_read_len)

        # read uzunluk filtresi
        ref_len = ref_lengths[ref]
        min_len = ref_len * (1 - min_len_tol)
        max_len = ref_len * (1 + len_tol)
        if not (min_len <= aln_read_len <= max_len):
            continue

        # hard clip offset
        q_offset = 0
        if aln.cigartuples:
            if aln.cigartuples[0][0] == 5:  # H
                q_offset += aln.cigartuples[0][1]

        pairs = aln.get_aligned_pairs(matches_only=True)
        if not pairs:
            continue

        q_start, r_start = pairs[0]
        prev_q, prev_r = pairs[0]
        blocks_temp = []

        for qpos, rpos in pairs[1:]:
            dq = qpos - prev_q
            dr = rpos - prev_r
            if dq > gap_tolerance or dr > gap_tolerance:
                block_len = prev_r - r_start + 1
                if block_len >= min_block_len:
                    blocks_temp.append((r_start, prev_r, q_start + q_offset, prev_q + q_offset))
                q_start, r_start = qpos, rpos
            prev_q, prev_r = qpos, rpos

        # son blok
        block_len = prev_r - r_start + 1
        if block_len >= min_block_len:
            blocks_temp.append((r_start, prev_r, q_start + q_offset, prev_q + q_offset))

        # boundary split ve annotation
        annotated_blocks = []
        for r1, r2, q1, q2 in blocks_temp:
            if r1 < boundary <= r2:
                len1 = boundary - r1
                if len1 >= min_block_len:
                    q_frac = int(q1 + (len1 * (q2 - q1) / (r2 - r1 + 1)))
                    annotated_blocks.append((r1, boundary-1, q1, q_frac, "transposon", aln_read_len))
                len2 = r2 - boundary + 1
                if len2 >= min_block_len:
                    q_frac2 = q_frac + 1 if len1 >= min_block_len else q1
                    annotated_blocks.append((boundary, r2, q_frac2, q2, "cargo", aln_read_len))
            else:
                annot = "transposon" if r1 < boundary else "cargo"
                annotated_blocks.append((r1, r2, q1, q2, annot, aln_read_len))

        temp_blocks[qname][ref].extend(annotated_blocks)

    bam.close()

    # 2) Read bazlı merge, overlap olursa ayrı bırak
    for qname, refs in temp_blocks.items():
        for ref, blist in refs.items():
            # tip ve ref_start’a göre sırala
            blist.sort(key=lambda x: (x[4], x[0]))
            merged = []
            for blk in blist:
                if not merged:
                    merged.append(blk)
                else:
                    prev = merged[-1]
                    # sadece aynı tip ve ardışık (prev[1]+1 == blk[0] ve prev[4]==blk[4]) ise merge
                    if blk[4] == prev[4] and blk[0] == prev[1] + 1:
                        merged[-1] = (
                            prev[0], max(prev[1], blk[1]), prev[2], max(prev[3], blk[3]), prev[4], read_lengths[qname]
                        )
                    else:
                        merged.append(blk)
            read_blocks[qname][ref].extend(merged)



    return read_blocks

'''


def collect_annotated_blocks_single_readlen(
    bam_file,
    gap_tolerance=5,
    min_block_len=300,
    min_len_tol=0.1,
    len_tol=0.1,
    merge_same_type=True,
    gap_merge=50  # read bazında merge için izin verilen boşluk
):
    """
    BAM içindeki alignments’ları read bazında toplayıp bloklara ayırır.
    Contig isimleri _transposon veya _cargo içeriyor.
    merge_same_type=True ise aynı tip blokları merge eder.
    gap_merge: read bazında ardışık bloklar arasındaki boşluk.
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except OSError as e:
        print(f"\033[91m[HATA] {bam_file} açılırken sorun: {e}. Dosya atlanıyor...\033[0m")
        return {}

    read_blocks = defaultdict(lambda: defaultdict(list))
    ref_lengths = {ref: bam.get_reference_length(ref) for ref in bam.references}
    temp_blocks = defaultdict(lambda: defaultdict(list))
    read_lengths = defaultdict(int)

    for aln in bam:
        if aln.is_unmapped:
            continue

        qname = aln.query_name
        ref = aln.reference_name

        # Contig tipi
        if re.search(r"_transposon(_\d+)?$", ref):
            annot = "transposon"
        elif ref.endswith("_cargo"):
            annot = "cargo"
        else:
            annot = "unknown"

        # read uzunluğu (hard clip dahil)
        aln_read_len = aln.query_length
        q_offset = 0
        if aln.cigartuples:
            if aln.cigartuples[0][0] == 5:
                aln_read_len += aln.cigartuples[0][1]
                q_offset += aln.cigartuples[0][1]
            if aln.cigartuples[-1][0] == 5:
                aln_read_len += aln.cigartuples[-1][1]

        read_lengths[qname] = max(read_lengths[qname], aln_read_len)

        # read uzunluk filtresi
        ref_len = ref_lengths[ref]
        min_len = ref_len * (1 - min_len_tol)
        max_len = ref_len * (1 + len_tol)
        if not (min_len <= aln_read_len <= max_len):
            continue

        # Alignment blokları
        pairs = aln.get_aligned_pairs(matches_only=True)
        if not pairs:
            continue

        q_start, r_start = pairs[0]
        prev_q, prev_r = pairs[0]
        blocks_temp = []

        for qpos, rpos in pairs[1:]:
            if qpos - prev_q > gap_tolerance or rpos - prev_r > gap_tolerance:
                block_len = prev_r - r_start + 1
                if block_len >= min_block_len:
                    blocks_temp.append((r_start, prev_r, q_start + q_offset, prev_q + q_offset, annot, aln_read_len))
                q_start, r_start = qpos, rpos
            prev_q, prev_r = qpos, rpos

        # son blok
        block_len = prev_r - r_start + 1
        if block_len >= min_block_len:
            blocks_temp.append((r_start, prev_r, q_start + q_offset, prev_q + q_offset, annot, aln_read_len))

        for blk in blocks_temp:
            temp_blocks[qname][annot].append(blk)

    bam.close()

    # Merge opsiyonel
    for qname, types in temp_blocks.items():
        for typ, blist in types.items():
            blist.sort(key=lambda x: x[2])  # q_start'a göre sırala

            if merge_same_type:
                merged = []
                for blk in blist:
                    if not merged:
                        merged.append(blk)
                    else:
                        prev = merged[-1]

                        if blk[4] == prev[4] and blk[2] <= prev[3] + gap_merge:
                            # Proximity merge: hangisi daha yakın, read mi ref mi?
                            read_gap = blk[2] - prev[3]
                            ref_gap  = blk[0] - prev[1]

                            if read_gap <= ref_gap:
                                # read bazlı merge
                                new_q_start = prev[2]
                                new_q_end   = max(prev[3], blk[3])
                                new_r_start = prev[0]
                                new_r_end   = max(prev[1], blk[1])
                            else:
                                # ref bazlı merge
                                new_q_start = prev[2]
                                new_q_end   = max(prev[3], blk[3])
                                new_r_start = prev[0]
                                new_r_end   = max(prev[1], blk[1])

                            merged[-1] = (new_r_start, new_r_end, new_q_start, new_q_end, prev[4], prev[5])
                        else:
                            merged.append(blk)
                read_blocks[qname][typ].extend(merged)
            else:
                read_blocks[qname][typ].extend(blist)

    return read_blocks


def write_cargo_transposon_reads_regex(blocks, out_file_ctc, out_file_partial):
    """
    Her read içindeki:
      - c->t->c patternlerini out_file_ctc dosyasına,
      - c->t veya t->c patternlerini out_file_partial dosyasına
    kaydeder.
    """
    with open(out_file_ctc, 'a') as f_ctc, open(out_file_partial, 'a') as f_partial:
        for read_id, refs in blocks.items():
            # Tüm blokları read bazında sırala
            all_blks = []
            for ref, blist in refs.items():
                all_blks.extend(blist)
            # q_start’a göre sırala
            all_blks.sort(key=lambda x: x[2])

            # annotation tiplerini al ve string oluştur
            pattern_str = ''.join(['c' if x[4] == 'cargo' else 't' for x in all_blks])

            # hangi dosyaya yazılacak
            if re.search(r'c.*t.*c', pattern_str):
                target_file = f_ctc
                pattern_name = 'c-t-c'
            elif re.search(r'c.*t|t.*c', pattern_str):
                target_file = f_partial
                pattern_name = 'partial'
            else:
                continue  # pattern yoksa atla

            # dosyaya yaz
            target_file.write(f">{read_id}\n")
            for r1, r2, q1, q2, annot, read_len in all_blks:
                target_file.write(f"{annot}: Ref {r1}-{r2} ({annot}) <=> Read {q1}-{q2} | ReadLen={read_len}\n")

            # terminal çıktısı
            print(f"\n{read_id} ({pattern_name})")
            for r1, r2, q1, q2, annot, read_len in all_blks:
                color = "\033[93m" if annot == "cargo" else "\033[91m"
                reset = "\033[0m"
                print(f"  {annot}: Ref {r1}-{r2} ({color}{annot}{reset}) <=> Read {q1}-{q2} | ReadLen={read_len}")



def bam_file_analyze(sorted_bam_file, out_file, out_file_partial):

    open(out_file, 'w').close()
    open(out_file_partial, 'w').close()
    len_tol = 5000      # %10 tolerans
    blocks = collect_annotated_blocks_single_readlen(
        sorted_bam_file,
        gap_tolerance=30,
        min_block_len=100,
        min_len_tol= 0.5,
        len_tol=len_tol
    )
    print(blocks)
    if not blocks:
        return None
    write_cargo_transposon_reads_regex(blocks, out_file, out_file_partial)


# === Kullanım ===


'''
if __name__ == "__main__":


    base_dir = '/media/lin-bio/back2/picota_IS26_test'
    out_file = '/media/lin-bio/back2/picota_IS26_test/annotated_reads_normalSRR22753363.txt'
    out_file_partial = '/media/lin-bio/back2/picota_IS26_test/annotated_reads_partialSRR22753363.txt'
    # dosyayı sıfırla
    open(out_file, 'w').close()
    # tüm BAM dosyalarını bul
    bam_files = []
    for root, dirs, files in os.walk(base_dir):
        for f in files:
            if f.endswith("SRR22753363_1_Tn4352-M20306_sorted.bam"): #change after demo
                bam_files.append(os.path.join(root, f))

    print(f"Found {len(bam_files)} BAM files.")
    for bam_file in bam_files:
        print('.', end='', flush=True)
        #print(f"\nProcessing: {bam_file}")
        #bam_file = '/media/lin-bio/back2/picota_IS26_test/circle_Control/Tn4352-M20306_ERR5902785_1/ERR5902785_1_Tn4352-M20306_sorted.bam'
        #boundary = 820      # transposon/cargo sınırı
        len_tol = 5000      # %10 tolerans
        blocks = collect_annotated_blocks_single_readlen(
            bam_file,
            gap_tolerance=30,
            min_block_len=100,
            min_len_tol= 0.5,
            len_tol=len_tol
        )
        print(blocks)
        # debug çıktı
    # debug çıktı (sadece tek tip blok olan readler)
    # debug çıktı (sadece hem transposon hem cargo olan readler)
        if not blocks:
            continue  # bozuk dosya ya da boş sonuç, atla
       
        for read, refs in blocks.items():
            for ref, blist in refs.items():
                block_types = set(annot for (_, _, _, _, annot, _) in blist)
                if "transposon" in block_types and "cargo" in block_types:
                    print(f"\n{read}")
                    for (r1, r2, q1, q2, annot, rlen) in blist:
                        print(f"  {ref}: Ref {r1}-{r2} ({annot}) <=> Read {q1}-{q2} | ReadLen={rlen}")
       
        # dosyaya yaz
        write_cargo_transposon_reads_regex(blocks, out_file, out_file_partial)
 '''