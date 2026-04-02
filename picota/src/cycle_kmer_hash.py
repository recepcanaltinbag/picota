from collections import defaultdict
import sys
import math

def reverse_complement(seq):
    complement_dict = {"A":"T", "T":"A", "G":"C", "C":"G", "a":"t", "t":"a", "g":"c", "c":"g"}
    return "".join([complement_dict[nucleotide] for nucleotide in reversed(seq)])


def get_kmer_hashes(seq, k):
    """Verilen bir sekans için k-mer'leri hashleyerek döner."""
    return set(seq[i:i + k] for i in range(len(seq) - k + 1))


# Progress bar fonksiyonu
def print_progress_bar(iteration, total, prefix='', suffix='', length=50, fill='█'):
    percent = f"{100 * (iteration / float(total)):.1f}"
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    sys.stdout.write(f'\r{prefix} |{bar}| {percent}% {suffix}')
    sys.stdout.flush()


# Yeni optimizasyonlu döngü
def filter_cycles_with_kmer_old(cycle_info_list, k_mer_sim, threshold_sim, name_prefix_cycle):
    
    
    kmer_hashes = []  # K-mer hashlerini saklamak için
    cycle_clear_list = []
    name_it = 0
    print(len(cycle_info_list))
    if len(cycle_info_list) == 0:
        print('Print empty Cycle List')
        return cycle_info_list
    print(cycle_info_list[0].name)


    for i, cycle_el in enumerate(cycle_info_list):
        reverse_ori = ''
        if cycle_el.reverseOriented:
            reverse_ori = 'reverseoriented_'
        
        # Mevcut cycle için k-mer hashlerini hesapla
        current_kmer_hashes = get_kmer_hashes(cycle_el.sequence, k_mer_sim)
        is_similar = False
        
        # Önceki tüm cycle'larla karşılaştırma
        for other_kmer_hashes in kmer_hashes:
            common_kmers = len(current_kmer_hashes.intersection(other_kmer_hashes))
            if (common_kmers / len(current_kmer_hashes)) * 100 >= threshold_sim:
                is_similar = True
                break

        if not is_similar:
            name_it += 1
            cycle_el.name = f"{name_prefix_cycle}_{reverse_ori}{name_it}"
            cycle_clear_list.append(cycle_el)
            kmer_hashes.append(current_kmer_hashes)  # Benzersiz olarak ekle

    # Progress bar'ı güncelle
        print_progress_bar(i, len(cycle_info_list), prefix='Processing:', suffix='Complete')
    
    print('\n')
    return cycle_clear_list


def filter_cycles_with_kmer(cycle_info_list, k_mer_sim, threshold_sim, name_prefix_cycle):
    """
    Benzer cycle'ları k-mer benzerliğine göre filtreler.

    Optimizasyon: inverted index (kmer → kabul edilen cycle indeksleri).
    Her yeni cycle sadece ortak k-mer paylaşan önceki cycle'larla karşılaştırılır.
    Ortalama karmaşıklık O(n) yerine O(n²) — büyük listelerde önemli fark.
    """
    if len(cycle_info_list) == 0:
        print('Print empty Cycle List')
        return cycle_info_list

    print(len(cycle_info_list))
    print(cycle_info_list[0].name)

    # accepted_combined[i] = i. kabul edilen cycle'ın (fwd | rc) k-mer seti
    accepted_combined = []
    # kmer_index: k-mer → [accepted cycle indeksleri]
    kmer_index = defaultdict(list)

    cycle_clear_list = []
    name_it = 0

    for i, cycle_el in enumerate(cycle_info_list):
        reverse_ori = 'reverseoriented_' if cycle_el.reverseOriented else ''

        seq = cycle_el.sequence
        rev_seq = reverse_complement(seq)

        fwd_kmers = get_kmer_hashes(seq, k_mer_sim)
        rc_kmers  = get_kmer_hashes(rev_seq, k_mer_sim)

        if not fwd_kmers:
            # Sekans k_mer_sim'den kısa — benzerlik hesaplanamaz, doğrudan ekle
            name_it += 1
            cycle_el.name = f"{name_prefix_cycle}_{reverse_ori}{name_it}"
            cycle_clear_list.append(cycle_el)
            print_progress_bar(i, len(cycle_info_list), prefix='Processing:', suffix='Complete')
            continue

        # Inverted index üzerinden aday cycle'ları bul
        all_query_kmers = fwd_kmers | rc_kmers
        candidate_indices = set()
        for km in all_query_kmers:
            candidate_indices.update(kmer_index[km])

        is_similar = False
        for cidx in candidate_indices:
            other = accepted_combined[cidx]
            common1 = len(fwd_kmers.intersection(other))
            common2 = len(rc_kmers.intersection(other))
            if (common1 / len(fwd_kmers)) * 100 >= threshold_sim or \
               (common2 / len(rc_kmers))  * 100 >= threshold_sim:
                is_similar = True
                break

        if not is_similar:
            name_it += 1
            cycle_el.name = f"{name_prefix_cycle}_{reverse_ori}{name_it}"
            cycle_clear_list.append(cycle_el)

            cidx = len(accepted_combined)
            combined = fwd_kmers | rc_kmers
            accepted_combined.append(combined)
            for km in combined:
                kmer_index[km].append(cidx)

        print_progress_bar(i, len(cycle_info_list), prefix='Processing:', suffix='Complete')

    print('\n')
    return cycle_clear_list
