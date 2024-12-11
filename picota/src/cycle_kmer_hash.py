from collections import defaultdict
import sys


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
def filter_cycles_with_kmer(cycle_info_list, k_mer_sim, threshold_sim, name_prefix_cycle):
    kmer_hashes = []  # K-mer hashlerini saklamak için
    cycle_clear_list = []
    name_it = 0
    print(len(cycle_info_list))
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