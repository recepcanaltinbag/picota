from Bio import Entrez, SeqIO
import os
import time

Entrez.email = "senin.emailin@example.com"  # Buraya mailini yaz

output_folder = "bacteria_xeno_proteins"
os.makedirs(output_folder, exist_ok=True)

search_terms = [
    "xenobiotic AND bacteria[Organism]",
    "cytochrome P450 AND bacteria[Organism]",
    "dehalogenase AND bacteria[Organism]",
    "monooxygenase AND bacteria[Organism]",
    "dioxygenase AND bacteria[Organism]"
]

BATCH_SIZE = 200  # Batch halinde çekilecek kayıt sayısı

def fetch_protein_ids(term):
    print(f"Aranıyor: {term}")
    handle = Entrez.esearch(db="protein", term=term, retmax=100000)  # Çok büyük limit
    record = Entrez.read(handle)
    handle.close()
    print(f" - Toplam {len(record['IdList'])} protein bulundu.")
    return record['IdList']

def fetch_and_save_protein_seqs_batch(protein_ids):
    for i in range(0, len(protein_ids), BATCH_SIZE):
        batch = protein_ids[i:i+BATCH_SIZE]
        print(f"Batch indiriliyor: {i+1} - {i+len(batch)} / {len(protein_ids)}")
        ids_str = ",".join(batch)
        try:
            handle = Entrez.efetch(db="protein", id=ids_str, rettype="fasta", retmode="text")
            records = SeqIO.parse(handle, "fasta")
            for rec in records:
                filename = f"{rec.id}.fasta"
                filepath = os.path.join(output_folder, filename)
                with open(filepath, "w") as f:
                    f.write(rec.format("fasta"))
                print(f"  {filename} kaydedildi.")
            handle.close()
        except Exception as e:
            print(f"Batch indirme hatası: {e}")
        time.sleep(1)  # NCBI isteği limiti için yavaşlat

all_protein_ids = set()

for term in search_terms:
    ids = fetch_protein_ids(term)
    all_protein_ids.update(ids)
    time.sleep(1)  # NCBI limitleri için bekle

print(f"\nToplam benzersiz protein sayısı: {len(all_protein_ids)}")

fetch_and_save_protein_seqs_batch(list(all_protein_ids))

print("\nTüm protein dizileri indirildi ve kaydedildi.")
