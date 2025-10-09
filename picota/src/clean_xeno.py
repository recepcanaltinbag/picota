from collections import Counter
from Bio import SeqIO

# ---- 1. Core enzymatic genes ----
enzymatic_core = [
    "hydrogenase", "oxygenase", "reductase", "dehydrogenase", "dioxygenase", "rieske", "ferrodoxin", "transporter", "regulator",
    "kinase", "transferase", "ligase", "lyase", "isomerase", "mutase"
]

# ---- 2. Xenobiotic / biodegradation enzymes----
xenobiotic_related = [
    # Halojenli bileşik bozunması
    "haloalkane dehalogenase",
    "haloacid dehalogenase",
    "chlorocatechol dioxygenase",
    "chlorobenzoate dioxygenase",

    # Aromatik halka kırıcılar
    "catechol 1,2-dioxygenase",
    "catechol 2,3-dioxygenase",
    "protocatechuate dioxygenase",
    "gentisate dioxygenase",
    "salicylate hydroxylase",
    "vanillate demethylase",

    # Pestisit / herbisit metabolizması
    "atrazine chlorohydrolase",
    "carbaryl hydrolase",
    "organophosphorus hydrolase",
    "parathion hydrolase",

    # Nitroaromatik ve patlayıcı bozunması
    "nitrobenzene nitroreductase",
    "dinitrotoluene dioxygenase",
    "trinitrotoluene reductase",

    # Genel detox enzimleri
    "aldo-keto reductase",
    "alcohol dehydrogenase",
    "formaldehyde dehydrogenase",
    "glyoxalase",
    "epoxide hydrolase",

    # Antioksidan / koruyucu
    "glutathione peroxidase",
    "thioredoxin reductase",
    "multicopper oxidase",
    "laccase"
]

# ---- 3. Hariç tutulanlar ----
exclude_keywords = [
    "transposase", "integrase", "recombinase", "resolvase"
]

# ---- 4. Input / Output ----
fasta_file = "picota/DBs/Xenobiotics/Xenobiotics.fasta"           # input fasta
out_file = "picota/DBs/Xenobiotics/Xenobiotics_classified.fasta"  # output fasta

enzymatic, xenobiotic, excluded, other = [], [], [], []

# ---- 5. Sequence sınıflandırma ----
for record in SeqIO.parse(fasta_file, "fasta"):
    header = record.description.lower()

    if any(excl in header for excl in exclude_keywords):
        excluded.append(record)
    elif any(kw in header for kw in enzymatic_core):
        enzymatic.append(record)
    elif any(kw in header for kw in xenobiotic_related):
        xenobiotic.append(record)
    else:
        other.append(record)

# ---- 6. Özet istatistikler ----
counts = {
    "enzymatic_core": len(enzymatic),
    "xenobiotic_related": len(xenobiotic),
    "excluded": len(excluded),
    "other": len(other)
}

enzymatic_subcounts = Counter(
    kw for rec in enzymatic for kw in enzymatic_core if kw in rec.description.lower()
)
xenobiotic_subcounts = Counter(
    kw for rec in xenobiotic for kw in xenobiotic_related if kw in rec.description.lower()
)

# ---- 7. Terminal çıktısı ----
print("[SUMMARY]")
print(f"enzymatic_core: {counts['enzymatic_core']}")
for k, v in enzymatic_subcounts.items():
    print(f"  - {k}: {v}")

print(f"xenobiotic_related: {counts['xenobiotic_related']}")
for k, v in xenobiotic_subcounts.items():
    print(f"  - {k}: {v}")

print(f"excluded: {counts['excluded']}")
print(f"other: {counts['other']}")

# ---- 8. FASTA çıktısı (sıralı) ----
with open(out_file, "w") as out:
    for rec in enzymatic + xenobiotic + excluded:
        out.write(f">{rec.description}\n{str(rec.seq)}\n")

