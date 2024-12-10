import re
import zstandard as zstd

def zst_file(zst_file_path, output_file_path):

    # Decompress and write the file
    with open(zst_file_path, "rb") as compressed_file:
        dctx = zstd.ZstdDecompressor()
        decompressed_data = dctx.decompress(compressed_file.read())

    with open(output_file_path, "wb") as output_file:
        output_file.write(decompressed_data)

    print(f"Decompressed data written to {output_file_path}")


def fasta_to_gfa(fasta_file, gfa_file):
    with open(fasta_file, "r") as fasta, open(gfa_file, "w") as gfa:
        # Write GFA header
        gfa.write("H\tVN:Z:1.0\n")

        node_id, sequence = None, ""
        for line in fasta:
            if line.startswith(">"):
                # Process previous node
                if node_id is not None:
                    # Write the segment with sequence
                    gfa.write(f"S\t{node_id}\t{sequence}\n")
                
                # Parse new header
                header = line.strip()
                match = re.match(r">SRR\d+_(\d+)", header)
                links = re.findall(r"L:([+-]):(\d+):([+-])", header)
                
                if match:
                    node_id = match.group(1)
                    sequence = ""  # Reset sequence buffer
                    # Write links (edges)
                    for direction, linked_node, linked_orientation in links:
                        gfa.write(f"L\t{node_id}\t{direction}\t{linked_node}\t{linked_orientation}\t0M\n")
            else:
                # Accumulate sequence
                sequence += line.strip()

        # Write the last node
        if node_id is not None:
            gfa.write(f"S\t{node_id}\t{sequence}\n")


# Paths for the input and output files
zst_file_path = "picota/test_data/SRR11362851.contigs.fa.zst"  # Input .zst file
fasta_path = "SRR643894.contigs.fa"  # Output uncompressed file
gfa_path = "assembly_graph_big.gfa"
# Example usage

zst_file(zst_file_path, fasta_path)
fasta_to_gfa(fasta_path, gfa_path)
print("GFA file created: assembly_graph.gfa")
