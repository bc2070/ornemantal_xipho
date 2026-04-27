import sys
import os

def parse_gff(file_path):
    """Parse the GFF file to extract gene/transcript ID and sequences."""
    sequences = {}
    with open(file_path, 'r') as file:
        current_id = None
        current_seq = []
        for line in file:
            if line.startswith(">"):
                if current_id and current_seq:
                    sequences[current_id] = ''.join(current_seq)
                current_seq = []
                # Extract transcript and gene ID
                parts = line.split()
                current_id = parts[0].split(":")[1] if ":" in parts[0] else parts[0][1:]
            else:
                current_seq.append(line.strip())
        # Add the last sequence
        if current_id and current_seq:
            sequences[current_id] = ''.join(current_seq)
    return sequences

def extract_sequences(gff_files, gene_ids_file):
    """Extract sequences from GFF files based on gene IDs."""
    # Read gene IDs
    with open(gene_ids_file, 'r') as f:
        gene_ids = [line.strip() for line in f if line.strip()]

    for gene_id in gene_ids:
        output_file = f"{gene_id}.cds.fa"
        with open(output_file, 'w') as out_f:
            for gff_file in gff_files:
                gff_name = os.path.basename(gff_file).split('.')[0]
                sequences = parse_gff(gff_file)
                if gene_id in sequences:
                    seq = sequences[gene_id]
                    out_f.write(f">{gene_id}, {gff_name}\n{seq}\n")
                # Skip if gene ID not found

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <gene_ids_file> <gff_file1> <gff_file2> ...")
        sys.exit(1)

    gene_ids_file = sys.argv[1]
    gff_files = sys.argv[2:]
    extract_sequences(gff_files, gene_ids_file)

