import sys


BIN_SIZE = 50000

def extract_id(attr_str):

    parts = attr_str.strip().split(';')
    for part in parts:
        part = part.strip()

        if part.startswith('ID='):

            full_id = part.split('=', 1)[1]
            

            if ':' in full_id:
                return full_id.split(':', 1)[1]
            else:
                return full_id
    return None

def load_genes_into_bins(gene_file):

    bins = {}
    print(f"Loading genes from {gene_file} ...", file=sys.stderr)
    
    with open(gene_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom = parts[0]

            try:
                g_start = int(parts[3]) - 1
                g_end = int(parts[4])
            except ValueError:
                continue
                
            g_len = g_end - g_start
            if g_len <= 0:
                continue
                

            gene_id = extract_id(parts[8])
            if not gene_id:
                continue 

            gene_tuple = (g_start, g_end, g_len, gene_id)


            if chrom not in bins:
                bins[chrom] = {}
            
            start_bin = g_start // BIN_SIZE
            end_bin = g_end // BIN_SIZE
            
            for b in range(start_bin, end_bin + 1):
                if b not in bins[chrom]:
                    bins[chrom][b] = []
                bins[chrom][b].append(gene_tuple)
                
    return bins

def calculate_overlap(r_start, r_end, g_start, g_end):

    overlap_start = max(r_start, g_start)
    overlap_end = min(r_end, g_end)
    return max(0, overlap_end - overlap_start)

def main():
    if len(sys.argv) < 3:
        print("Usage: python annotate_ids.py <gene_gff_file> <region_bed_file> [output_file]")
        sys.exit(1)

    gene_file = sys.argv[1]
    region_file = sys.argv[2]
    

    out_f = sys.stdout
    if len(sys.argv) > 3:
        out_f = open(sys.argv[3], 'w')


    gene_bins = load_genes_into_bins(gene_file)


    print(f"Processing regions from {region_file} ...", file=sys.stderr)
    
    with open(region_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            
            if len(parts) < 3:
                continue
                
            chrom = parts[0]
            try:
                r_start = int(parts[1])
                r_end = int(parts[2])
            except ValueError:
                continue

            found_ids = set()

            if chrom in gene_bins:
                start_bin = r_start // BIN_SIZE
                end_bin = r_end // BIN_SIZE
                
                candidate_genes = set()
                for b in range(start_bin, end_bin + 1):
                    if b in gene_bins[chrom]:
                        for gene in gene_bins[chrom][b]:
                            candidate_genes.add(gene)
                
                for g_start, g_end, g_len, g_id in candidate_genes:
                    overlap_len = calculate_overlap(r_start, r_end, g_start, g_end)
                    

                    if (overlap_len / g_len) > 0.8:
                        found_ids.add(g_id)


            id_str = ",".join(sorted(list(found_ids))) if found_ids else "."
            

            out_f.write(f"{line}\t{id_str}\n")

    if out_f is not sys.stdout:
        out_f.close()
        print("Done.", file=sys.stderr)

if __name__ == "__main__":
    main()
