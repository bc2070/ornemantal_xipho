import sys

def parse_id_from_gff(attr_str):
 
    parts = attr_str.strip().split(';')
    for part in parts:
        part = part.strip()
        if part.startswith('ID='):
            full_id = part.split('=', 1)[1]
            if ':' in full_id:
                return full_id.split(':', 1)[1]
            return full_id
    return None

def load_gene_lengths(gff_file):

    lengths = {}
    print(f"Loading gene lengths from {gff_file}...", file=sys.stderr)
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            

            gene_id = parse_id_from_gff(parts[8])
            if not gene_id:
                continue
            
            try:

                length = int(parts[4]) - int(parts[3])
                if length < 0: 
                    length = 0
                lengths[gene_id] = length
            except ValueError:
                continue
    return lengths

def load_dn_values(dn_file):

    dn_data = {}
    species_set = set()
    
    print(f"Loading dN values from {dn_file}...", file=sys.stderr)
    with open(dn_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            
            gene_id = parts[0]
            species = parts[1]
            val_str = parts[2]
            

            try:
                val = float(val_str)
            except ValueError:
                continue
                
            if gene_id not in dn_data:
                dn_data[gene_id] = {}
            
            dn_data[gene_id][species] = val
            species_set.add(species)
            
    return dn_data, sorted(list(species_set))

def main():
    if len(sys.argv) < 4:
        print("Usage: python calculate_d_value.py <region_gene_file> <dn_file> <gff_file> [output_file]")
        sys.exit(1)

    region_file = sys.argv[1]
    dn_file = sys.argv[2]
    gff_file = sys.argv[3]
    
    output_stream = sys.stdout
    if len(sys.argv) > 4:
        output_stream = open(sys.argv[4], 'w')


    gene_lengths = load_gene_lengths(gff_file)
    dn_data, species_list = load_dn_values(dn_file)

    print(f"Processing regions...", file=sys.stderr)
    

    header = ["chrom", "start", "end", "genes"] + [f"d_{sp}" for sp in species_list]
    output_stream.write("\t".join(header) + "\n")

    with open(region_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            
            # (chrom, start, end, gene_str)
            if len(parts) < 4:
                continue

            chrom = parts[0]
            try:
                r_start = int(parts[1])
                r_end = int(parts[2])
            except ValueError:
                continue
                
            gene_str = parts[3]
            

            fragment_len = r_end - r_start
            if fragment_len == 0:
                fragment_len = 1 

            genes = []
            if gene_str != ".":
                genes = gene_str.split(',')


            d_values = []
            
            for sp in species_list:
                weighted_sum = 0.0
                
                for gene_id in genes:

                    if gene_id in gene_lengths and \
                       gene_id in dn_data and \
                       sp in dn_data[gene_id]:
                        
                        length = gene_lengths[gene_id]
                        dn = dn_data[gene_id][sp]
                        
    
                        weighted_sum += (dn * length)
                

                d_val = weighted_sum / fragment_len

                d_values.append(f"{d_val:.8f}")


            out_line = parts[:4] + d_values
            output_stream.write("\t".join(map(str, out_line)) + "\n")

    if output_stream is not sys.stdout:
        output_stream.close()
    print("Done.", file=sys.stderr)

if __name__ == "__main__":
    main()
