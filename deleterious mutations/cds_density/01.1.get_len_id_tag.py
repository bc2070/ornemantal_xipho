import sys
import os

def read_regions(filename):
    regions = []
    if not os.path.exists(filename):
        return regions
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
                regions.append((chrom, start, end))
            except ValueError:
                continue
    return regions

def read_gff(filename):
    genes = {} # {transcript_id: (chrom, start, end)}
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type == 'mRNA': 
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                attributes = parts[8]

                transcript_id = None
                for attr in attributes.split(';'):
                    if attr.startswith('ID=transcript:'):
                        transcript_id = attr.split(':')[1]
                        break
                    elif attr.startswith('ID='): 
                        transcript_id = attr.split('=')[1]
                        break

                if transcript_id:
                    genes[transcript_id] = (chrom, start, end)
    return genes

def read_cds_lengths(filename):

    cds_lengths = {} # {transcript_id: length}
    if not os.path.exists(filename):
        return cds_lengths
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                cds_lengths[parts[0]] = int(parts[1])
    return cds_lengths

def calculate_overlap(gene_start, gene_end, region_start, region_end):

    overlap_start = max(gene_start, region_start)
    overlap_end = min(gene_end, region_end)
    return max(0, overlap_end - overlap_start)

def main():
    if len(sys.argv) != 8:
        print("Usage: python script.py <gene_pos_gff> <cds_len_txt> <output_csv> <target_regions_dir> <inversion_bed> <N_bed> <low_posterior_bed>")
        sys.exit(1)

    gene_pos_gff = sys.argv[1]
    cds_len_txt = sys.argv[2]
    output_csv = sys.argv[3]
    target_regions_dir = sys.argv[4]
    inversion_bed = sys.argv[5]
    N_bed = sys.argv[6]
    low_posterior_bed = sys.argv[7]


    print(f"Reading gene positions from {gene_pos_gff}...")
    genes_info = read_gff(gene_pos_gff)
    print(f"Found {len(genes_info)} genes.")


    print(f"Reading CDS lengths from {cds_len_txt}...")
    cds_lengths = read_cds_lengths(cds_len_txt)
    print(f"Found {len(cds_lengths)} CDS lengths.")


    target_gene_labels = {} # {transcript_id: [label1, label2, ...]}
    
    if not os.path.isdir(target_regions_dir):
        print(f"Error: {target_regions_dir} is not a directory.")
        sys.exit(1)


    bed_files = sorted([f for f in os.listdir(target_regions_dir) if f.endswith('.bed')])
    
    if not bed_files:
        print(f"No .bed files found in {target_regions_dir}")

    for bed_file in bed_files:
        prefix = os.path.splitext(bed_file)[0] 
        full_path = os.path.join(target_regions_dir, bed_file)
        
        print(f"Processing target region file: {bed_file} (Label: {prefix})...")
        regions = read_regions(full_path)


        regions_by_chrom = {}
        for r_chrom, r_start, r_end in regions:
            if r_chrom not in regions_by_chrom:
                regions_by_chrom[r_chrom] = []
            regions_by_chrom[r_chrom].append((r_start, r_end))

        for transcript_id, (gene_chrom, gene_start, gene_end) in genes_info.items():
            if gene_chrom not in regions_by_chrom:
                continue
                
            gene_length = gene_end - gene_start + 1
            if gene_length <= 0:
                continue

            for r_start, r_end in regions_by_chrom[gene_chrom]:
                overlap = calculate_overlap(gene_start, gene_end, r_start, r_end)

                if overlap / gene_length >= 0.8:
                    if transcript_id not in target_gene_labels:
                        target_gene_labels[transcript_id] = []
                    if prefix not in target_gene_labels[transcript_id]:
                        target_gene_labels[transcript_id].append(prefix)
                    break


    print("Processing exclusion regions...")
    exclusion_beds = {
        "Inversion": inversion_bed,
        "N_Region": N_bed,
        "Low_Posterior": low_posterior_bed
    }
    excluded_genes = set()

    for exc_type, bed_file in exclusion_beds.items():
        if not os.path.exists(bed_file):
            print(f"Warning: Exclusion file {bed_file} ({exc_type}) not found. Skipping.")
            continue
        print(f"Processing exclusion file: {bed_file}...")
        regions = read_regions(bed_file)
        

        exc_by_chrom = {}
        for r_chrom, r_start, r_end in regions:
            if r_chrom not in exc_by_chrom:
                exc_by_chrom[r_chrom] = []
            exc_by_chrom[r_chrom].append((r_start, r_end))

        for transcript_id, (gene_chrom, gene_start, gene_end) in genes_info.items():
            if gene_chrom not in exc_by_chrom:
                continue
            gene_length = gene_end - gene_start + 1
            for r_start, r_end in exc_by_chrom[gene_chrom]:
                overlap = calculate_overlap(gene_start, gene_end, r_start, r_end)
                if overlap / gene_length >= 0.8:
                    excluded_genes.add(transcript_id)
                    break


    print("Generating final output...")
    final_results = []

    for transcript_id, labels in target_gene_labels.items():
        if transcript_id not in excluded_genes:
            cds_len = cds_lengths.get(transcript_id, 'N/A')

            sorted_labels = sorted(labels)
            final_results.append((transcript_id, ';'.join(sorted_labels), cds_len))


    with open(output_csv, 'w') as outfile:
        outfile.write("GeneID,Labels,CDS_Length\n")
        for gene_id, labels, cds_len in final_results:
            outfile.write(f"{gene_id},\"{labels}\",{cds_len}\n")

    print(f"Results saved to {output_csv}")
    print(f"Total genes in output: {len(final_results)}")

if __name__ == "__main__":
    main()
