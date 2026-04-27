import sys
import os
import shutil


BAD_CHARS = {'N', 'n', '*'}

def parse_fasta(file_path):
    """read FASTA，back (header, sequence)"""
    header = None
    seq_parts = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield header, "".join(seq_parts).upper()
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
        if header and seq_parts:
            yield header, "".join(seq_parts).upper()

def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_clean_phy.py file1.fasta file2.fasta ...")
        sys.exit(1)

    input_files = sys.argv[1:]
    file_names = []
    output_filenames = []
    output_handles = []

    print(f"Processing {len(input_files)} files...")
    print(f"Filtering out columns containing: {BAD_CHARS}")

    
    for f_path in input_files:
        base_name = os.path.basename(f_path)
        name_root, _ = os.path.splitext(base_name)
        file_names.append(name_root)
        out_name = f"{name_root}_clean.tmp"
        output_filenames.append(out_name)
        output_handles.append(open(out_name, 'w'))

    fasta_iterators = [parse_fasta(f) for f in input_files]

    total_sites_raw = 0
    kept_sites_final = 0

    
    while True:
        try:
            records = [next(it) for it in fasta_iterators]
        except StopIteration:
            break
        
        headers = [r[0] for r in records]
        seqs = [r[1] for r in records]
        
        
        current_chr_len = len(seqs[0])
        if any(len(s) != current_chr_len for s in seqs):
            print(f"[ERROR] Input FASTA files are not aligned at {headers[0]}!")
            print(f"Lengths detected: {[len(s) for s in seqs]}")
            sys.exit(1)

        total_sites_raw += current_chr_len
        
        valid_columns = []
        
        
        for bases_tuple in zip(*seqs):
            
            is_clean = True
            for b in bases_tuple:
                if b in BAD_CHARS:
                    is_clean = False
                    break
            
            if is_clean:
                valid_columns.append(bases_tuple)
        
        kept_count = len(valid_columns)
        kept_sites_final += kept_count
        
        ratio = (kept_count / current_chr_len * 100) if current_chr_len > 0 else 0
        print(f"Chr: {headers[0]} | Raw: {current_chr_len} | Kept: {kept_count} | Ratio: {ratio:.2f}%")

        if not valid_columns:
            continue
            
       
        filtered_seqs = zip(*valid_columns)
        for handle, filtered_seq_tuple in zip(output_handles, filtered_seqs):
            handle.write("".join(filtered_seq_tuple))

    
    for h in output_handles:
        h.close()

    global_ratio = (kept_sites_final / total_sites_raw * 100) if total_sites_raw > 0 else 0
    print("-" * 50)
    print(f"Global Raw Sites : {total_sites_raw}")
    print(f"Global Kept Sites: {kept_sites_final}")
    print(f"Global Keep Ratio: {global_ratio:.2f}%")
    print("-" * 50)

    
    phy_file = "all_samples_clean.phy"
    print(f"Generating final PHY file: {phy_file}")
    
    with open(phy_file, 'w') as phy_out:
        
        phy_out.write(f"{len(input_files)} {kept_sites_final}\n") 
        for name, tmp_file in zip(file_names, output_filenames):
            phy_out.write(f"{name}\t")
            with open(tmp_file, 'r') as f_in:
                shutil.copyfileobj(f_in, phy_out)
            phy_out.write("\n")

            os.remove(tmp_file)

    print("Done.")

if __name__ == "__main__":
    main()
