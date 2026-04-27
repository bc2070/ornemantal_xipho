import sys
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_phy.py input.phy")
        sys.exit(1)

    input_file = sys.argv[1]
    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_fixed{ext}"

    sequences = []
    
    print(f"Reading {input_file}...")

    with open(input_file, 'r') as f:
        first_line = f.readline()
        

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue
            
            parts = line.split(maxsplit=1)
            
            if len(parts) < 2:
                print(f"Warning: Line {line_num} seems malformed (no sequence found). Skipping.")
                continue
            
            name = parts[0]
            seq = parts[1]
            
            seq = seq.replace(" ", "").replace("\t", "")
            
            sequences.append((name, seq))

    if not sequences:
        print("Error: No sequences found in the file.")
        sys.exit(1)

    num_seqs = len(sequences)
    
    ref_len = len(sequences[0][1])
    
    print(f"Detected: {num_seqs} sequences.")
    print(f"Base length (from first seq): {ref_len}")

    for name, seq in sequences:
        if len(seq) != ref_len:
            print(f"\n[ERROR] Sequence length mismatch detected!")
            print(f"Sample '{sequences[0][0]}' length: {ref_len}")
            print(f"Sample '{name}' length: {len(seq)}")
            print("Phylip (.phy) files required all sequences to be of equal length.")
            print("Please check your data processing step (e.g. removal of *).")
            sys.exit(1)

    with open(output_file, 'w') as f_out:
        f_out.write(f"{num_seqs} {ref_len}\n")
        
        for name, seq in sequences:
            f_out.write(f"{name}\t{seq}\n")

    print(f"\nSuccess! Fixed file saved to: {output_file}")
    print(f"Header updated to: {num_seqs} {ref_len}")

if __name__ == "__main__":
    main()
