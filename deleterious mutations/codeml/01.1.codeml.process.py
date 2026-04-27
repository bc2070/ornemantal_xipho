import sys

def process_sequences(input_file):
    sequence_count = 0
    total_base_count = 0
    output_lines = []
    
    # Define stop codons
    stop_codons = ['TAG', 'TAA', 'TGA']

    def replace_stop_codons(sequence):
        """Replace stop codons with '???' in the sequence."""
        for stop_codon in stop_codons:
            sequence = sequence.replace(stop_codon, '???')
        return sequence

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Increment sequence count
                sequence_count += 1
                # Extract and simplify the sequence ID to only include the species name, and remove ">"
                identifier = line.split(",")[-1].strip()  # Get the part after the last comma
                output_lines.append(identifier)  # Add species name without ">"
            else:
                # Process the sequence
                trimmed_sequence = line[:-3]  # Remove the last 3 bases
                trimmed_sequence = replace_stop_codons(trimmed_sequence)  # Replace stop codons
                total_base_count += len(trimmed_sequence)
                output_lines.append(trimmed_sequence)

    # Calculate average base count
    average_base_count = total_base_count // sequence_count if sequence_count > 0 else 0

    # Prepare the summary line
    summary_line = f"{sequence_count}\t{average_base_count}"

    # Print the output
    print(summary_line)
    for output_line in output_lines:
        print(output_line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    process_sequences(input_file)

