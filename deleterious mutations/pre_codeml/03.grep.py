import os

def extract_sequences(file_path, output_dir, patterns, output_names):

    with open(file_path, "r") as f:
        lines = f.readlines()

    for pattern, output_name in zip(patterns, output_names):
        selected_sequences = []
        for i, line in enumerate(lines):
            if line.startswith(">") and any(species in line.upper() for species in pattern):
                selected_sequences.append(line)  
                selected_sequences.append(lines[i + 1].upper()) 

        if selected_sequences:
            output_path = os.path.join(output_dir, output_name)
            with open(output_path, "w") as out_file:
                out_file.writelines(selected_sequences)
            print(f"file: {output_path}")


def process_directory(base_dir, sub_dir, patterns, output_names):
    target_dir = os.path.join(base_dir, sub_dir)
    if not os.path.exists(target_dir):
        print(f"dir {target_dir} no exist，skip")
        return

    for folder in os.listdir(target_dir):
        folder_path = os.path.join(target_dir, folder)
        if os.path.isdir(folder_path):
            for file_name in os.listdir(folder_path):
                if file_name.endswith(".cds.fa"):
                    file_path = os.path.join(folder_path, file_name)
                    extract_sequences(
                        file_path=file_path,
                        output_dir=folder_path,
                        patterns=patterns,
                        output_names=output_names
                    )


if __name__ == "__main__":
    base_directory = "/data2/projects/jgeng/peak_hm_GO/genespace/Xma/outgroup_prot/cds_outgroups/new_allregion"  
    # o1o2 
    o1o2_patterns = [
        ["GAMBUSIA_AFFINIS", "POECILIA_FORMOSA", "XIPHOPHORUS_HELLERII"],
        ["GAMBUSIA_AFFINIS", "POECILIA_FORMOSA", "XIPHOPHORUS_MACULATUS"],
        ["GAMBUSIA_AFFINIS", "POECILIA_FORMOSA", "XIPHOPHORUS_VARIATUS"]
    ]
    o1o2_output_names = ["o1o2_hellerii.cds",
                         "o1o2_maculatus.cds",
                         "o1o2_variatus.cds"]
    process_directory(base_directory, "o1o2", o1o2_patterns, o1o2_output_names)

    # single_o1 
    single_o1_patterns = [
        ["GAMBUSIA_AFFINIS", "XIPHOPHORUS_HELLERII"],
        ["GAMBUSIA_AFFINIS", "XIPHOPHORUS_MACULATUS"],
        ["GAMBUSIA_AFFINIS", "XIPHOPHORUS_VARIATUS"]
    ]
    single_o1_output_names = ["gambusia_hellerii.cds",
                              "gambusia_maculatus.cds",
                              "gambusia_variatus.cds"]
    process_directory(base_directory, "single_o1", single_o1_patterns, single_o1_output_names)

    # single_o2 
    single_o2_patterns = [
        ["POECILIA_FORMOSA", "XIPHOPHORUS_HELLERII"],
        ["POECILIA_FORMOSA", "XIPHOPHORUS_MACULATUS"],
        ["POECILIA_FORMOSA", "XIPHOPHORUS_VARIATUS"]
    ]
    single_o2_output_names = ["poecilia_hellerii.cds",
                              "poecilia_maculatus.cds",
                              "poecilia_variatus.cds"]
    process_directory(base_directory, "single_o2", single_o2_patterns, single_o2_output_names)

