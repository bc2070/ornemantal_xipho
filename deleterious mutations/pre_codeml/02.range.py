import os
import sys
import shutil

def classify_file(file_path):

    with open(file_path, "r") as f:
        lines = f.readlines()

    new_lines = [line.upper() if not line.startswith(">") else line for line in lines]
    with open(file_path, "w") as f:
        f.writelines(new_lines)

    gambusia = None
    poecilia = None

    for i, line in enumerate(new_lines):
        if line.startswith(">") and "GAMBUSIA_AFFINIS" in line.upper():
            gambusia = new_lines[i + 1].strip()  
        elif line.startswith(">") and "POECILIA_FORMOSA" in line.upper():
            poecilia = new_lines[i + 1].strip()


    if gambusia is None and poecilia is None:
        target_dir = "none"
    else:

        gambusia_all_n = gambusia is not None and all(base == "N" for base in gambusia)
        poecilia_all_n = poecilia is not None and all(base == "N" for base in poecilia)


        if gambusia_all_n and not poecilia_all_n:
            target_dir = "single_o2"
        elif poecilia_all_n and not gambusia_all_n:
            target_dir = "single_o1"
        elif not gambusia_all_n and not poecilia_all_n:
            target_dir = "o1o2"
        elif gambusia_all_n and poecilia_all_n:
            target_dir = "none"
        else:
            print(f"Error: unknown {file_path}")
            return

    base_name = os.path.basename(file_path).replace(".cds.fa", "")
    sub_dir = os.path.join(target_dir, base_name)

    os.makedirs(sub_dir, exist_ok=True)
    shutil.move(file_path, os.path.join(sub_dir, os.path.basename(file_path)))
    print(f"file {file_path} classified {sub_dir}/")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python classify_sequences.py <file>")
        sys.exit(1)

    file_path = sys.argv[1]
    classify_file(file_path)

