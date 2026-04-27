import sys
import pandas as pd


def read_input():
    data = []
    is_header = True      for line in sys.stdin:
        fields = line.strip().split()
        if is_header:
            is_header = False              continue
        if len(fields) >= 8:
            chrom, start, end, _, frag_type, gene_id, species, dn_value = fields[:8]
            try:
                dn_value = float(dn_value)
            except ValueError:
                dn_value = None
            data.append([chrom, int(start), int(end), frag_type, gene_id, species, dn_value])
    return pd.DataFrame(data, columns=["chrom", "start", "end", "type", "gene_id", "species", "dn_value"])


def calculate_dn_averages(df):
    results = []
    major_sources = {"hm": "hellerii", "ys": "maculatus", "mv": "maculatus"}
    minor_sources = {"hm": "maculatus", "ys": "variatus", "mv": "variatus"}

    for (chrom, start, end, frag_type), group in df.groupby(["chrom", "start", "end", "type"]):
        if frag_type not in major_sources or frag_type not in minor_sources:
            continue  

        major_source = major_sources[frag_type]
        minor_source = minor_sources[frag_type]


        major_dn_values = group[(group["species"] == major_source) & group["dn_value"].notna()]["dn_value"]
        minor_dn_values = group[(group["species"] == minor_source) & group["dn_value"].notna()]["dn_value"]


        valid_genes = group.groupby("gene_id").filter(
            lambda g: major_source in g["species"].values and minor_source in g["species"].values
        )

        if valid_genes.empty:
            continue  

        major_avg = valid_genes[(valid_genes["species"] == major_source)]["dn_value"].mean()
        minor_avg = valid_genes[(valid_genes["species"] == minor_source)]["dn_value"].mean()

        results.append([frag_type, f"{chrom}\t{start}\t{end}", round(major_avg, 5), round(minor_avg, 5)])

    return results


def main():
    df = read_input()
    if df.empty:
        sys.stderr.write("The input file is empty or has an incorrect format.\\n")
        return

    results = calculate_dn_averages(df)

    sys.stdout.write("\t".join(["type", "seg", "avg_dn_major", "avg_dn_minor"]) + "\n")
    for row in results:
        sys.stdout.write("\t".join(map(str, row)) + "\n")

if __name__ == "__main__":
    main()

