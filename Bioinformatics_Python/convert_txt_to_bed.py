#!/usr/bin/env python3

input_txt = r"C:\Dinesh_Yadav_2026\Learning_GWAS_January\Learning_Data\LiftOver\try\try.txt"
output_bed = r"C:\Dinesh_Yadav_2026\Learning_GWAS_January\Learning_Data\LiftOver\try\BRC_Dec2024.bed"

with open(input_txt, "r") as fin, open(output_bed, "w") as fout:
    for line in fin:
        line = line.strip()
        if not line:
            continue  # skip empty lines

        # split on any whitespace (space or tab)
        cols = line.split()

        # BED requires at least: chr start end
        if len(cols) < 3:
            continue

        # write TAB-separated BED
        fout.write("\t".join(cols) + "\n")
