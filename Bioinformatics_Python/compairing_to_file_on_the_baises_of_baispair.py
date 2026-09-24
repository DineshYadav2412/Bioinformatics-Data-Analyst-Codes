import os

file1 = r"C:\Dinesh_Yadav_2026\New_Data_19_May_2026\Beadchip_illumina_project\Human_omini_express\1_Human_omini_express_Full_QC\Beagle_Imputation\Human_omini_chrX.bim"
file2 = r"C:\Dinesh_Yadav_2026\New_Data_19_May_2026\Beadchip_illumina_project\All_Beagle_Imputation_reference_files_Are_Here\Beagle_Reference_file\1000G_ChrX\New folder\New1000Genome\chrX_ref.bim"

output_dir = r"C:\Dinesh_Yadav_2026\New_Data_19_May_2026\Beadchip_illumina_project\All_Beagle_Imputation_reference_files_Are_Here\Beagle_Reference_file\1000G_ChrX\comapred_data\comapaired_with_Human_omini__1000Genome_ChrX"
os.makedirs(output_dir, exist_ok=True)

# Read column 2 from file 1
set1 = set()
with open(file1) as f:
    for line in f:
        cols = line.strip().split()
        if len(cols) >= 2:
            set1.add(cols[1])

# Read column 2 from file 2
set2 = set()
with open(file2) as f:
    for line in f:
        cols = line.strip().split()
        if len(cols) >= 2:
            set2.add(cols[1])

# Common
common = set1 & set2

# Not common (present in only one file)
uncommon = set1 ^ set2

# Save common
with open(os.path.join(output_dir, "common.txt"), "w") as f:
    for x in sorted(common):
        f.write(x + "\n")

# Save uncommon
with open(os.path.join(output_dir, "uncommon.txt"), "w") as f:
    for x in sorted(uncommon):
        f.write(x + "\n")

print("Common:", len(common))
print("Uncommon:", len(uncommon))








# compare_bp_lists.py

# gsa_file = r"C:\Dinesh_Yadav_2026\New_Data_19_May_2026\Beadchip_illumina_project\All_Beagle_Imputation_reference_files_Are_Here\Beagle_Reference_file\1000G_ChrX\New folder\1000Genomens_refe_bp_list.txt"
# confluence_file = r"C:\Dinesh_Yadav_2026\New_Data_19_May_2026\Beadchip_illumina_project\All_Beagle_Imputation_reference_files_Are_Here\Beagle_Reference_file\1000G_ChrX\New folder\GSA_Confluence_bp_list.txt"

# common_output = r"C:\Dinesh_Yadav_2026\New_Data_19_May_2026\Beadchip_illumina_project\All_Beagle_Imputation_reference_files_Are_Here\Beagle_Reference_file\1000G_ChrX\New folder\common_bp.txt"
# not_common_output = r"C:\Dinesh_Yadav_2026\New_Data_19_May_2026\Beadchip_illumina_project\All_Beagle_Imputation_reference_files_Are_Here\Beagle_Reference_file\1000G_ChrX\New folder\not_common_bp.txt"

# # Read files
# with open(gsa_file, "r") as f:
#     gsa = set(line.strip() for line in f if line.strip())

# with open(confluence_file, "r") as f:
#     confluence = set(line.strip() for line in f if line.strip())

# # Common positions
# common = sorted(gsa & confluence)

# # Not common in either file (symmetric difference)
# not_common = sorted(gsa ^ confluence)

# # Save outputs
# with open(common_output, "w") as f:
#     f.write("\n".join(common))

# with open(not_common_output, "w") as f:
#     f.write("\n".join(not_common))

# print(f"Common BPs: {len(common)}")
# print(f"Not Common BPs: {len(not_common)}")
# print(f"Saved: {common_output}")
# print(f"Saved: {not_common_output}")