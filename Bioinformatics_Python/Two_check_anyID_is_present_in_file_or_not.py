import pandas as pd

file_path = r"C:\Dinesh_Yadav_2026\GBC_PheWAS\Results\New_extract_columns\Check_For_Clumping_Data\G2_N_1875_withHLA_META-with_G1_O_48_withHLA_BothBeagle.meta_with-BETA.csv_With_BETA-SE-OR_n_ConfIntervals_95_with-AllelFreqs-CombinedData_WithrsID.csv"   # change this
target_snp = "rs1059542"

found = False

for chunk in pd.read_csv(file_path, chunksize=300000, dtype=str):
    # Clean SNP column
    chunk["SNP"] = chunk["SNP"].str.strip()
    
    # Check if SNP exists
    if target_snp in chunk["SNP"].values:
        print(f"✅ SNP {target_snp} FOUND in file")
        found = True
        break

if not found:
    print(f"❌ SNP {target_snp} NOT FOUND in file")