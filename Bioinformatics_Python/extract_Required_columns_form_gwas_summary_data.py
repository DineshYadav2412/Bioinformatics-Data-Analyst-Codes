import pandas as pd

# -------------------------------------------------------------
# Script Name : prepare_gwas_for_pics2.py
#
# Purpose:
# Reads the original GWAS summary statistics file, extracts only
# the columns required for PICS2 fine-mapping analysis, renames
# them to standard names, removes incomplete and duplicate SNPs,
# and saves a cleaned CSV file.
# -------------------------------------------------------------

# Read the original GWAS summary statistics file
df = pd.read_csv(
    "G2_N_1875_withHLA_META-with_G1_O_48_withHLA_BothBeagle.meta_with-BETA.csv_With_BETA-SE-OR_n_ConfIntervals_95_with-AllelFreqs-CombinedData_WithrsID.csv"
)

# Select only the columns required for PICS2 analysis
gwas = df[
    [
        "CHR",                  # Chromosome
        "BP",                   # Base pair position
        "rsID",                 # SNP identifier
        "A1_MinorAllele",       # Effect (minor) allele
        "A2_MajorAllele",       # Other (major) allele
        "P",                    # P-value
        "BETA",                 # Effect size
        "SE",                   # Standard error
        "A1_MinorAlleleFreq",   # Effect allele frequency
        "NCHROBS",              # Sample size
    ]
].copy()

# Rename columns to the standard names expected by PICS2
gwas.columns = [
    "CHR",
    "BP",
    "SNP",
    "A1",
    "A2",
    "P",
    "BETA",
    "SE",
    "EAF",
    "N",
]

# Remove rows with missing values in essential columns
gwas = gwas.dropna(
    subset=["CHR", "BP", "SNP", "P", "BETA", "SE"]
)

# Remove duplicate SNPs, keeping the first occurrence
gwas = gwas.drop_duplicates(subset="SNP")

# Save the cleaned GWAS dataset for downstream PICS2 analysis
gwas.to_csv("GWAS_PICS2_Input.csv", index=False)

# Display basic summary
print("Rows:", len(gwas))
print(gwas.head())