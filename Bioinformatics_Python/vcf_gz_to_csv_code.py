import pandas as pd

# Read VCF file
df = pd.read_csv(
    "ukb-d-30640_irnt.vcf.gz",
    compression="gzip",
    comment="#",
    sep="\t",
    header=None
)

# Add column names manually
df.columns = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
    "ukb-d-30640_irnt"
]

# Save as CSV
df.to_csv(
    "ukb-d-30640_irnt.csv",
    index=False
)

print(df.head())
print("VCF converted successfully")