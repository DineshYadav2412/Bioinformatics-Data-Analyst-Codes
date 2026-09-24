import pandas as pd

# =========================================================
# VCF/CSV Processing Script
# Author: Dinesh Yadav
# Purpose:
# Split ES:SE:LP:ID values into separate columns
# and save clean output into CSV
# =========================================================

# ---------- INPUT / OUTPUT FILES ----------
INPUT_FILE = "ebi-a-GCST005523.csv"
OUTPUT_FILE = "final_output.csv"

try:
    # ---------- READ FILE ----------
    # dtype=str preserves exact values without rounding
    df = pd.read_csv(INPUT_FILE, dtype=str, low_memory=False)

    print(f"\nFile Loaded Successfully")
    print(f"Total Rows   : {len(df)}")
    print(f"Total Columns: {len(df.columns)}")

    # ---------- GET LAST COLUMN ----------
    last_col = df.columns[-1]

    print(f"\nProcessing Column: {last_col}")

    # ---------- CLEAN COLUMN ----------
    # Remove unwanted text if present
    df[last_col] = (
        df[last_col]
        .str.replace("ES:SE:LP:ID", "", regex=False)
        .str.strip()
    )

    # ---------- SPLIT VALUES ----------
    split_cols = df[last_col].str.split(":", expand=True)

    # ---------- VALIDATION ----------
    if split_cols.shape[1] < 4:
        raise ValueError(
            "Expected 4 values (ES:SE:LP:ID) but found fewer."
        )

    # ---------- CREATE NEW COLUMNS ----------
    df["ES"] = split_cols[0]
    df["SE"] = split_cols[1]
    df["LP"] = split_cols[2]
    df["Variant_ID"] = split_cols[3]

    # ---------- FINAL REQUIRED COLUMNS ----------
    required_columns = [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "INFO",
        "ES",
        "SE",
        "LP",
        "Variant_ID"
    ]

    # Check missing columns
    missing_cols = [col for col in required_columns if col not in df.columns]

    if missing_cols:
        raise ValueError(f"Missing Columns: {missing_cols}")

    final_df = df[required_columns]

    # ---------- SAVE OUTPUT ----------
    final_df.to_csv(OUTPUT_FILE, index=False)

    print(f"\nOutput Saved Successfully")
    print(f"Output File: {OUTPUT_FILE}")

    # ---------- PREVIEW ----------
    print("\nPreview:")
    print(final_df.head())

except FileNotFoundError:
    print(f"\nERROR: File not found -> {INPUT_FILE}")

except Exception as e:
    print(f"\nERROR: {e}")