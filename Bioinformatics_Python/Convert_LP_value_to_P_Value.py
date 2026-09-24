import pandas as pd
import numpy as np

# =========================================================
# INPUT / OUTPUT FILES
# =========================================================

INPUT_FILE = "ebi-a-GCST005523_final_output.csv"
OUTPUT_FILE = "ebi-a-GCST005523_final_output_final_output_with_pvalue.csv"

try:
    # =====================================================
    # READ FILE
    # =====================================================

    df = pd.read_csv(INPUT_FILE, dtype=str, low_memory=False)

    print("\nFile Loaded Successfully")

    # =====================================================
    # CONVERT LP COLUMN TO NUMERIC
    # =====================================================

    lp_numeric = pd.to_numeric(df["LP"], errors="coerce")

    # =====================================================
    # CALCULATE P-VALUE
    # Formula:
    # P = 10^(-LP)
    # =====================================================

    p_values = 10 ** (-lp_numeric)

    # =====================================================
    # ADD NEW COLUMN AFTER LP
    # =====================================================

    df.insert(
        df.columns.get_loc("LP") + 1,
        "P_VALUE",
        p_values
    )

    # =====================================================
    # OPTIONAL:
    # FORMAT P_VALUE IN SCIENTIFIC NOTATION
    # =====================================================

    df["P_VALUE"] = df["P_VALUE"].apply(lambda x: f"{x:.10e}")

    # =====================================================
    # SAVE OUTPUT
    # =====================================================

    df.to_csv(OUTPUT_FILE, index=False)

    print(f"\nOutput Saved Successfully")
    print(f"Output File: {OUTPUT_FILE}")

    # =====================================================
    # PREVIEW
    # =====================================================

    print("\nPreview:")
    print(df.head())

except Exception as e:
    print(f"\nERROR: {e}")