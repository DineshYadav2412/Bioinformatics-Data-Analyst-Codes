import pandas as pd

input_file = "G2_N_1875_withHLA_META-with_G1_O_48_withHLA_BothBeagle.meta_with-BETA.csv_With_BETA-SE-OR_n_ConfIntervals_95_with-AllelFreqs-CombinedData_WithrsID.csv"

# Read only N column
df = pd.read_csv(input_file, usecols=["N"])

# Convert to numeric
df["N"] = pd.to_numeric(df["N"], errors="coerce")

# Get values NOT equal to 2
other_values = df[df["N"] != 2]["N"].dropna().unique()

if len(other_values) == 0:
    print("✅ All values in N column are 2")
else:
    print("❌ Found values other than 2:")
    print(other_values)