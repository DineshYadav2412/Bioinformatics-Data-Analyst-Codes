import pandas as pd
import requests
import time
import re
import os
from tqdm import tqdm

# ==========================================
# 1. FILE PATHS
# ==========================================

input_file = "rsid_extract_data.xlsx"

output_file = "rsid_extract_data_with_dbSNP_GRCh37.xlsx"

# Ensembl GRCh37 REST API
BASE_URL = "https://grch37.rest.ensembl.org"

# ==========================================
# 2. READ EXCEL FILE
# ==========================================

df = pd.read_excel(input_file)

print("Total rows:", len(df))
print("Excel columns:", df.columns.tolist())

if "CHR_BP" not in df.columns:
    raise ValueError("CHR_BP column not found!")

# ==========================================
# 3. EXTRACT CHROMOSOME AND BP
# ==========================================

# Supports examples:
# 6:32623550
# chr6:32623550
# 6_32623550

parsed = df["CHR_BP"].astype(str).str.extract(
    r"^\s*(?:chr)?(\d+|X|Y|MT|M)[:_](\d+)\s*$",
    flags=re.IGNORECASE
)

df["CHR"] = parsed[0].str.upper()
df["BP"] = pd.to_numeric(parsed[1], errors="coerce")

df["CHR"] = df["CHR"].replace({"M": "MT"})

# ==========================================
# 4. PREPARE API SESSION
# ==========================================

session = requests.Session()

session.headers.update({
    "Content-Type": "application/json",
    "Accept": "application/json"
})

# ==========================================
# 5. FUNCTION TO EXTRACT dbSNP rsID
# ==========================================

def get_rsid(chrom, bp):

    if pd.isna(chrom) or pd.isna(bp):
        return None

    chrom = str(chrom).replace("chr", "")
    pos = int(bp)

    url = (
        f"{BASE_URL}/overlap/region/human/"
        f"{chrom}:{pos}-{pos}"
    )

    params = {"feature": "variation"}

    for attempt in range(5):

        try:
            response = session.get(
                url,
                params=params,
                timeout=30
            )

            if response.status_code == 429:
                time.sleep(2 * (attempt + 1))
                continue

            response.raise_for_status()

            data = response.json()

            # Keep records at exact BP
            ids = sorted({
                item["id"]
                for item in data
                if item.get("start") == pos
                and item.get("id", "").startswith("rs")
            })

            return ";".join(ids) if ids else None

        except requests.exceptions.RequestException:
            time.sleep(2 * (attempt + 1))

    return None

# ==========================================
# 6. LOOK UP ALL COORDINATES
# ==========================================

# Cache repeated coordinates to avoid
# querying the same position multiple times.

unique_coords = (
    df[["CHR", "BP"]]
    .dropna()
    .drop_duplicates()
)

rsid_lookup = {}

for row in tqdm(
    unique_coords.itertuples(index=False),
    total=len(unique_coords),
    desc="Looking up GRCh37 variants"
):

    chrom = row.CHR
    bp = int(row.BP)

    key = (chrom, bp)

    rsid_lookup[key] = get_rsid(chrom, bp)

    # Be polite to the public API
    time.sleep(0.12)

    # Save lookup progress periodically
    if len(rsid_lookup) % 250 == 0:
        pd.DataFrame(
            [
                {"CHR": k[0], "BP": k[1], "dbSNP_rsid": v}
                for k, v in rsid_lookup.items()
            ]
        ).to_csv("dbSNP_lookup_checkpoint.csv", index=False)

# ==========================================
# 7. ADD rsIDs TO EXCEL DATA
# ==========================================

df["dbSNP_rsid"] = [
    rsid_lookup.get((chrom, int(bp)))
    if pd.notna(chrom) and pd.notna(bp)
    else None
    for chrom, bp in zip(df["CHR"], df["BP"])
]

# ==========================================
# 8. COMPARE WITH EXISTING variant_id
# ==========================================

def compare_ids(existing, extracted):

    if pd.isna(extracted) or not str(extracted).strip():
        return "No dbSNP match"

    if pd.isna(existing) or not str(existing).strip():
        return "Existing ID missing"

    old_ids = {
        x.strip().lower()
        for x in str(existing).split(";")
    }

    new_ids = {
        x.strip().lower()
        for x in str(extracted).split(";")
    }

    if old_ids & new_ids:
        return "Matched"

    return "Different / Check"


if "variant_id" in df.columns:

    df["Match_Status"] = [
        compare_ids(old, new)
        for old, new in zip(
            df["variant_id"],
            df["dbSNP_rsid"]
        )
    ]

# ==========================================
# 9. EXPORT FINAL EXCEL
# ==========================================

df.to_excel(output_file, index=False)

print("\nExtraction completed!")
print("Output saved:", os.path.abspath(output_file))

if "Match_Status" in df.columns:
    print("\nMatch summary:")
    print(df["Match_Status"].value_counts(dropna=False))