
import os
import time
import requests
import pandas as pd

# --------------------------------------------------
# 1. FILE PATHS
# --------------------------------------------------

folder = (
    r"D:\COMMON_FOLDER_FOR_TEJAS_DINESH"
    r"\BMC_Lead_SNP_buccal_mucosa_cancer_in_India_and_multi_ancestry_meta_analysis_ebiomedicine"
)

input_file = os.path.join(folder, "Without Trate names.xlsx")
output_file = os.path.join(folder, "Ensembl_Gene_Annotations.xlsx")

# --------------------------------------------------
# 2. READ YOUR EXCEL FILE
# --------------------------------------------------

df = pd.read_excel(input_file)

if "trait" not in df.columns:
    raise ValueError(f"Column 'trait' not found. Found: {df.columns.tolist()}")

gene_ids = (
    df["trait"]
    .dropna()
    .astype(str)
    .str.strip()
    .str.replace(r"\.\d+$", "", regex=True)
    .drop_duplicates()
    .tolist()
)

print(f"Total unique Ensembl IDs: {len(gene_ids)}")

# --------------------------------------------------
# 3. QUERY THE OFFICIAL ENSEMBL REST API
# --------------------------------------------------

url = "https://rest.ensembl.org/lookup/id"

headers = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

results = []
batch_size = 100

with requests.Session() as session:

    for start in range(0, len(gene_ids), batch_size):
        batch = gene_ids[start:start + batch_size]

        for attempt in range(3):
            try:
                response = session.post(
                    url,
                    headers=headers,
                    json={"ids": batch},
                    params={"species": "homo_sapiens"},
                    timeout=60
                )

                if response.status_code == 429 or response.status_code >= 500:
                    print("Ensembl temporarily unavailable. Retrying...")
                    time.sleep(3 * (attempt + 1))
                    continue

                response.raise_for_status()
                data = response.json()
                break

            except (requests.RequestException, ValueError) as error:
                if attempt == 2:
                    print("Batch request failed:", error)
                    data = {}
                else:
                    time.sleep(3 * (attempt + 1))

        for gene_id in batch:
            annotation = data.get(gene_id, {})

            if not isinstance(annotation, dict):
                annotation = {}

            symbol = annotation.get("display_name", "")
            description = annotation.get("description", "")

            results.append({
                "Ensembl_Gene_ID": gene_id,
                "Gene_Symbol": symbol,
                "Gene_Name": description,
                "Biotype": annotation.get("biotype", ""),
                "Chromosome": annotation.get("seq_region_name", ""),
                "Status": (
                    "Found" if symbol or description
                    else "Not found in Ensembl response"
                )
            })

# --------------------------------------------------
# 4. SAVE THE RESULTS
# --------------------------------------------------

result_df = pd.DataFrame(results)
result_df.to_excel(output_file, index=False)

print("\nCompleted!")
print("Output file:", output_file)
print(result_df.to_string(index=False))

print("\nLookup status:")
print(result_df["Status"].value_counts())
