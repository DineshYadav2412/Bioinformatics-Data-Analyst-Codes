import pandas as pd

# given a file containing the CHR and BP data of SNPs, this code will use the cytoBand information file from UCSC genome browser and find out which cytoband the SNPs lies within & return that information

cyto = str(input("Enter the full path to the cytoband information file, containing the CHR and From-BP, To-BP & CytoBand columns\t:\t"))
toFind = str(input("Enter the full path to the .CSV file containing the CHR and BP values (those exact column names) for the SNPs for which you want the cytoband\t:\t"))

output_path = str(input("Enter the full path including the output file name with the .csv extension for the output file\t:\t"))

cyto_df = pd.read_csv(cyto)
find_df = pd.read_csv(toFind)


find_df["CytoBand"] = ""
for i, row in find_df.iterrows():
    chr = row["CHR"]
    bp = row["BP"]
    row_found = cyto_df[
        (cyto_df["CHR"] == chr) & 
        (cyto_df["From-BP"] <= bp) & 
        (cyto_df["To-BP"] >= bp)
        ]
    if len(row_found) > 0:
        cy = row_found["CytoBand"].iloc[0]
    else:
        cy = None
    find_df.at[i, "CytoBand"] = cy
    
    
find_df.to_csv(output_path, index=False)

