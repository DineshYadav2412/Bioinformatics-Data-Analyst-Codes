setwd("D:/COMMON_FOLDER_FOR_TEJAS_DINESH/Traits_selected_by_Gib_30_sep_2026/1_finn_b_K11_GALLBILPANC")

library(data.table)

# --------------------------------------------------
# 1. Read VCF
# --------------------------------------------------

read_data <- fread(
  "finn-b-K11_GALLBILPANC.vcf.gz",
  skip = "#CHROM",
  header = TRUE,
  sep = "\t"
)

print(dim(read_data))

# --------------------------------------------------
# 2. Rename chromosome column
# --------------------------------------------------

setnames(read_data, "#CHROM", "CHROM")

# --------------------------------------------------
# 3. Identify the last/sample column
# --------------------------------------------------

sample_col <- names(read_data)[10]

# --------------------------------------------------
# 4. Split the sample column
# --------------------------------------------------

split_values <- tstrsplit(
  read_data[[sample_col]],
  ":",
  fill = NA,
  fixed = TRUE
)

# --------------------------------------------------
# 5. Take only first 4 values
# ES = BETA
# SE = SE
# LP = LOGP
# AF = Allele_Frequency
# --------------------------------------------------

split_df <- data.table(
  BETA = split_values[[1]],
  SE = split_values[[2]],
  LOGP = split_values[[3]],
  Allele_Frequency = split_values[[4]]
)

# --------------------------------------------------
# 6. Add extracted columns
# --------------------------------------------------

read_data[, BETA := split_df$BETA]
read_data[, SE := split_df$SE]
read_data[, LOGP := split_df$LOGP]
read_data[, Allele_Frequency := split_df$Allele_Frequency]

# --------------------------------------------------
# 7. Remove unnecessary FORMAT/sample columns
# --------------------------------------------------

read_data[, FORMAT := NULL]
read_data[, (sample_col) := NULL]

# --------------------------------------------------
# 8. Convert numeric columns
# --------------------------------------------------

read_data[, BETA := as.numeric(BETA)]
read_data[, SE := as.numeric(SE)]
read_data[, LOGP := as.numeric(LOGP)]
read_data[, Allele_Frequency := as.numeric(Allele_Frequency)]

# --------------------------------------------------
# 9. Save CSV
# --------------------------------------------------

fwrite(
  read_data,
  "finn_b_K11_GALLBILPANC.csv"
)

# --------------------------------------------------
# 10. Check final output
# --------------------------------------------------

print(dim(read_data))
print(names(read_data))
print(head(read_data))

