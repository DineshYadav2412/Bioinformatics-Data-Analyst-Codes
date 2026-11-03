setwd("E:/LC_phewas_27_oct_2026/LC_data_27_oct_2026/EAF_PROCESS")


library(data.table)
library(dplyr)
# Step 1: Read both datasets
file1 <- fread("3_indai_DuplisCleaned_Logistic_All_10PCs_Age_Gen_ChewYN_OR.assoc.logistic.csv_With_BETA-SE-OR_n_ConfIntervals_95_output_with_CHRBP.csv")
file2 <- fread("ebi-a-GCST004747_COLUMNS_FIXED_FOR_META.csv_With_BETA-SE-OR_n_ConfIntervals_95.csv")
dim(file1)
dim(file2)
# Step 2: Convert CHR and BP to same data type
file1[, CHR := as.character(CHR)]
file2[, CHR := as.character(CHR)]

file1[, BP := as.numeric(BP)]
file2[, BP := as.numeric(BP)]

# Step 3: Find common CHR and BP and merge datasets
common_data <- merge(
  file1,
  file2,
  by = c("CHR", "BP"),
  all = FALSE,
  suffixes = c("", "_mckay")
)

# Step 4: Check dimensions
cat("File 1 dimensions:", dim(file1), "\n")
cat("File 2 dimensions:", dim(file2), "\n")
cat("Common merged data dimensions:", dim(common_data), "\n")

# Step 5: Save common merged dataset
fwrite(common_data, "Common_CHR_BP_Merged.csv")

# Step 6: Count common SNPs
cat("Total common CHR-BP combinations:", 
    uniqueN(common_data, by = c("CHR", "BP")), "\n")

head(common_data)
names(common_data)

######################################################################################################################################################################
# For match allels
######################################################################################################################################################################


# Step 1: Read common SNP dataset
data <- fread("Common_CHR_BP_Merged.csv")

# Step 2: Compare alleles row by row
data[, Allele_Status := fifelse(
  A1_MINORALLELE == A1_MINORALLELE_mckay &
    A2_MAJORALLELE == A2_MAJORALLELE_mckay,
  "Match",
  "Unmatch"
)]

# Step 3: Count Match and Unmatch
match_count <- sum(data$Allele_Status == "Match", na.rm = TRUE)
unmatch_count <- sum(data$Allele_Status == "Unmatch", na.rm = TRUE)

# Step 4: Print results
cat("Total SNPs:", nrow(data), "\n")
cat("Match:", match_count, "\n")
cat("Unmatch:", unmatch_count, "\n")

# Step 5: Percentage
cat("Match percentage:", 
    round(match_count / nrow(data) * 100, 2), "%\n")

cat("Unmatch percentage:", 
    round(unmatch_count / nrow(data) * 100, 2), "%\n")

# Step 6: Save dataset with allele status
fwrite(data, "Common_SNPs_Allele_Comparison.csv")


ALLELE_DATA <- fread("Common_SNPs_Allele_Comparison.csv")
dim(ALLELE_DATA)
head(ALLELE_DATA)


result <- ALLELE_DATA %>%
  filter(CHR == 1 & BP == 768253)

View(result)

######################################################################################################################################################################
# For match allels
######################################################################################################################################################################


library(data.table)

# Step 1: Read dataset
data <- fread("Common_SNPs_Allele_Comparison.csv")

# Step 2: Identify Unmatch rows
unmatch_rows <- which(data$Allele_Status == "Unmatch")

cat("Total Unmatch before flipping:", length(unmatch_rows), "\n")

# Step 3: Flip alleles for Unmatch rows
temp <- data$A1_MINORALLELE_mckay[unmatch_rows]

data$A1_MINORALLELE_mckay[unmatch_rows] <-
  data$A2_MAJORALLELE_mckay[unmatch_rows]

data$A2_MAJORALLELE_mckay[unmatch_rows] <- temp

# Step 4: Change BETA sign
data$BETA_mckay[unmatch_rows] <-
  data$BETA_mckay[unmatch_rows] * -1

# Step 5: Calculate reciprocal of OR
data$OR_mckay[unmatch_rows] <-
  1 / data$OR_mckay[unmatch_rows]

# Step 6: Flip allele frequency
data$A1_MINORALLELEFREQ_mckay[unmatch_rows] <-
  1 - data$A1_MINORALLELEFREQ_mckay[unmatch_rows]

# Step 7: Recheck allele status
data[, Allele_Status := fifelse(
  A1_MINORALLELE == A1_MINORALLELE_mckay &
  A2_MAJORALLELE == A2_MAJORALLELE_mckay,
  "Match",
  "Unmatch"
)]

# Step 8: Count final Match and Unmatch
cat("\nFinal allele comparison:\n")
print(data[, .N, by = Allele_Status])

# Step 9: Save updated dataset
fwrite(data, "Common_SNPs_Alleles_Flipped.csv")

cat("\nUpdated dataset saved successfully!\n")



######################################################################################################################################################################
# For flip
######################################################################################################################################################################


# Step 1: Read dataset
data <- fread("Common_SNPs_Allele_Comparison.csv")

# Step 2: Identify Unmatch rows
unmatch_rows <- which(data$Allele_Status == "Unmatch")

cat("Total Unmatch before flipping:", length(unmatch_rows), "\n")

# Step 3: Flip alleles for Unmatch rows
temp <- data$A1_MINORALLELE_mckay[unmatch_rows]

data$A1_MINORALLELE_mckay[unmatch_rows] <-
  data$A2_MAJORALLELE_mckay[unmatch_rows]

data$A2_MAJORALLELE_mckay[unmatch_rows] <- temp

# Step 4: Change BETA sign
data$BETA_mckay[unmatch_rows] <-
  data$BETA_mckay[unmatch_rows] * -1

# Step 5: Calculate reciprocal of OR
data$OR_mckay[unmatch_rows] <-
  1 / data$OR_mckay[unmatch_rows]

# Step 6: Flip allele frequency
data$A1_MINORALLELEFREQ_mckay[unmatch_rows] <-
  1 - data$A1_MINORALLELEFREQ_mckay[unmatch_rows]

# Step 7: Recheck allele status
data[, Allele_Status := fifelse(
  A1_MINORALLELE == A1_MINORALLELE_mckay &
    A2_MAJORALLELE == A2_MAJORALLELE_mckay,
  "Match",
  "Unmatch"
)]

# Step 8: Count final Match and Unmatch
cat("\nFinal allele comparison:\n")
print(data[, .N, by = Allele_Status])

# Step 9: Save updated dataset
fwrite(data, "Common_SNPs_Alleles_Flipped.csv")

cat("\nUpdated dataset saved successfully!\n")


######################################################################################################################################################################
# Delete Unmatch data
######################################################################################################################################################################
data <- fread("Common_SNPs_Alleles_Flipped.csv")

unmatch_count <- sum(data$Allele_Status == "Unmatch", na.rm = TRUE)
cat("Unmatch Count:", unmatch_count, "\n")

data <- data[data$Allele_Status != "Unmatch", ]
table(data$Allele_Status, useNA = "ifany")
dim(data)

fwrite(data, "Common_SNPs_Alleles_Flipped_removed_4194.csv")



######################################################################################################################################################################
# Calculating Indian GWAS results Effect Allele Sample Count (EASC-Ind)
###########################################################################################################################################


# Step 1: Read updated dataset
data <- fread("Common_SNPs_Alleles_Flipped_removed_4194.csv")

# Step 2: Create new column EASC_Ind
data[, EASC_McKay := 9859 * A1_MINORALLELEFREQ_mckay]

# Step 3: Check first 10 rows
head(data[, .(A1_MINORALLELEFREQ_mckay, EASC_McKay)], 10)
head(data)
# Step 4: Save updated dataset
fwrite(data, "Common_SNPs_With_EASC_Ind_&_EASC_McKay.csv")

cat("EASC_Ind column created successfully!\n")




######################################################################################################################################################################
# Meta-Effect Allele Samples Count (MEASC) = EASC-Ind + EASC-McKay
###########################################################################################################################################

# Step 1: Read updated dataset
data <- fread("Common_SNPs_With_EASC_Ind_&_EASC_McKay.csv")

# Step 2: Create new column EASC_Ind
data[, Meta_Effect_Allele_Samples_Count := EASC_Ind + EASC_McKay]

# Step 3: Check first 10 rows
head(data[, .(EASC_Ind, EASC_McKay, Meta_Effect_Allele_Samples_Count)], 10)
head(data)
# Step 4: Save updated dataset
fwrite(data, "Common_SNPs_With_Meta_Effect_Allele_Samples_Count_MEASC.csv")


######################################################################################################################################################################
# Meta-Effect Allele Samples Count (MEASC) = EASC-Ind + EASC-McKay
###########################################################################################################################################

# Step 1: Read updated dataset
data <- fread("Common_SNPs_With_Meta_Effect_Allele_Samples_Count_MEASC.csv")

# Step 2: Create new column EASC_Ind
data[, Meta_EAF := Meta_Effect_Allele_Samples_Count / 15197]

# Step 3: Check first 10 rows
head(data[, .(EASC_Ind, EASC_McKay, Meta_Effect_Allele_Samples_Count, Meta_EAF)], 10)
head(data)
# Step 4: Save updated dataset
fwrite(data, "Common_SNPs_With_Meta_Effect_Allele_Frequency.csv")


Final_data <- fread("Common_SNPs_With_Meta_Effect_Allele_Frequency.csv")
dim(Final_data)
head(Final_data)




######################################################################################################################################################################
# Meta_Effect_Allele_Frequency column from File1 into File2 by matching CHR and BP.
###########################################################################################################################################



library(data.table)

setwd("D:/COMMON_FOLDER_FOR_TEJAS_DINESH/COMMON_FILES_THAT_USE_MOST_OF_THE_TIME/LC_DATA_ALLELE_FREQUENCY")

# Read both files
File1 <- fread("ebi-a-GCST004747_With_rsid.csv")
head(File1)

#library(data.table)

# Read file
#File1 <- fread("ebi-a-GCST004747_COLUMNS_FIXED_FOR_META.csv_With_BETA-SE-OR_n_ConfIntervals_95.csv")

# Rename SNP to rsid
#setnames(File1, "SNP", "rsid")

# Save updated file
#fwrite(File1, "ebi-a-GCST004747_With_rsid.csv")

# Check
#head(File1)


File2 <- fread("LC-Meta_McKay_AllSamples_AgeGen10PCsChewYN_CommonSNPsOnly_csv_With_BETA-SE-OR_n_ConfIntervals_95_with_Effect_allele_frequecy.csv")
dim(File2)
head(File2)
# Select only the required columns from File1
File1_subset <- File1[, .(CHR, BP, rsid)]

# Match CHR and BP and add Meta_Effect_Allele_Frequency to File2
File2_updated <- merge(
  File2,
  File1_subset,
  by = c("CHR", "BP"),
  all.x = TRUE
)

# Check results
head(File2_updated)
dim(File2_updated)

# Count matched and unmatched rows
sum(!is.na(File2_updated$rsid))
sum(is.na(File2_updated$rsid))

# Save updated File2
fwrite(File2_updated, "LC-Meta_McKay_AllSamples_AgeGen10PCsChewYN_CommonSNPsOnly.csv_With_BETA-SE-OR_n_ConfIntervals_95_with_Effect_allele_frequecy_And_rsid.csv")
