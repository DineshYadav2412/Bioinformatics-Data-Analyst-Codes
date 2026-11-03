library(data.table)

setwd("D:/COMMON_FOLDER_FOR_TEJAS_DINESH/COMMON_FILES_THAT_USE_MOST_OF_THE_TIME/LC_DATA_ALLELE_FREQUENCY/Final_Alignment")
library(data.table)

## ---------------------------------------------------------------
## Inputs
## ---------------------------------------------------------------
ind_file  <- "3_INDAI_DuplisCleaned_Logistic_All_10PCs_Age_Gen_ChewYN_OR.assoc.logistic.csv_With_BETA-SE-OR_n_ConfIntervals_95.csv"
mck_file  <- "ebi-a-GCST004747_COLUMNS_FIXED_FOR_META.csv_With_BETA-SE-OR_n_ConfIntervals_95.csv"
N_IND <- 5338
N_MCK <- 9859

ind <- fread(ind_file)
mck <- fread(mck_file)

## Keep only the columns needed (Indian: A1 is the effect allele, A2 the other)
ind <- ind[, .(CHR = as.character(CHR), BP = as.integer(BP), IND_SNP = SNP,
               IND_A1 = toupper(A1), IND_A2 = toupper(A2_MAJORALLELE),
               IND_EAF = A1_MINORALLELEFREQ,
               IND_BETA = BETA, IND_OR = OR, IND_P = P)]

mck <- mck[, .(CHR = as.character(CHR), BP = as.integer(BP), MCK_SNP = SNP,
               MCK_A1 = toupper(A1_MINORALLELE), MCK_A2 = toupper(A2_MAJORALLELE),
               MCK_EAF = A1_MINORALLELEFREQ,
               MCK_BETA = BETA, MCK_SE = SE, MCK_OR = OR, MCK_P = P,
               MCK_BETA_L = BETA_CONFINT_95_LOWER, MCK_BETA_U = BETA_CONFINT_95_UPPER,
               MCK_OR_L = OR_CONFINT_95_LOWER, MCK_OR_U = OR_CONFINT_95_UPPER)]

## Duplicates are kept as they are (no removal)
cat("Duplicated CHR:BP - Indian:", sum(duplicated(ind[, .(CHR, BP)])),
    " McKay:", sum(duplicated(mck[, .(CHR, BP)])), "\n")

## ---------------------------------------------------------------
## Step 1: Common SNPs
## ---------------------------------------------------------------
## Merge on CHR + BP (rsIDs can differ between files).
## allow.cartesian = TRUE: a position duplicated in either file gives one row per pairing
common <- merge(ind, mck, by = c("CHR", "BP"), allow.cartesian = TRUE)
cat("Common SNPs:", nrow(common), "\n")

## ---------------------------------------------------------------
## Step 2: Matched SNPs (same allele pair, in either order)
## ---------------------------------------------------------------
common[, same_order := (IND_A1 == MCK_A1 & IND_A2 == MCK_A2)]
common[, swapped    := (IND_A1 == MCK_A2 & IND_A2 == MCK_A1)]
matched <- common[same_order | swapped]
cat("Matched SNPs:", nrow(matched), "\n")

## Flag strand-ambiguous (palindromic) SNPs: A/T and C/G.
## For these, a same-order or swapped match cannot be confirmed as the same strand
## without looking at frequencies. Flagged here; set DROP_PALINDROMIC <- TRUE to remove.
matched[, palindromic := paste0(pmin(IND_A1, IND_A2), pmax(IND_A1, IND_A2)) %in% c("AT", "CG")]
cat("Palindromic among matched:", sum(matched$palindromic), "\n")
DROP_PALINDROMIC <- FALSE
if (DROP_PALINDROMIC) matched <- matched[palindromic == FALSE]

## ---------------------------------------------------------------
## Step 4: Align McKay to Indian alleles
## ---------------------------------------------------------------
aligned <- copy(matched)
aligned[, flipped := swapped]

fl <- aligned$flipped
## NOTE: copy() is essential. data.table's := edits columns in place, so a plain
## `tmp <- aligned$MCK_A2` would be silently modified when MCK_A2 is overwritten.
## swap alleles
tmp <- copy(aligned$MCK_A2)
aligned[fl, MCK_A2 := MCK_A1]
aligned[fl, MCK_A1 := tmp[fl]]
## flip effect sizes / frequency
aligned[fl, MCK_BETA := -1 * MCK_BETA]
aligned[fl, MCK_OR   := 1 / MCK_OR]
aligned[fl, MCK_EAF  := 1 - MCK_EAF]
## flip confidence intervals (bounds swap as well)
bl <- copy(aligned$MCK_BETA_L); bu <- copy(aligned$MCK_BETA_U)
ol <- copy(aligned$MCK_OR_L);   ou <- copy(aligned$MCK_OR_U)
aligned[fl, MCK_BETA_L := -bu[fl]]
aligned[fl, MCK_BETA_U := -bl[fl]]
aligned[fl, MCK_OR_L := 1 / ou[fl]]
aligned[fl, MCK_OR_U := 1 / ol[fl]]

## sanity check: alleles must now agree everywhere
stopifnot(all(aligned$IND_A1 == aligned$MCK_A1), all(aligned$IND_A2 == aligned$MCK_A2))
cat("Flipped:", sum(fl), " Unchanged:", sum(!fl), "\n")

## ---------------------------------------------------------------
## Steps 6-8: Effect allele sample counts and Meta EAF
## ---------------------------------------------------------------
aligned[, EASC_IND := N_IND * IND_EAF]
aligned[, EASC_MCK := N_MCK * MCK_EAF]
aligned[, MEASC    := EASC_IND + EASC_MCK]
aligned[, META_EAF := MEASC / (N_IND + N_MCK)]

## ---------------------------------------------------------------
## Output (effect allele = Indian A1, aligned across both studies)
## ---------------------------------------------------------------
out <- aligned[, .(CHR, BP, IND_SNP, MCK_SNP,
                   EFFECT_ALLELE = IND_A1, OTHER_ALLELE = IND_A2,
                   IND_EAF, MCK_EAF_ALIGNED = MCK_EAF,
                   EASC_IND, EASC_MCK, MEASC, META_EAF,
                   IND_BETA, MCK_BETA_ALIGNED = MCK_BETA,
                   IND_OR, MCK_OR_ALIGNED = MCK_OR,
                   IND_P, MCK_P, palindromic, flipped)]

fwrite(out, "Final_Aligned_Meta_EAF_Indian_McKay.csv")
cat("Written:", nrow(out), "SNPs to Aligned_Meta_EAF_Indian_McKay.csv\n")
head(out)
