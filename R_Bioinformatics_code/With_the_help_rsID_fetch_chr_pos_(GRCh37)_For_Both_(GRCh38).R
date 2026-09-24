#=============================================================
# ADD CYTOBAND TO SNP POSITION FILES
# GRCh37 + GRCh38
#=============================================================

#-------------------------------------------------------------
# 1. LOAD PACKAGES
#-------------------------------------------------------------

library(data.table)
library(dplyr)


#-------------------------------------------------------------
# 2. SET WORKING DIRECTORY
#-------------------------------------------------------------

setwd("D:/COMMON_FOLDER_FOR _TEJAS_DINESH/tryal_cytoband")


#=============================================================
# 3. FILE NAMES
#=============================================================

# SNP position files
grch37_file <- "SNP_Position_GRCh37.csv"
grch38_file <- "SNP_Position_GRCh38.csv"

# Cytoband files
cyto37_file <- "cytoBand_hg19_GRCh37.txt"
cyto38_file <- "cytoBand_HG38.txt"

# Output files
output37 <- "SNP_Position_GRCh37_with_Cytoband.csv"
output38 <- "SNP_Position_GRCh38_with_Cytoband.csv"


#=============================================================
# 4. CHECK THAT FILES EXIST
#=============================================================

required_files <- c(
  grch37_file,
  grch38_file,
  cyto37_file,
  cyto38_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "The following file(s) were not found:",
      paste(missing_files, collapse = ", ")
    )
  )
}


#=============================================================
# 5. READ SNP POSITION FILES
#=============================================================

cat("\n=============================================\n")
cat("Reading SNP position files...\n")
cat("=============================================\n")

grch37 <- fread(
  grch37_file
)

grch38 <- fread(
  grch38_file
)


#=============================================================
# 6. READ CYTOBAND FILES
#=============================================================

cat("\n=============================================\n")
cat("Reading cytoband files...\n")
cat("=============================================\n")

cyto37 <- fread(
  cyto37_file,
  header = FALSE
)

cyto38 <- fread(
  cyto38_file,
  header = FALSE
)


#=============================================================
# 7. SET CYTOBAND COLUMN NAMES
#=============================================================

setnames(
  cyto37,
  c(
    "Chromosome",
    "Start",
    "End",
    "Cytoband",
    "Stain"
  )
)

setnames(
  cyto38,
  c(
    "Chromosome",
    "Start",
    "End",
    "Cytoband",
    "Stain"
  )
)


#=============================================================
# 8. PREPARE CYTOBAND DATA
#=============================================================

cyto37 <- cyto37 %>%
  mutate(
    Chromosome = as.character(Chromosome),
    Start = as.numeric(Start),
    End = as.numeric(End),
    Cytoband = as.character(Cytoband)
  )

cyto38 <- cyto38 %>%
  mutate(
    Chromosome = as.character(Chromosome),
    Start = as.numeric(Start),
    End = as.numeric(End),
    Cytoband = as.character(Cytoband)
  )


#=============================================================
# 9. FUNCTION TO ADD CYTOBAND
#=============================================================

add_cytoband <- function(snp_data, cyto_data) {
  
  #-----------------------------------------------------------
  # Make a copy so original data is not changed
  #-----------------------------------------------------------
  
  snp_data <- copy(snp_data)
  
  
  #-----------------------------------------------------------
  # Make sure required columns exist
  #-----------------------------------------------------------
  
  required_columns <- c(
    "rsID",
    "Chromosome",
    "Position"
  )
  
  missing_columns <- required_columns[
    !required_columns %in% names(snp_data)
  ]
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste(
        "Missing column(s) in SNP file:",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  
  #-----------------------------------------------------------
  # Convert chromosome to character
  #-----------------------------------------------------------
  
  snp_data$Chromosome <- as.character(
    snp_data$Chromosome
  )
  
  
  #-----------------------------------------------------------
  # Remove existing "chr" if present
  #-----------------------------------------------------------
  
  snp_data$Chromosome <- sub(
    "^chr",
    "",
    snp_data$Chromosome,
    ignore.case = TRUE
  )
  
  
  #-----------------------------------------------------------
  # Add "chr" because cytoband files use chr1, chr2, etc.
  #-----------------------------------------------------------
  
  snp_data$Chromosome <- paste0(
    "chr",
    snp_data$Chromosome
  )
  
  
  #-----------------------------------------------------------
  # Convert position to numeric
  #-----------------------------------------------------------
  
  snp_data$Position <- as.numeric(
    snp_data$Position
  )
  
  
  #-----------------------------------------------------------
  # Create Cytoband column
  #-----------------------------------------------------------
  
  snp_data$Cytoband <- NA_character_
  
  
  #===========================================================
  # Find cytoband for every SNP
  #===========================================================
  
  for (i in seq_len(nrow(snp_data))) {
    
    chr <- snp_data$Chromosome[i]
    pos <- snp_data$Position[i]
    
    
    #---------------------------------------------------------
    # Skip if chromosome or position is missing
    #---------------------------------------------------------
    
    if (
      is.na(chr) ||
      is.na(pos)
    ) {
      next
    }
    
    
    #---------------------------------------------------------
    # Find matching cytoband
    #
    # Start <= Position < End
    #---------------------------------------------------------
    
    hit <- cyto_data[
      cyto_data$Chromosome == chr &
        cyto_data$Start <= pos &
        pos < cyto_data$End,
    ]
    
    
    #---------------------------------------------------------
    # Assign cytoband
    #---------------------------------------------------------
    
    if (nrow(hit) > 0) {
      
      snp_data$Cytoband[i] <- hit$Cytoband[1]
      
    }
  }
  
  
  #-----------------------------------------------------------
  # Remove "chr" again from final SNP chromosome column
  #-----------------------------------------------------------
  
  snp_data$Chromosome <- sub(
    "^chr",
    "",
    snp_data$Chromosome
  )
  
  
  return(snp_data)
}


#=============================================================
# 10. ADD GRCh37 CYTOBANDS
#=============================================================

cat("\n=============================================\n")
cat("Adding GRCh37 cytobands...\n")
cat("=============================================\n")

grch37_final <- add_cytoband(
  snp_data = grch37,
  cyto_data = cyto37
)


#=============================================================
# 11. ADD GRCh38 CYTOBANDS
#=============================================================

cat("\n=============================================\n")
cat("Adding GRCh38 cytobands...\n")
cat("=============================================\n")

grch38_final <- add_cytoband(
  snp_data = grch38,
  cyto_data = cyto38
)


#=============================================================
# 12. SAVE GRCh37
#=============================================================

cat("\n=============================================\n")
cat("Saving GRCh37...\n")
cat("=============================================\n")

fwrite(
  grch37_final,
  output37
)


#=============================================================
# 13. SAVE GRCh38
#=============================================================

cat("\n=============================================\n")
cat("Saving GRCh38...\n")
cat("=============================================\n")

fwrite(
  grch38_final,
  output38
)


#=============================================================
# 14. SHOW GRCh37 RESULTS
#=============================================================

cat("\n=============================================\n")
cat("GRCh37 RESULTS\n")
cat("=============================================\n\n")

print(grch37_final)


#=============================================================
# 15. SHOW GRCh38 RESULTS
#=============================================================

cat("\n=============================================\n")
cat("GRCh38 RESULTS\n")
cat("=============================================\n\n")

print(grch38_final)


#=============================================================
# 16. SUMMARY
#=============================================================

cat("\n=============================================\n")
cat("SUMMARY\n")
cat("=============================================\n")

cat(
  "\nTotal GRCh37 SNPs:",
  nrow(grch37_final),
  "\n"
)

cat(
  "GRCh37 SNPs with cytoband:",
  sum(!is.na(grch37_final$Cytoband)),
  "\n"
)

cat(
  "GRCh37 SNPs without cytoband:",
  sum(is.na(grch37_final$Cytoband)),
  "\n"
)


cat(
  "\nTotal GRCh38 SNPs:",
  nrow(grch38_final),
  "\n"
)

cat(
  "GRCh38 SNPs with cytoband:",
  sum(!is.na(grch38_final$Cytoband)),
  "\n"
)

cat(
  "GRCh38 SNPs without cytoband:",
  sum(is.na(grch38_final$Cytoband)),
  "\n"
)


#=============================================================
# 17. SHOW SNPs WITH MISSING CYTOBAND
#=============================================================

cat("\n=============================================\n")
cat("GRCh37 SNPs WITHOUT CYTOBAND\n")
cat("=============================================\n")

print(
  grch37_final[
    is.na(grch37_final$Cytoband),
  ]
)


cat("\n=============================================\n")
cat("GRCh38 SNPs WITHOUT CYTOBAND\n")
cat("=============================================\n")

print(
  grch38_final[
    is.na(grch38_final$Cytoband),
  ]
)


#=============================================================
# 18. OUTPUT FILE LOCATIONS
#=============================================================

cat("\n=============================================\n")
cat("FILES CREATED\n")
cat("=============================================\n")

cat(
  "\nGRCh37:\n",
  file.path(getwd(), output37),
  "\n"
)

cat(
  "\nGRCh38:\n",
  file.path(getwd(), output38),
  "\n"
)

cat("\nDONE.\n")

