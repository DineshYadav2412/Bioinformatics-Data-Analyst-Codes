#=============================================
# Get SNP Positions (GRCh37)
#=============================================

library(httr)
library(jsonlite)
library(dplyr)

# Set working directory
setwd("D:/COMMON_FOLDER_FOR _TEJAS_DINESH/tryal_cytoband")

#---------------------------------------
# SNP list
#---------------------------------------

snps <- c(
  "rs2291428","rs2290846","rs1800961","rs601338","rs708686",
  "rs34851490","rs1169288","rs13280055","rs174567","rs11012737",
  "rs2469991","rs1935","rs17240268","rs12004","rs55971546",
  "rs11641445","rs17138478","rs2292553","rs12968116","rs11887534",
  "rs212100","rs12633863","rs4148808","rs6471717","rs686030",
  "rs1260326","rs2070959","rs28929474","rs56398830",
  "rs756082276","rs756935975","rs45575636"
)

#---------------------------------------
# Function to query GRCh37 Ensembl
#---------------------------------------

get_snp <- function(rsid){
  
  url <- paste0(
    "https://grch37.rest.ensembl.org/variation/human/",
    rsid
  )
  
  r <- GET(
    url,
    add_headers("Content-Type" = "application/json")
  )
  
  if(status_code(r) != 200){
    
    return(data.frame(
      rsID = rsid,
      Chromosome = NA,
      Position = NA,
      Cytoband = NA,
      stringsAsFactors = FALSE
    ))
    
  }
  
  x <- fromJSON(content(r, "text", encoding = "UTF-8"))
  
  chromosome <- NA
  position <- NA
  
  if(!is.null(x$mappings) && nrow(x$mappings) > 0){
    
    chromosome <- x$mappings$seq_region_name[1]
    position <- x$mappings$start[1]
    
  }
  
  data.frame(
    rsID = rsid,
    Chromosome = chromosome,
    Position = position,
    Cytoband = NA,
    stringsAsFactors = FALSE
  )
}

#---------------------------------------
# Query all SNPs
#---------------------------------------

result <- bind_rows(lapply(snps, get_snp))

#---------------------------------------
# Save output
#---------------------------------------

write.csv(
  result,
  "SNP_Position_GRCh37.csv",
  row.names = FALSE
)

print(result)
