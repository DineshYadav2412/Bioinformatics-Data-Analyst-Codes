setwd("C:/Users/User/Downloads/Selected_Traits/ukb_b_19953")

library(data.table)

#-----------------------------------------------------
# Read GWAS file (exposure dataset)
#-----------------------------------------------------
result_g <- fread("ukb_b_19953.csv")
dim(result_g)
head(result_g)
cat(colnames(result_g),sep = "\n")
head(result_g)
#----------------------------------------------------------
#convert in P
#----------------------------------------------------------
result_g$P <- 10^(-result_g$LOGP)
write.csv(result_g, "ukb_b_19953_Pvalue.csv", row.names = FALSE)
head(result_g)
dim(result_g)