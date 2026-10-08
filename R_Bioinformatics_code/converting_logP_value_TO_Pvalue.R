
setwd("C:/Users/Sharayu/Desktop/New_Analysis/bbj_a_73")

library(data.table)

#-----------------------------------------------------
# Read GWAS file (exposure dataset)
#-----------------------------------------------------
result_g <- fread("bbj_a_73.csv")
dim(result_g)
cat(colnames(result_g),sep = "\n")
head(result_g)
#----------------------------------------------------------
#convert in P
#----------------------------------------------------------
result_g$P <- 10^(-result_g$LOGP)
write.csv(result_g, "bbj_a_73_Pvalue.csv", row.names = FALSE)
head(result_g)
dim(result_g)


data_new <- fread("bbj_a_73_Pvalue.csv")
dim(data_new)
names(data_new)

data_new <- data_new[data_new$P < 5e-8, ]
dim(data_new)
write.csv(data_new,"bbj_a_73_Pvalue_5e_8.csv")


