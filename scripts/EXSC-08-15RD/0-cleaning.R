###
# cleaning.R
###

setwd("~/Documents/OMI/Metabolon/")
library(readxl)
df = read.csv("./rawData/EXSC-08-15RD/clean-CSV/EXSC-08-15RD-Xenobiotic-and-Unknown-REMOVED.csv", header = TRUE, stringsAsFactors = FALSE)




patients.df = df[,grep("REFERENCE.VS", names(df))]
row.names(patients.df) = df$Biochemical.Name



###
#
#  TRANSPOSE DATA
#
### 
# first remember the names
n <- row.names(patients.df)
# transpose all but the first column (name)
biochem.df <- as.data.frame(t(patients.df))
colnames(biochem.df) <- n
names(biochem.df) = make.names(names(biochem.df))

biochem.df$sample = sapply(row.names(biochem.df), FUN = function(x) strsplit(x, split = "REFERENCE.VS.")[[1]][2])

###
#
# MERGE WITH META DATA
#
###
df.keys = as.data.frame(read_excel("./rawData/EXSC-08-15RD/20150722_Metabolon_CFS_HC_Study_001.xlsx", sheet = 1, col_names = TRUE))
names(df.keys) = make.names(names(df.keys))
df.keys$Tube.ID.. = toupper(df.keys$Tube.ID..)

df.cfs = merge(df.keys, biochem.df, by.x = "Tube.ID..", by.y = "sample", all = FALSE)

#df.cfs[is.na(df.cfs)] = 0
for (i in which(sapply(df.cfs, is.numeric))) {
  df.cfs[is.na(df.cfs[, i]), i] <- mean(df.cfs[, i],  na.rm = TRUE)
}
df.cfs.clean = df.cfs[,-which(sapply(df.cfs,anyNA))]
pca.result <- prcomp(df.cfs.clean[,-c(1:9)])

pca.df = as.data.frame(pca.result$x, scale. = FALSE)
pca.df$cohort = df.cfs$Cohorts
pca.df$sample = df.cfs$Tube.ID..
pca.df$gender = df.cfs$Gender
ggplot(pca.df, aes(x = PC1, y = PC2, colour = cohort)) + geom_point()
