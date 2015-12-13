##
## rough-PCA.R
##

setwd("~/Documents/OMI/Metabolon/")
library(readxl)

df = as.data.frame(read_excel("./processedData/Metabolon_Statistics_Report_2 X006-00-00PM__R-friendly.xlsx", sheet = 1, col_names = TRUE))
pathway.keys = as.data.frame(read_excel("./processedData/BiochemicalKeys-2.xlsx", sheet = 1, col_names = TRUE))

names(pathway.keys) = make.names(names(pathway.keys))
names(df) = make.names(names(df))
##
## remove pyruvate since it is duplicated for no apparent reason!
##
df[which(df$Biochemical.Name %in% "pyruvate"),]

##
## record 179 is all missing so we remove
##

df = df[-179, ]
pathway.keys = unique(pathway.keys)
df.clean = merge(pathway.keys, df[,-c(1,2)], by = "Biochemical.Name")

patients.df = df.clean[,grep("REFERENCE.VS", names(df.clean))]

row.names(patients.df) = df.clean$Biochemical.Name

# first remember the names
n <- row.names(patients.df)

# transpose all but the first column (name)
biochem.df <- as.data.frame(t(patients.df))
colnames(biochem.df) <- n

names(biochem.df) = make.names(names(biochem.df))

biochem.df$subject = sapply(row.names(biochem.df), FUN = function(x) strsplit(x, split = "002_")[[1]][2])

ggplot(biochem.df, aes(x = xanthine, y = taurodeoxycholate, label = subject)) + geom_text() + geom_hline(yintercept=0)  + geom_vline(xintercept=0)
ggplot(biochem.df, aes(x = sphingomyelin..d18.1.20.0..d16.1.22.0.., y = sphingomyelin..d18.1.21.0..d17.1.22.0..d16.1.23.0.., label = subject)) + geom_text() + geom_hline(yintercept=0) + geom_vline(xintercept=0)

ggplot(biochem.df, aes(x = xanthine, y = sphingomyelin..d18.1.21.0..d17.1.22.0..d16.1.23.0.., label = subject)) + geom_text() + geom_hline(yintercept=0) + geom_vline(xintercept=0)

ggplot(biochem.df, aes(x = kynurenine, y = histidine, label = subject)) + geom_text() + geom_hline(yintercept=0) + geom_vline(xintercept=0)
ggplot(biochem.df, aes(x = arachidonate..20.4n6., y = X1..1.enyl.palmitoyl..2.arachidonoyl.GPC..P.16.0.20.4.., label = subject)) + geom_text() + geom_hline(yintercept=0) + geom_vline(xintercept=0)
##
# REPLACE MISSING VALUES WITH MINIMUM OF EACH COLUMN

for(i in 1:ncol(biochem.df)){
  if (length(biochem.df[is.na(biochem.df[,i]), i])>0){
    biochem.df[is.na(biochem.df[,i]), i] <- min(biochem.df[,i], na.rm = TRUE)
  }
}

library(checkmate)

biochem.clean.df = biochem.df[,-which(sapply(biochem.df,anyInfinite))]
pca.result <- prcomp(na.omit(biochem.clean.df))

top.loadings = lapply(as.data.frame(pca.result$rotation), FUN = function(x) row.names(pca.result$rotation)[order(abs(x), decreasing = TRUE)][c(1:5)])

pca.df = as.data.frame(pca.result$x)
