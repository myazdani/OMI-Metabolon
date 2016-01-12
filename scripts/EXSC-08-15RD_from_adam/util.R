###
# HC-variance-analysis.R
###
setwd("~/Documents/OMI/Metabolon/")

df.counts = read.delim("./rawData/ECSC-08-15-from-Adam/clean-metabolite-counts.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
pathways.keys = read.delim("./rawData/ECSC-08-15-from-Adam/metabolite-names-key.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
df.pathways = merge(pathways.keys, df.counts)

df.super.pathway = df.pathways
df = df.super.pathway[,c(1,13:ncol(df.super.pathway))]

biochem.keys = data.frame(BIOCHEMICAL = df$BIOCHEMICAL, BIOCHEM.name = make.names(df$BIOCHEMICAL))

## Tranpose data
# first remember the names
row.names(df) = df$BIOCHEMICAL
df = df[,-1]
n <- row.names(df)
# transpose all but the first column (name)
biochem.df <- as.data.frame(t(df), stringsAsFactors = FALSE)
colnames(biochem.df) <- n
names(biochem.df) = make.names(names(biochem.df))

biochem.df[,c(1:ncol(biochem.df))] = lapply(biochem.df, as.numeric)

biochem.df$patient.sample = make.names(row.names(biochem.df))
biochem.df$sample.type = sapply(biochem.df$patient.sample, FUN = function(x) strsplit(x, "CFS.")[[1]][2])
biochem.df$sample.type[is.na(biochem.df$sample.type)] = "HC"

biochem.df.pre = subset(biochem.df, sample.type != "Post")
biochem.df.post = subset(biochem.df, sample.type != "Pre")

biochem.df.pre[,-c(ncol(biochem.df)-1, ncol(biochem.df))] = scale(log(biochem.df[,-c(ncol(biochem.df)-1, ncol(biochem.df))]))

## turn to long format
library(reshape)
df.m = reshape::melt(biochem.df, id.vars = c("patient.sample", "sample.type"), factorsAsStrings = FALSE)
df.m$variable = as.character(df.m$variable)