###
# t-tests.R
###
setwd("~/Documents/OMI/Metabolon/")
df.counts = read.delim("./rawData/ECSC-08-15-from-Adam/clean-metabolite-counts.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
pathways.keys = read.delim("./rawData/ECSC-08-15-from-Adam/metabolite-names-key.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
df.pathways = merge(pathways.keys, df.counts)


Z_SCALE = FALSE

df.super.pathway = subset(df.pathways, SUPER_PATHWAY != "NA")

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

if (Z_SCALE){
  biochem.df[,-c(ncol(biochem.df)-1, ncol(biochem.df))] = scale(log10(1 + biochem.df[,-c(ncol(biochem.df)-1, ncol(biochem.df))]))
  biochem.hc = subset(biochem.df, sample.type == "HC")
  mean.hc = lapply(biochem.hc[,-c(ncol(biochem.hc)-1, ncol(biochem.hc))], FUN = function(x) mean(log10(x), na.rm = TRUE))
  sd.hc = lapply(biochem.hc[,-c(ncol(biochem.hc)-1, ncol(biochem.hc))], FUN = function(x) sd(log10(x), na.rm = TRUE))
  biochem.df[,-c(ncol(biochem.df)-1, ncol(biochem.df))] = (log10(biochem.df[,-c(ncol(biochem.df)-1, ncol(biochem.df))]) - mean.hc)/sd.hc
}

## turn to long format
library(reshape)
df.m = reshape::melt(biochem.df, id.vars = c("patient.sample", "sample.type"), factorsAsStrings = FALSE)
df.m$variable = as.character(df.m$variable)

if (Z_SCALE == FALSE){
  df.m$value = log10(1+df.m$value)
}


#########################################################
#
# ttests on HC vs Pre
#
#########################################################

df.pre.HC = subset(df.m, sample.type != "Post")
df.pre.HC$sample.type = factor(df.pre.HC$sample.type, levels = c("HC", "Pre"))


## analysis for all metabolites
library(dplyr)
library(broom)

my_ttest = function(input.data){
  if (length(which(is.na(input.data$value))) > 30){
    return(data.frame())
  } else {
    return(tidy(t.test(value ~ sample.type, alternative = "two.sided", data = input.data)))
  }
}


df.pre.HC %>%
  group_by(variable) %>%
  do(my_ttest(.)) %>%
  as.data.frame() -> ttest.results

pre.HC.ttest = merge(ttest.results, biochem.keys, by.x = "variable", by.y = "BIOCHEM.name", all.x = TRUE)


library(ggplot2)
library(plotly)
ggplot(pre.HC.ttest, aes(x = estimate, y = -log10(p.value), label = BIOCHEMICAL)) + geom_text() + ggtitle("HC vs Pre on raw log scaled data") -> p

#ggplotly(p, filename = "metabolon/HC-pre-raw-log", fileopt = "overwrite")
plotly_POST(p, filename = "metabolon/HC-pre-raw-log", fileopt = "overwrite")

#########################################################
#
# ttests on HC vs Post
#
#########################################################
df.post.HC = subset(df.m, sample.type != "Pre")
df.post.HC$sample.type = factor(df.post.HC$sample.type, levels = c("HC", "Post"))


df.post.HC %>%
  group_by(variable) %>%
  do(my_ttest(.)) %>%
  as.data.frame() -> ttest.results

post.HC.ttest = merge(ttest.results, biochem.keys, by.x = "variable", by.y = "BIOCHEM.name", all.x = TRUE)

ggplot(post.HC.ttest, aes(x = estimate, y = -log10(p.value), label = BIOCHEMICAL)) + geom_text() + ggtitle("HC vs Post on raw log scaled data") -> p
plotly_POST(p, filename = "metabolon/HC-post-raw-log", fileopt = "overwrite")


#########################################################
#
# ttests on Pre vs Post
#
#########################################################
df.pre.post = subset(df.m, sample.type != "HC")
df.pre.post$sample.type = factor(df.pre.post$sample.type, levels = c("Pre", "Post"))


df.pre.post %>%
  group_by(variable) %>%
  do(my_ttest(.)) %>%
  as.data.frame() -> ttest.results

pre.post.ttest = merge(ttest.results, biochem.keys, by.x = "variable", by.y = "BIOCHEM.name", all.x = TRUE)

ggplot(pre.post.ttest, aes(x = estimate, y = -log10(p.value), label = BIOCHEMICAL)) + geom_text() + ggtitle("Pre vs Post on raw log scaled data") -> p
plotly_POST(p, filename = "metabolon/pre-post-raw-log", fileopt = "overwrite")
