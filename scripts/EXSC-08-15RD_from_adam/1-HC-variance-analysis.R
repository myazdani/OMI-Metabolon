###
# HC-variance-analysis.R
###
setwd("~/Documents/OMI/Metabolon/")

df = read.delim("./rawData/ECSC-08-15-from-Adam/clean-metabolite-counts.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
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

## turn to long format
library(reshape)
df.m = reshape::melt(biochem.df, id.vars = c("patient.sample", "sample.type"), factorsAsStrings = FALSE)
df.m$variable = as.character(df.m$variable)

#########################################################
#
#  count number of NA's
#
#########################################################
library(dplyr)
df.m %>%
  group_by(patient.sample, sample.type) %>%
  summarise(num.missing = length(which(is.na(value)))) %>%
  as.data.frame() -> na.stats

na.stats$patient.sample = paste0(na.stats$sample.type, ".", na.stats$patient.sample)
library(ggplot2)
ggplot(na.stats, aes(x = num.missing, y = patient.sample, colour = sample.type)) + geom_point()

#########################################################
#
#  compute variances
#
#########################################################

df.m %>%
  group_by(sample.type, variable) %>%
  summarise(median.metabolite.na.rm = median(log10(1+value), na.rm = TRUE),
            mad.metabolite.na.rm = mad(log10(1+value), na.rm = TRUE),
            median.metabolite = median(log10(1+value), na.rm = FALSE),
            mad.metabolite = mad(log10(1+value), na.rm = FALSE)) %>%
  as.data.frame() -> df.stats

hc.stats = subset(df.stats, sample.type == "HC")
hc.stats.ordered = hc.stats[order(hc.stats$mad.metabolite), ]
hc.least.variant = hc.stats.ordered$variable[c(1:50)]

#########################################################
#
#  make slope graphs
#
#########################################################
## remove helath and only keep least variance matabolites amongst HC
top10.df = subset(df.m, variable %in% hc.least.variant)
top10.df = merge(top10.df, biochem.keys, by.x = "variable", by.y = "BIOCHEM.name", all.x = TRUE)
top10.df$sample.type = factor(top10.df$sample.type, levels = c("HC", "Pre", "Post"))
ggplot(top10.df, aes(factor(BIOCHEMICAL), log10(value))) + geom_boxplot(aes(fill = sample.type)) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1))



cfs.df = subset(top10.df, sample.type != "HC")
cfs.df$patient.sample = sapply(cfs.df$patient.sample, FUN = function(x) strsplit(x, split = "\\.CFS")[[1]][1])

cfs.df$sample.type = factor(cfs.df$sample.type, levels = c("Pre", "Post"))

loner.patients = c("ID.9", "ID.78", "ID.311", "ID.730")
ggplot(subset(cfs.df, !(patient.sample %in% loner.patients)), aes(x = sample.type, y = value, colour = BIOCHEMICAL, group= BIOCHEMICAL)) +
  geom_line() + geom_point(shape=21, fill = "white") + xlab("")  + facet_wrap(~patient.sample) + ylab("") -> p