###
# HC-variance-analysis.R
###
setwd("~/Documents/OMI/Metabolon/")

df.counts = read.delim("./rawData/ECSC-08-15-from-Adam/clean-metabolite-counts.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
pathways.keys = read.delim("./rawData/ECSC-08-15-from-Adam/metabolite-names-key.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
df.pathways = merge(pathways.keys, df.counts)

USING_SUPER_PATHWAY = TRUE

Z_SCALE = TRUE

if (USING_SUPER_PATHWAY){
  pathways.list = split(df.pathways, f = df.pathways$SUPER_PATHWAY)
} else{ #use sub_pathway and only for lipids
  lipids = subset(df.pathways, SUPER_PATHWAY == "Lipid")
  pathways.list = split(lipids, f = lipids$SUB_PATHWAY)
}

print(names(pathways.list))

for (i in c(1:length(pathways.list))){
  subpathway = names(pathways.list)[i]
  print(paste("working on", subpathway))

  df.super.pathway = pathways.list[[subpathway]]
  
  if (nrow(df.super.pathway) >= 10) {
    num.metabolites = 10
  } else{
    num.metabolites = nrow(df.super.pathway)
  }
  
  
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
    #biochem.df[,-c(ncol(biochem.df)-1, ncol(biochem.df))] = scale(log10(1 + biochem.df[,-c(ncol(biochem.df)-1, ncol(biochem.df))]))
    biochem.hc = subset(biochem.df, sample.type == "HC")
    mean.hc = lapply(biochem.hc[,-c(ncol(biochem.hc)-1, ncol(biochem.hc))], FUN = function(x) mean(log10(x), na.rm = TRUE))
    sd.hc = lapply(biochem.hc[,-c(ncol(biochem.hc)-1, ncol(biochem.hc))], FUN = function(x) sd(log10(x), na.rm = TRUE))
    biochem.df[,-c(ncol(biochem.df)-1, ncol(biochem.df))] = (log10(biochem.df[,-c(ncol(biochem.df)-1, ncol(biochem.df))]) - mean.hc)/sd.hc
  }
  
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
  if (Z_SCALE == FALSE){
    df.m$value = log10(1+value)
  }
  
  df.m %>%
    group_by(sample.type, variable) %>%
    summarise(median.metabolite.na.rm = median(value, na.rm = TRUE),
              mad.metabolite.na.rm = mad(value, na.rm = TRUE),
              median.metabolite = median(value, na.rm = FALSE),
              mad.metabolite = mad(value, na.rm = FALSE)) %>%
    as.data.frame() -> df.stats
  
  hc.stats = subset(df.stats, sample.type == "HC")
  hc.stats.ordered = hc.stats[order(hc.stats$mad.metabolite), ]
  hc.least.variant = hc.stats.ordered$variable[c(1:num.metabolites)]
  
  #########################################################
  #
  #  make slope graphs
  #
  #########################################################
  ## remove helath and only keep least variance matabolites amongst HC
  top10.df = subset(df.m, variable %in% hc.least.variant)
  top10.df = merge(top10.df, biochem.keys, by.x = "variable", by.y = "BIOCHEM.name", all.x = TRUE)
  top10.df$sample.type = factor(top10.df$sample.type, levels = c("HC", "Pre", "Post"))
  ggplot(top10.df, aes(factor(BIOCHEMICAL), value)) + 
    geom_boxplot(aes(fill = sample.type)) + xlab("") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 24)) + ggtitle(subpathway) -> p
  
  ggsave(filename = paste0("./figures/ECSC-08-15-from-Adam/superpathway-metabolties/z-scaled/least-variance/" , subpathway, "-least-variant-HC-boxplot.png"), 
         plot = p, width = 26, height = 18)
  
  
  
  cfs.df = subset(top10.df, sample.type != "HC")
  cfs.df$patient.sample = sapply(cfs.df$patient.sample, FUN = function(x) strsplit(x, split = "\\.CFS")[[1]][1])
  
  cfs.df$sample.type = factor(cfs.df$sample.type, levels = c("Pre", "Post"))
  
  loner.patients = c("ID.9", "ID.78", "ID.311", "ID.730")
  ggplot(subset(cfs.df, !(patient.sample %in% loner.patients)), aes(x = sample.type, y = value, colour = BIOCHEMICAL, group= BIOCHEMICAL)) +
    geom_line() + geom_point(shape=21, fill = "white") + xlab("")  + 
    facet_wrap(~patient.sample) + ggtitle(subpathway) + theme(legend.text=element_text(size=24)) -> p
  
  ggsave(filename = paste0("./figures/ECSC-08-15-from-Adam/superpathway-metabolties/z-scaled/least-variance/" , subpathway, "-least-variant-HC-slopegraph.png"), 
         plot = p, width = 26, height = 18)

}