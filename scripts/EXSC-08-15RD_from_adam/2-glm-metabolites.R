###
# HC-glm-analysis.R
###
setwd("~/Documents/OMI/Metabolon/")
df.counts = read.delim("./rawData/ECSC-08-15-from-Adam/clean-metabolite-counts.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
pathways.keys = read.delim("./rawData/ECSC-08-15-from-Adam/metabolite-names-key.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
df.pathways = merge(pathways.keys, df.counts)
pathways.list = split(df.pathways, f = df.pathways$SUPER_PATHWAY)

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
  
  ## turn to long format
  library(reshape)
  df.m = reshape::melt(biochem.df, id.vars = c("patient.sample", "sample.type"), factorsAsStrings = FALSE)
  df.m$variable = as.character(df.m$variable)
  
  #########################################################
  #
  # run glm per metabolites
  #
  #########################################################
  
  df.pre.HC = subset(df.m, sample.type != "Post")
  df.pre.HC$sample.type = factor(df.pre.HC$sample.type, levels = c("HC", "Pre"))
  
  
  ## analysis for all metabolites
  library(dplyr)
  library(broom)
  na.omit(df.pre.HC) %>%
    group_by(variable) %>%
    do(tidy(glm(sample.type ~ value, family = "binomial", .))) %>%
    as.data.frame() -> linear_classifiers
  
  
  na.omit(df.pre.HC) %>%
    group_by(variable) %>%
    do(glance(glm(sample.type ~ value, family = "binomial", .))) %>%
    as.data.frame() -> linear_classifiers.errors
  
  intercepts = subset(linear_classifiers, term == "(Intercept)")
  intercepts.sorted = intercepts[order(intercepts$p.value), ]
  slopes = subset(linear_classifiers, term == "value")
  slopes.sorted = slopes[order(slopes$p.value), ]
  aic.sorted = linear_classifiers.errors$variable[order(linear_classifiers.errors$AIC)[c(1:num.metabolites)]]
  
  
  top10.df = subset(df.m, variable %in% intercepts.sorted$variable[c(1:num.metabolites)])
  #top10.df = subset(df.m, variable %in% aic.sorted)
  
  #########################################################
  #
  #  make slope graphs
  #
  #########################################################
  ## remove helath and only keep least variance matabolites amongst HC
  top10.df = merge(top10.df, biochem.keys, by.x = "variable", by.y = "BIOCHEM.name", all.x = TRUE)
  top10.df$sample.type = factor(top10.df$sample.type, levels = c("HC", "Pre", "Post"))
  
  ggplot(top10.df, aes(factor(BIOCHEMICAL), log10(value))) + 
    geom_boxplot(aes(fill = sample.type)) + xlab("") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14)) + ggtitle(subpathway) -> p
  
  ggsave(filename = paste0("./figures/ECSC-08-15-from-Adam/superpathway-metabolties/glm-intercept/" , subpathway, "-HC-vs-pre-glm-boxplots.png"), 
         plot = p, width = 26, height = 18)
  
  
  cfs.df = subset(top10.df, sample.type != "HC")
  cfs.df$patient.sample = sapply(cfs.df$patient.sample, FUN = function(x) strsplit(x, split = "\\.CFS")[[1]][1])
  
  cfs.df$sample.type = factor(cfs.df$sample.type, levels = c("Pre", "Post"))
  
  loner.patients = c("ID.9", "ID.78", "ID.311", "ID.730")
  ggplot(subset(cfs.df, !(patient.sample %in% loner.patients)), aes(x = sample.type, y = value, colour = BIOCHEMICAL, group= BIOCHEMICAL)) +
    geom_line() + geom_point(shape=21, fill = "white") + xlab("")  + facet_wrap(~patient.sample) + ylab("") + scale_y_log10() -> p
  
  ggsave(filename = paste0("./figures/ECSC-08-15-from-Adam/superpathway-metabolties/glm-intercept/" , subpathway, "-HC-vs-pre-glm-slopegraph.png"), 
         plot = p, width = 26, height = 18)

}