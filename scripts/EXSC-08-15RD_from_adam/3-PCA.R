###
# PCA.R
###
setwd("~/Documents/OMI/Metabolon/")

df.counts = read.delim("./rawData/ECSC-08-15-from-Adam/clean-metabolite-counts.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
#########################################################
#
# !!!!!!!!!!!!!!!!!!REMOVE NA'S HERE!!!!!!!!!!!!!!!!!!!!!
#
#########################################################
#df.counts = na.omit(df.counts)
df.counts[,-1] = lapply(df.counts[,-1], as.numeric)
for(i in c(2:ncol(df.counts))){
  df.counts[is.na(df.counts[,i]),i ] = min(df.counts[,i], na.rm = TRUE)
}


pathways.keys = read.delim("./rawData/ECSC-08-15-from-Adam/metabolite-names-key.tsv", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
df.pathways = merge(pathways.keys, df.counts)

USING_SUPER_PATHWAY = FALSE

if (USING_SUPER_PATHWAY){
  pathways.list = split(df.pathways, f = df.pathways$SUPER_PATHWAY)
} else{ #use sub_pathway and only for lipids
  lipids = subset(df.pathways, SUPER_PATHWAY == "Lipid")
  pathways.list = split(lipids, f = lipids$SUB_PATHWAY)
}


print(names(pathways.list))
for (i in c(1:length(pathways.list))){
  subpathway = names(pathways.list)[i]
  if (subpathway == "Energy"){
    next
  }
  print(paste("working on", subpathway))
  
  df.super.pathway = pathways.list[[subpathway]]
  
  df = df.super.pathway[,c(1,13:ncol(df.super.pathway))]
  biochem.keys = data.frame(BIOCHEMICAL = df$BIOCHEMICAL, BIOCHEM.name = make.names(df$BIOCHEMICAL))
  
  if(nrow(df) < 6){
    next
  }
  
  num.loadings.per.pc = 3
  num.pcs = 2
  
  
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
  
  add.patient.data = function(x.df){
    x.df$patient.sample = make.names(row.names(x.df))
    x.df$sample.type = sapply(x.df$patient.sample, FUN = function(x) strsplit(x, "CFS.")[[1]][2])
    x.df$sample.type[is.na(x.df$sample.type)] = "HC"
    return(x.df)
  }
  
  
  #########################################################
  #
  # PCA
  #
  #########################################################
  pca.result = prcomp(log(biochem.df), scale = FALSE)
  
  top.loadings = lapply(as.data.frame(pca.result$rotation), FUN = function(x) row.names(pca.result$rotation)[order(abs(x), decreasing = TRUE)][c(1:num.loadings.per.pc)])
  pca.df = add.patient.data(as.data.frame(pca.result$x))
  ggplot(pca.df, aes(x = PC1, y = PC2, colour = sample.type)) + geom_point() + ggtitle("PCA") + ggtitle(subpathway) -> p
  ggsave(filename = paste0("./figures/ECSC-08-15-from-Adam/lipids-subpathways/PCA/" , subpathway, "-min-imputation.png"), 
         plot = p, width = 12, height = 9)
  #########################################################
  #
  # MDS
  #
  #########################################################
  library(MASS)
  d = dist(log(biochem.df), method = "maximum")
  fit <- cmdscale(d, k=2, eig = TRUE) # k is the number of dim
  mds.df = add.patient.data(as.data.frame(fit$points))
  ggplot(mds.df, aes(x = V1, y = V2, colour = sample.type)) + geom_point() + ggtitle("MDS") -> p
  
  #########################################################
  #
  #  make slope graphs
  #
  #########################################################
  ## turn to long format
  library(reshape)
  df.m = reshape::melt(add.patient.data(biochem.df), id.vars = c("patient.sample", "sample.type"), factorsAsStrings = FALSE)
  df.m$variable = as.character(df.m$variable)
  
  top10.df = subset(df.m, variable %in% unique(unname(unlist(head(top.loadings, num.pcs)))))
  top10.df = merge(top10.df, biochem.keys, by.x = "variable", by.y = "BIOCHEM.name", all.x = TRUE)
  top10.df$sample.type = factor(top10.df$sample.type, levels = c("HC", "Pre", "Post"))
  
  ggplot(top10.df, aes(factor(BIOCHEMICAL), log10(value))) + 
    geom_boxplot(aes(fill = sample.type)) + xlab("") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14)) + ggtitle(subpathway) -> p
  
  ggsave(filename = paste0("./figures/ECSC-08-15-from-Adam/lipids-subpathways/PCA/" , subpathway, "-min-imputation-boxplots.png"), 
         plot = p, width = 26, height = 18)
  
  cfs.df = subset(top10.df, sample.type != "HC")
  cfs.df$patient.sample = sapply(cfs.df$patient.sample, FUN = function(x) strsplit(x, split = "\\.CFS")[[1]][1])
  
  cfs.df$sample.type = factor(cfs.df$sample.type, levels = c("Pre", "Post"))
  
  loner.patients = c("ID.9", "ID.78", "ID.311", "ID.730")
  ggplot(subset(cfs.df, !(patient.sample %in% loner.patients)), aes(x = sample.type, y = value, colour = BIOCHEMICAL, group= BIOCHEMICAL)) +
    geom_line() + geom_point(shape=21, fill = "white") + xlab("")  + facet_wrap(~patient.sample) + ylab("") + scale_y_log10() -> p
  
  ggsave(filename = paste0("./figures/ECSC-08-15-from-Adam/lipids-subpathways/PCA/" , subpathway, "-min-imputation-slopegraph.png"), 
         plot = p, width = 26, height = 18)
}
