###
# heatmaps.R
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

USING_SUPER_PATHWAY = TRUE

if (USING_SUPER_PATHWAY){
  pathways.list = split(df.pathways, f = df.pathways$SUPER_PATHWAY)
} else{ #use sub_pathway and only for lipids
  lipids = subset(df.pathways, SUPER_PATHWAY == "Lipid")
  pathways.list = split(lipids, f = lipids$SUB_PATHWAY)
}

library(gplots)
print(names(pathways.list))
for (i in c(1:length(pathways.list))){
  subpathway = names(pathways.list)[i]

  print(paste("working on", subpathway))
  
  df.super.pathway = pathways.list[[subpathway]]
  
  df = df.super.pathway[,c(1,13:ncol(df.super.pathway))]
  row.names(df) = df$BIOCHEMICAL
  df = df[,-1]
  
  color.var = colnames(df)
  sample.type = sapply(color.var, FUN = function(x) strsplit(x, "CFS.")[[1]][2])
  sample.type[is.na(sample.type)] = "HC"
  sample.type = replace(sample.type, which(sample.type == "HC"), "black")
  sample.type = replace(sample.type, which(sample.type == "Pre"), "magenta")
  sample.type = replace(sample.type, which(sample.type == "Post"), "deepskyblue")
  pdf(paste0("./figures/ECSC-08-15-from-Adam/heatmaps/log-col-scaled/", subpathway, "-col-scaled-heatmap.pdf"))
  heatmap.2(as.matrix(log2(df)), trace = "none", ColSideColors = sample.type, scale = "column",
            margins = c(8,12), cexRow = .3, cexCol = .4)
  legend("topright",      
         legend = c("HC", "Pre", "Post"),
         col = c("black", "magenta", "deepskyblue"), 
         lty= 1,             
         lwd = 5,           
         cex=.7
  )
  dev.off()
}
