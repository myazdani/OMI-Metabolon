###
# merge-HZ-z-patient-z.R
###
setwd("~/Documents/OMI/Metabolon/")

library(readxl)
library(gplots)
HC.z = read_excel("./rawData/FromBryan/HC Z Scores.xlsx")
names(HC.z) = make.names(names(HC.z))

patient.z = read_excel("./rawData/ECSC-08-15-from-Adam/metabolites-z-scores-sheet-excel.xlsx")
names(patient.z) = make.names(names(patient.z))

# patient.df = patient.z[,-1]
# row.names(patient.df) = row.names(patient.z)
# 
# sample.type = colnames(patient.df)
# sample.type = replace(sample.type, grep("Pre", sample.type), "magenta")
# sample.type = replace(sample.type, grep("Post", sample.type), "deepskyblue")
# for(i in c(1:ncol(patient.df))){
#   patient.df[is.na(patient.df[,i]),i ] = mean(patient.df[,i], na.rm = TRUE)
# }
# pdf("Pre-Post-z-scale-heatmap.pdf")
# heatmap.2(as.matrix(patient.df), trace = "none", ColSideColors = sample.type,
#           margins = c(8,12), cexRow = .3, cexCol = .4)
# dev.off()

hc.names = HC.z$COMPOUND
p.names = patient.z$Biochemical.Name

clean.up.func = function(x){
  #removed.stuff = gsub(" ", "", x, fixed = TRUE)
  removed.stuff = gsub("\\s", "", x)
  return(make.names(tolower(removed.stuff), unique = TRUE))
}

length(intersect(clean.up.func(hc.names), clean.up.func(p.names)))

HC.z$HC.clean.compound.name = clean.up.func(HC.z$COMPOUND)
patient.z$patient.clean.compound.name = clean.up.func(patient.z$Biochemical.Name)

inner.join = merge(HC.z, patient.z, by.x ="HC.clean.compound.name", by.y = "patient.clean.compound.name")
outer.join = merge(HC.z, patient.z, by.x ="HC.clean.compound.name", by.y = "patient.clean.compound.name", all  = TRUE)

###
### heatmap of inner join
###

row.names(inner.join) = inner.join$Biochemical.Name
HC.df = inner.join[,grep("Meta", names(inner.join))]
pre.df = inner.join[,grep("Pre", names(inner.join))]
post.df = inner.join[,grep("Post", names(inner.join))]

df = cbind(HC.df, pre.df, post.df)
row.names(df) = row.names(HC.df)
#########################################################
#
# !!!!!!!!!!!!!!!!!!REMOVE NA'S HERE!!!!!!!!!!!!!!!!!!!!!
#
#########################################################
df = as.data.frame(lapply(df, as.numeric))
row.names(df) = row.names(HC.df)
for(i in c(1:ncol(df))){
  df[is.na(df[,i]),i ] = mean(df[,i], na.rm = TRUE)
}


# remove top subjets with largest metabolite values
sapply(df, max) -> max.stuff
names(max.stuff[order(max.stuff, decreasing = TRUE)])[c(1:10)] -> top.HC

df.clean = df[,-which(names(df) %in% top.HC)]
row.names(df.clean)

sample.type = colnames(df.clean)
sample.type = replace(sample.type, grep("Meta", sample.type), "black")
sample.type = replace(sample.type, grep("Pre", sample.type), "magenta")
sample.type = replace(sample.type, grep("Post", sample.type), "deepskyblue")


pdf("z-scale-heatmap.pdf")
heatmap.2(as.matrix(df.clean), trace = "none", ColSideColors = sample.type, 
          margins = c(8,20), cexRow = .1, cexCol = .4)
legend("topright",      
       legend = c("HC", "Pre", "Post"),
       col = c("black", "magenta", "deepskyblue"), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)
dev.off()


df.unknown.removed = df.clean[-grep("X - ", row.names(df.clean)), ]


pdf("z-scale-heatmap-unknowns-removed.pdf")
heatmap.2(as.matrix(df.unknown.removed), trace = "none", ColSideColors = sample.type, 
          margins = c(8,20), cexRow = .1, cexCol = .4)
legend("topright",      
       legend = c("HC", "Pre", "Post"),
       col = c("black", "magenta", "deepskyblue"), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)
dev.off()

pdf("z-scale-heatmap-unknowns-removed-log-scaled.pdf")
heatmap.2(as.matrix(log2(11+df.unknown.removed)), trace = "none", ColSideColors = sample.type, 
          margins = c(8,20), cexRow = .1, cexCol = .4)
legend("topright",      
       legend = c("HC", "Pre", "Post"),
       col = c("black", "magenta", "deepskyblue"), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)
dev.off()

# ## subjects qq-plots
# for(i in c(1:ncol(df))){
#   pdf(paste(names(df)[i], "qqplot.pdf"))
#   qqnorm(df[,i], main = names(df)[i]); qqline(df[,i])
#   dev.off()
# }
# 
# ## chemicals qq-plots
# for(i in c(1:nrow(df))){
#   pdf(paste("./figures/z-scores-qq-plots/chemicals/", make.names(row.names(df)[i]), "qqplot.pdf"))
#   qqnorm(df[i,], main = row.names(df)[i]); qqline(df[i,])
#   dev.off()
# }
