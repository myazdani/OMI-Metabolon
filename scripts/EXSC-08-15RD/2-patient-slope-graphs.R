###
# slope-graphs.R
###

setwd("~/Documents/OMI/Metabolon/")
library(readxl)
df = read.csv("./rawData/EXSC-08-15RD/clean-CSV/EXSC-08-15RD-Xenobiotic-and-Unknown-REMOVED.csv", header = TRUE, stringsAsFactors = FALSE)

patients.df = df[,grep("REFERENCE.VS", names(df))]
row.names(patients.df) = df$Biochemical.Name

#########################################################
#
#  TRANSPOSE DATA
#
#########################################################
# first remember the names
n <- row.names(patients.df)
# transpose all but the first column (name)
biochem.df <- as.data.frame(t(patients.df))
colnames(biochem.df) <- n
names(biochem.df) = make.names(names(biochem.df))

biochem.df$sample = sapply(row.names(biochem.df), FUN = function(x) strsplit(x, split = "REFERENCE.VS.")[[1]][2])

#########################################################
#
# MERGE WITH META DATA
#
#########################################################
df.keys = as.data.frame(read_excel("./rawData/EXSC-08-15RD/20150722_Metabolon_CFS_HC_Study_001.xlsx", sheet = 1, col_names = TRUE))
names(df.keys) = make.names(names(df.keys))
df.keys$Tube.ID.. = toupper(df.keys$Tube.ID..)

df.cfs = merge(df.keys, biochem.df, by.x = "Tube.ID..", by.y = "sample", all = FALSE)

#########################################################
#
# find largest differences
#
#########################################################
library(reshape)
df.m = melt(df.cfs, id.vars = names(df.cfs)[c(1:9)])


patient.df = unique(df.m[,c("Cohorts", "Patient.ID..")])
print(as.numeric(names(which(table(patient.df$Patient.ID..) < 2))))
## remove subject 9 and 730 since they do not have values in both pre and post
df.m = subset(df.m, Patient.ID.. != 9 & Patient.ID.. != 730 & Patient.ID.. != 78 & Patient.ID.. != 311)


library(dplyr)
df.m %>%
  group_by(Patient.ID.., variable) %>%
  summarise(abs.diff = abs(value[1] - value[2])) %>%
  as.data.frame() -> df.diff

p = ggplot(df.diff, aes(x = variable, y = abs.diff, label = variable)) + geom_text() + facet_wrap(~Patient.ID..) +  theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                                     axis.title.x=element_blank(),
                                                                                                     axis.title.y=element_blank(),legend.position="none",
                                                                                                     panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                                                     panel.grid.minor=element_blank(),plot.background=element_blank())

ggplotly(p, filename = "/metabolon/pre-post-metabolyte-diff", fileopt = "overwrite")

top.10.selector = function(df){
  df.sorted = df[order(df$abs.diff, decreasing = TRUE)[c(1:20)],]
  return(df.sorted)
}


top10.diff = ddply(df.diff, .(Patient.ID..), top.10.selector)

p = ggplot(top10.diff, aes(x = variable, y = abs.diff, label = variable)) + geom_text() + facet_wrap(~Patient.ID..) +  theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                                                          axis.title.x=element_blank(),
                                                                                                                          axis.title.y=element_blank(),legend.position="none",
                                                                                                                          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                                                                          panel.grid.minor=element_blank(),plot.background=element_blank())
ggplotly(p, filename = "/metabolon/pre-post-metabolyte-top10-diff", fileopt = "overwrite")

top.10.lists = function(df){
  return(df$variable)
}

patient.lists = dlply(top10.diff, .(Patient.ID..), top.10.lists )
df.m$Cohorts = factor(df.m$Cohorts, levels = c("CFS Pre", "CFS Post"))
for (i in c(1:length(names(patient.lists)))){
  df.m.patient = subset(df.m, Patient.ID.. == as.numeric(names(patient.lists)[i]) & variable %in% patient.lists[[i]] )
  df.m.patient$variable = as.factor(df.m.patient$variable)
  ggplot(df.m.patient, aes(x = Cohorts, y = value, colour = variable, group= variable)) +
    geom_line() + geom_point(shape=21, fill = "white") + xlab("")  + 
    ggtitle(paste("Patient", names(patient.lists)[i]))-> p
  png(filename = paste0("figures/patient-top-20-M",names(patient.lists)[i], ".png"))
  print(p)
  dev.off()
}
