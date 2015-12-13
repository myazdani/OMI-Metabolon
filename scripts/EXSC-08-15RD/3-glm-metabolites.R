###
# glm-metabolites.R
###
setwd("~/Documents/OMI/Metabolon/")
library(readxl)
library(dplyr)
library(broom)
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
# melt data (put in long format)
#
#########################################################
library(reshape)
df.m = melt(df.cfs, id.vars = names(df.cfs)[c(1:9)])


#########################################################
#
# run glm per metabolites
#
#########################################################

## remove subject 9 and 730 since they do not have values in both pre and post
df.m = subset(df.m, Patient.ID.. != 9 & Patient.ID.. != 730 & Patient.ID.. != 78 & Patient.ID.. != 311)


df.m$Cohorts = as.factor(df.m$Cohorts)

## sample analysis for single metabolite
df.m %>%
  filter(variable %in% c("glycine") ) -> glycine
mod = glm(Cohorts ~ value,, family = "binomial", data = glycine)

## analysis for all metabolites
na.omit(df.m) %>%
  group_by(variable) %>%
  do(tidy(glm(Cohorts ~ value, family = "binomial", .))) %>%
  as.data.frame() -> linear_classifiers

na.omit(df.m) %>%
  group_by(variable) %>%
  do(glance(glm(Cohorts ~ value, family = "binomial", .))) %>%
  as.data.frame() -> linear_classifiers.stats


intercepts = subset(linear_classifiers, term == "(Intercept)")
slopes = subset(linear_classifiers, term == "value")
#top10.candidates = as.character(slopes$variable[order(slopes$p.value)[c(1:20)]])
top10.candidates = as.character(intercepts$variable[order(intercepts$p.value)[c(1:20)]])

top10.df = subset(df.m, variable %in% top10.candidates)
ggplot(top10.df, aes(x = variable, y = value, colour = Cohorts)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("") + ggtitle("Top intercept terms")


patient.list = unique(top10.df$Patient.ID..)
top10.df$variable = as.factor(top10.df$variable)
df.m$Cohorts = factor(df.m$Cohorts, levels = c("CFS Pre", "CFS Post"))
for (i in c(1:length(patient.list))){
  df.m.patient = subset(top10.df, Patient.ID.. == patient.list[i] )
  ggplot(df.m.patient, aes(x = Cohorts, y = value, colour = variable, group= variable)) +
    geom_line() + geom_point(shape=21, fill = "white") + xlab("")  + 
    ggtitle(paste("Patient", patient.list[i]))-> p
  png(filename = paste0("figures/patient-top-20-M",patient.list[i], ".png"))
  print(p)
  dev.off()
}