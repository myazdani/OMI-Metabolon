###
# confidence-intervales.R
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

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


dfwc_between <- summarySE(data=df.cfs, measurevar="cystine", groupvars="Cohorts", na.rm=FALSE, conf.interval=.95)


ggplot(dfwc_between, aes(x=Cohorts, y=cystine, group=1)) +
  geom_line() +
  geom_errorbar(width=.1, aes(ymin=cystine-ci, ymax=cystine+ci), colour="red") +
  geom_errorbar(width=.1, aes(ymin=cystine-ci, ymax=cystine+ci), data=dfwc_between) +
  geom_point(shape=21, size=3, fill="white")

df.cfs$Patient.ID.. = as.factor(df.cfs$Patient.ID..)
ggplot(df.cfs, aes(x = Cohorts, y = cystine, colour = Patient.ID.., group= Patient.ID..)) +
  geom_line() + geom_point(shape=21, fill = "white")

## outlier subjects are Patient.ID.. == 730 and Patient.ID.. == 9