###
# clean-data.R
# take the sheet from excel (manually exported as a CSV) and clean up the header row 
# and replace blanks with NA
# save the cleaned up data as clean-metabolite-counts.csv
###
setwd("~/Documents/OMI/Metabolon/")
library(readxl)

raw.df = as.data.frame(read_excel("./rawData/ECSC-08-15-from-Adam/metabolite-counts-sheet-from-excel.xlsx", col_names = TRUE))

##
## create clean header row
##
names(raw.df)[which(names(raw.df) == "PATIENT_ID"):ncol(raw.df)] = paste0("ID.", names(raw.df)[which(names(raw.df) == "PATIENT_ID"):ncol(raw.df)])
meta.header = make.names(names(raw.df))
real.header = raw.df[1, ]
names(raw.df) = paste0(names(raw.df), ".", real.header)
raw.df = raw.df[-1,]
##
## select relevant columns
##
clean.df = raw.df[, c(2, 13:ncol(raw.df))]
names(clean.df)[1] = "BIOCHEMICAL"
names(clean.df) = make.names(names(clean.df))
##
## ensure data is numeric type clean header row
## (replace balnks with NA)
##
clean.df[,-1] = lapply(clean.df[,-1], as.numeric)

write.table(clean.df, file = "./rawData/ECSC-08-15-from-Adam/clean-metabolite-counts.tsv", row.names = FALSE, quote = FALSE, sep = "\t")
