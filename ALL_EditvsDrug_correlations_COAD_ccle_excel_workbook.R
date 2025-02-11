## all correlation data into a file.
#install.packages("openxlsx")
library(openxlsx)
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot")
file = list.files()
file = file[grepl(x = file, pattern = "AG and TC Variants Only",fixed = T)]
drugnames = substr(x = file, start = 1, stop = (unlist(gregexpr(text = file, pattern = "AG and TC Variants Only",fixed = T))-2))
drugnames = gsub(":", " ", drugnames)
drugnames = substr(drugnames,1,stop = ifelse(test = nchar(drugnames)>31,yes = 29,no = nchar(drugnames)))
drugnames = vapply(X =  drugnames, function(X) substr(X,1,stop = ifelse(test = nchar(drugnames)>31,yes = 31,no = nchar(drugnames))),FUN.VALUE = c("a"))
wb <- createWorkbook()
drugnames=make.unique(tolower(drugnames),sep = ".")
for(i in 1:length(file)){
  f=read.csv(file[i], header = T, row.names = 1,stringsAsFactors = F, check.names = F)
  addWorksheet(wb,sheetName = drugnames[i])
  writeData(x = f,wb = wb, sheet = drugnames[i],colNames = T, rowNames = T)
}
saveWorkbook(wb, file = "090719_ALL_EditvsDrug_correlations_COAD_ccle.xlsx", overwrite = TRUE)
for(i in 1:length(file)){
drugnames[23]
is.character('asasas')
class(drugnames)
drugnames[221]
aa=drugnames[drugnames=="erlotinib PLX-4032 (2 1 mol mol"]
