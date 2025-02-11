library(officer)
library(flextable)

install.packages("flextable",dependencies = T)
install.packages("magick",dependencies = T)
library(flextable)

word = read_docx()
  body_end_section_continuous(x = word)
  setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot")
  file = list.files()
  file = file[grepl(x = file, pattern = "AG and TC Variants Only",fixed = T)]
  drugnames = substr(x = file, start = 1, stop = (unlist(gregexpr(text = file, pattern = "AG and TC Variants Only",fixed = T))-2))
  drugnames = gsub(":", " ", drugnames)
  drugnames = substr(drugnames,1,stop = ifelse(test = nchar(drugnames)>31,yes = 29,no = nchar(drugnames)))
  drugnames = vapply(X =  drugnames, function(X) substr(X,1,stop = ifelse(test = nchar(drugnames)>31,yes = 31,no = nchar(drugnames))),FUN.VALUE = c("a"))
  drugnames=make.unique(tolower(drugnames),sep = ".")
  for(i in 1:5){
    f=read.csv(file[i], header = T, row.names = 1,stringsAsFactors = F, check.names = F)
    aa=flextable (f)
    aa = autofit(aa)
    body_add_flextable(x = word,split = F,value = aa)
  }
  body_end_section_landscape(word)
print(word, target = "test.docx")

