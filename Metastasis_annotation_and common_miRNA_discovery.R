
##################### for liver annotation to miRBAse 21
setwd("/home/ahadli/Desktop/Secil/")
source("https://bioconductor.org/biocLite.R")
biocLite("miRNAmeConverter")
matrix = read.csv("/home/ahadli/Desktop/Secil/liver_metas_sig_forannot.csv", row.names=1, header = T)
pos = regexpr(':', rownames(matrix))
matrix$newAnnot = substr(rownames(matrix),1, pos-1)
library(miRNAmeConverter)
nc = MiRNANameConverter() # Create MiRNANameConverter object
for(i in 1:nrow(matrix)){
  if(paste0(translateMiRNAName(nc,matrix$newAnnot[i])$v21.0, "a") == "a"){
    matrix$newAnnot[i] = "a"
  }
  else{
    matrix$newAnnot[i] = translateMiRNAName(nc,matrix$newAnnot[i])$v21.0
  }
}
write.csv(matrix,"liver_metas_sig_forannot.csv")

######################## for lymph
setwd("/home/ahadli/Desktop/Secil/")
source("https://bioconductor.org/biocLite.R")
biocLite("miRNAmeConverter")
matrix = read.csv("lymph_node_met_sig_forannot.csv", row.names=1, header = T)
pos = regexpr(':', rownames(matrix))
matrix$newAnnot = substr(rownames(matrix),1, pos-1)
library(miRNAmeConverter)
nc = MiRNANameConverter() # Create MiRNANameConverter object
for(i in 1:nrow(matrix)){
  if(paste0(translateMiRNAName(nc,matrix$newAnnot[i])$v21.0, "a") == "a"){
    matrix$newAnnot[i] = "a"
  }
  else{
    matrix$newAnnot[i] = translateMiRNAName(nc,matrix$newAnnot[i])$v21.0
  }
}
write.csv(matrix,"lymph_node_met_sig_forannot.csv")

#####################3 finding pre and postop similar between metastasis and pre-postop
preop= read.csv("preop-postop_samplewise.csv", header = T, row.names = 1, stringsAsFactors = F)
liv_met = read.csv("liver_metas_sig_forannot.csv", header = T, row.names = 1, stringsAsFactors = F)
lym_met = read.csv("lymph_node_met_sig_forannot.csv", header = T, row.names = 1)

common_bet_liver = c()
for(i in 1:nrow(liv_met)){
  for (j in 1:nrow(preop)){
    if((preop$new.annotation[j]!="a")|(liv_met$newAnnot[i]!="a")){
      if(liv_met$newAnnot[i]==preop$new.annotation[j]){
        common_bet_liver = c(common_bet_liver, preop$new.annotation[j])
      }
    }
  }
}
pos_da=which(liv_met$newAnnot %in% common_bet_liver)
common_bet_liver_mat = liv_met[pos_da,]
#######################################Node
common_bet_lym = c()
for(i in 1:nrow(lym_met)){
  for (j in 1:nrow(preop)){
    if((preop$new.annotation[j]!="a")|(lym_met$newAnnot[i]!="a")){
      if(lym_met$newAnnot[i]==preop$new.annotation[j]){
        common_bet_lym = c(common_bet_lym, preop$new.annotation[j])
      }
    }
  }
}

pos_da2=which(lym_met$newAnnot %in% common_bet_lym)
common_bet_lym_mat = lym_met[pos_da2,]
