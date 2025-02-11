# common between febbits and its comparison to qPCR
common_febbits =  read.csv("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project/18_common_febbits.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
rownames(common_febbits)
g61741 =  read.csv("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project/GSE61741/significant_ttest_all_1.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
g24709 = read.csv("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project/GSE24709/significant_ttest_all_1.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
qPCR = read.csv("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project/GSE64591/GSE64591_normalized.csv", header = T, stringsAsFactors = F, check.names = F)

#### annotate the qpcr dataset
matrix = qPCR
matrix$newAnnot = substr(matrix$ID_REF, 1 , ifelse(test = regexpr("#", matrix$ID_REF)>0, yes = (regexpr("#", matrix$ID_REF)-1), no = (nchar(matrix$ID_REF)-7))) 
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
qPCR = matrix
#### 18 common between 61741
temp = data.frame(row.names = g61741$newAnnot)
rownames(g61741) = g61741$newAnnot
common_g61741 = g61741[rownames(common_febbits),]
setwd("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project")
write.csv(common_g61741, "commons_between_febbits_gse61741.csv")

###### 18 common between 24709
temp = data.frame(row.names = g24709$newAnnot)
rownames(g24709) = g24709$newAnnot
common_g24709 = g24709[rownames(common_febbits),]
setwd("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project")
write.csv(common_g24709, "commons_between_febbits_gse24709.csv")
############ commons between qpcr
c_qpcr=data.frame()
for( i in 1:nrow(common_febbits)){
  bool = rownames(common_febbits)[i] == qPCR$newAnnot
  dt = qPCR[bool,]
  c_qpcr = rbind(c_qpcr,dt)
}

write.csv(c_qpcr, "commons_between_febbits_gse64591_qPCR.csv")
