setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/a-to-i_all_func_annotated")
files = list.files(getwd())
all_files = list()
#dir.create(path = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/A-to-I GIREMI Colon CCLE"))
for ( i in 1:length(files)){
  #dir.create(path = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/filtered files/",files[i]))
  mat = read.csv(files[i], header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  mat = mat[mat$ifSNP==0,]
  edit_ratio = mat$major_ratio
  for(j in 1:nrow(mat)){
    if(mat$reference_base[j] == mat$major_base[j]){
      edit_ratio[j] = 1-edit_ratio[j]
    }
  }
  mat$editing_ratio=edit_ratio
  rownames(mat) = paste0(rownames(mat),"|", mat$RNAE_t, "|", mat$gene,"||||", mat$func_anno)
  #write.csv(mat,paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/filtered files/",files[i],"/",files[i]))
  all_files[[files[i]]]= mat
}
r = ls()[-1]
rm(c("r"))
rnms = c()
for ( i in 1:length(all_files)){
  rnms = c(rnms, rownames(all_files[[i]]))
}
bool = !duplicated(rnms)
rnms_filtered = rnms[bool]
edit_mat = data.frame(row.names = rnms_filtered)
for ( i in 1:length(all_files)){
  mat = all_files[[i]]
  mat = mat[rownames(edit_mat), c(24)]
  for(j in 1:length(mat)){
    if(is.na(mat[j])){
      mat[j] = 0
    }
  }
  edit_mat = cbind(edit_mat, mat)
}
colnames(edit_mat) = names(all_files)
sum(is.na(edit_mat))
write.csv(edit_mat, "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/EditingRatio_all_cell_lines.csv")

edit_mat2 = as.data.frame(t(edit_mat))
df_f <- edit_mat2[,apply(edit_mat2, 2, var, na.rm=TRUE) != 0]
edit_pca2 = prcomp(df_f )
summary(edit_pca2)
str(edit_pca2)
# install.packages("devtools")
# library(devtools)
# install_github("vqv/ggbiplot")
#library(ggbiplot)
ggbiplot(edit_pca2,choices = c(58,) ,var.axes = F)

install.packages("rgl")
install.packages("pca3d")
library(pca3d)
pca3d(edit_pca2 , components = c(1,2,3))
pca2d(edit_pca2 , components = c(1,3))
