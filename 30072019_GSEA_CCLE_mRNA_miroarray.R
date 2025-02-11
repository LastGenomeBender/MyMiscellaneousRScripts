##### read and edit the mRNA array
ccle_mrna_array=read.delim("/home/ahadli/Desktop/CCLE_Expression_Entrez_2012-09-29.gct",fill = T,stringsAsFactors = F,sep = "\t",header = F)
ccle_mrna_array = ccle_mrna_array[c(-1,-2),]
arr = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
name_as_in_ccle = as.vector(as.character(arr[2,c(3:ncol(arr))]))
ccle_mrna_selected = ccle_mrna_array[c(1,2,which(ccle_mrna_array[1,]%in%name_as_in_ccle))]
meta = arr[c(1,2)]
arr = arr[c(-1,-2),c(-1,-2)]
arr = as.data.frame(sapply(arr, as.numeric))
arr = log2(arr+1)
names(name_as_in_ccle) = colnames(arr)
edit = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
edit_selected = edit[,(colnames(arr)[1:ncol(arr)])]
AHR_edit_rppa_selected = edit_selected[grepl(pattern ="7|17384354",x =  rownames(edit_selected),fixed = T),]
AHR_nonedit = colnames(AHR_edit_rppa_selected[,AHR_edit_rppa_selected==0])
AHR_edit = setdiff(colnames(AHR_edit_rppa_selected),y = AHR_nonedit)
whcih
colnames(ccle_mrna_selected) = ccle_mrna_selected[1,]
ccle_mrna_selected_ordered = ccle_mrna_selected[c("Name","Description",name_as_in_ccle[AHR_edit], name_as_in_ccle[AHR_nonedit])]
write.table(ccle_mrna_selected_ordered, "ccle_mrna_AHRedit_vs_AHRnonedit.txt",sep = "\t")
  getwd(
)
