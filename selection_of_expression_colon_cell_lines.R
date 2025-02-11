expr_ccle = read.table( "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_RNAseq_genes_rpkm_20180929.gct", check.names = F, stringsAsFactors = F,sep = "\t",fill=T)
expnms2 = as.character(expr_ccle[2,])
expnms2_b1 = grepl(x = expnms2, pattern = "_FIBROBLAST", fixed = T) 
expnms2 = expnms2[expnms2_b1]
expnms=read.table("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/cell_line_names_expression.csv", header = T, row.names = 1, check.names = F, stringsAsFactors = F,sep = "\t")
expnms = as.data.frame(t(expnms))
edit = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
editnms = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/cell_line_names.csv", header = T, stringsAsFactors = F, row.names = 1, check.names = F)
editnms = editnms[colnames(edit),]
class(gregexpr(text="aaadaad", "d", fixed = T)[[1]])
vec_expnms = c(rownames(expnms))
underline_indices = regexpr(text = vec_expnms, pattern = "_", fixed = T)
raw_expnms = substr( rownames(expnms), 1,(underline_indices -1))
raw_expnms=gsub("[^[:alnum:] ]", "", raw_expnms)
raw_expnms = toupper(raw_expnms)
raw_editnms = gsub("[^[:alnum:] ]", "", editnms)
raw_editnms = toupper(raw_editnms)

#rownames(expnms) = raw_expnms
commons=intersect(raw_expnms, raw_editnms)
indices = match(x = commons, table = raw_expnms) + 2
indices2 = match(x = commons, table = raw_editnms)

my_colon_expr = expr_ccle[,c(1,2,indices)]
add_edit_nms = editnms[indices2]
my_colon_expr = my_colon_expr[-1,]
colnms = colnames(edit)[indices2]
head(my_colon_expr[1,])
getwd()
colnames(my_colon_expr)[3:ncol(my_colon_expr)] = colnms
my_colon_expr = rbind(c(".",".",add_edit_nms),my_colon_expr )
write.csv(my_colon_expr, "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv")

