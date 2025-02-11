#edit_mat = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
cell_lines_nms = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/cell_line_names.csv",row.names = 1, header = T, stringsAsFactors = F, check.names = F)
#rownames(cell_lines_nms) = colnames(edit_mat)
#write.csv(cell_lines_nms,"/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/cell_line_names.csv")
rownames(cell_lines_nms)
cell_lines_nms = cell_lines_nms$x
depth = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/bam_files_depth.csv",row.names = 1, header = T, stringsAsFactors = F, check.names = F)
depth = depth$x 
depth = substr(depth, 12, nchar(depth))
options(digits = 15)
depth = as.numeric(depth)
print(depth)
names(depth) = cell_lines_nms
edit_mat_TC = edit_mat[grepl("\\|TC\\|" , rownames(edit_mat))|grepl("\\|AG\\|", rownames(edit_mat)) ,] #"\\|AG\\|" for AG
edit_mat_bool_TC = edit_mat_TC>0
edit_mat_bool2_TC = colSums(edit_mat_bool_TC)
edit_mat_bool2_TC
setwd("/home/ahadli/Desktop")
write.csv(edit_mat_bool2_TC, "AG and TC Editings.csv")
write.csv(depth, "CCLE colon cancer RNA-seq files depths.csv")
names(edit_mat_bool2_TC) = cell_lines_nms
axis(1, xaxp=c(2,1, 19), las=2)
scatter.smooth(depth,edit_mat_bool2_TC,span = 10, xlab = "Coverage", ylab = "Editing Number", main ="AG and TC Editing")
depth=depth[order(depth)]
edit_mat_bool2_TC = edit_mat_bool2_TC[names(depth)]
R_vec = c()
p_vec = c()
for(i in 1:(length(depth)-1)){
  temp_depth = depth[i+1:length(depth)]
  temp_edit = edit_mat_bool2_TC[i+1:length(edit_mat_bool2_TC)]
  R=cor.test(temp_depth,temp_edit, method = "pearson")$estimate
  p=cor.test(temp_depth,temp_edit,method = "pearson")$p.value
  p_vec=c(p_vec,p)
  R_vec = c(R_vec, R)
}

plot(c(1:length(p_vec)), p_vec)
# edit_mat_AG = edit_mat[grepl("\\|AG\\|", rownames(edit_mat)),] #"\\|AG\\|" for AG
# edit_mat_bool_AG = edit_mat_AG>0
# edit_mat_bool2_AG = colSums(edit_mat_bool_AG)
# edit_mat_bool2_AG
# names(edit_mat_bool2_AG) = cell_lines_nms
# scatter.smooth(depth,edit_mat_bool2_AG,span = 10, xlab = "Coverage", ylab = "SNV(AG) Number", main ="SNV (AG)")
