setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE")
ccle_ADAR = read.delim("/media/ahadli/Data/Farid/AOG_lab/CCLE-rna-seq/CCLE_RNAseq_genes_rpkm_20180929.gct", sep = "\t")
head(ccle_ADAR$Name)
sum(ccle_ADAR$Description == "ADAR")
which(ccle_ADAR$Description%in%"ADAR")
ccle_ADAR = ccle_ADAR[3238,]
rownames(ccle_ADAR) = ccle_ADAR$Description
ccle_ADAR=ccle_ADAR[-1]
head(ccle_ADAR)
ccle_intestine =ccle_ADAR[grepl("large_intestine", colnames(ccle_ADAR), ignore.case = T)]
colnames(ccle_intestine)
regexpr("_",colnames(ccle_intestine),fixed = F)
colnames(ccle_intestine) = substr(colnames(ccle_intestine), 1, (regexpr("_",colnames(ccle_intestine),fixed = F)-1))
colnames(ccle_intestine)
ccle_intestine_ADAR_log2 = log2(ccle_intestine)
depth_ccle = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/bam_files_depth.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
rownames(depth_ccle)
cellLineNames = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/cell_line_names.csv", header = T,row.names = 1, check.names = F, stringsAsFactors = F)
cellLineNames$x = gsub("[^[:alnum:]]","",cellLineNames$x)
cellLineNames$x = toupper(cellLineNames$x)
cellLineNames$x[55] = "SNUC1"
cellLineNames$x

mtch=match(cellLineNames$x,colnames(ccle_intestine_ADAR_log2))
mtch = mtch[!is.na(mtch)]
selected_ccle_intestine_ADAR_log2 = ccle_intestine_ADAR_log2[mtch]
rownames(depth_ccle) = cellLineNames$x
selected_depth_ccle = depth_ccle[mtch,]
edit_ratio = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/EditingRatio_all_cell_lines.csv", header = T,row.names = 1, check.names = F, stringsAsFactors = F)
edit_ratio = edit_ratio[mtch]
edit_ratio_bool = edit_ratio > 0
edit_AG_number=colSums(edit_ratio_bool)
names(edit_AG_number) = colnames(selected_ccle_intestine_ADAR_log2)
selected_depth_ccle=as.numeric(gsub("Average = ","",selected_depth_ccle))
names(selected_depth_ccle) = colnames(selected_ccle_intestine_ADAR_log2)
selected_depth_ccle = selected_depth_ccle[order(selected_depth_ccle)]

r_values = c()
p_values = c()
pdf(file = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/ADAR vs Editing_perdepth/AG_snv/cor_plot.pdf" )
for(i in 1:52){
  l_depth = selected_depth_ccle[i:54]
  l_ADAR = selected_ccle_intestine_ADAR_log2[names(l_depth)]
  l_ADAR = l_ADAR[order(as.numeric(l_ADAR[1,]))]
  l_AG_edit = edit_AG_number[names(l_ADAR)]
  r_values = c(r_values,cor.test(as.numeric(l_ADAR[1,]), l_AG_edit, method = "pearson")$estimate)
  p_values = c(p_values,cor.test(as.numeric(l_ADAR[1,]), l_AG_edit,method =  "pearson")$p.value)
  scatter.smooth(as.numeric(l_ADAR[1,]), l_AG_edit,span = 100,main = paste0("bigger than or equal to ",l_depth[1]," Pearson R = ",r_values[length(r_values)]),xlab = "ADAR expression (log2(RPKM))",ylab = "SNV_AG")
}
dev.off()
scatter.smooth(selected_depth_ccle[1:52], r_values,degree = 2, xlab = "Depth Cutoff (bigger than or equal to)", ylab = "Pearson R")


########### with output of giremi
edit_ratio_giremi = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/EditingRatio_all_cell_lines.csv", header = T,row.names = 1, check.names = F, stringsAsFactors = F)
edit_ratio_giremi = edit_ratio_giremi[mtch]
edit_ratio_giremi = edit_ratio_giremi[grepl("|AG|",rownames(edit_ratio_giremi),fixed = T),]
edit_ratio_giremi_bool = edit_ratio_giremi > 0
edit_AG_number_giremi=colSums(edit_ratio_giremi_bool)
depth_ccle = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/bam_files_depth.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
rownames(depth_ccle) = cellLineNames$x
selected_depth_ccle = depth_ccle[mtch,]
names(edit_AG_number_giremi) = colnames(selected_ccle_intestine_ADAR_log2)
selected_depth_ccle=as.numeric(gsub("Average = ","",selected_depth_ccle))
names(selected_depth_ccle) = colnames(selected_ccle_intestine_ADAR_log2)
selected_depth_ccle = selected_depth_ccle[order(selected_depth_ccle)]

r_values_giremi = c()
p_values_giremi = c()
pdf(file = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/ADAR vs Editing_perdepth/giremi_AG_only/cor_plot.pdf" )
for(i in 1:52){
  l_depth = selected_depth_ccle[i:54]
  l_ADAR = selected_ccle_intestine_ADAR_log2[names(l_depth)]
  l_ADAR = l_ADAR[order(as.numeric(l_ADAR[1,]))]
  l_AG_edit = edit_AG_number_giremi[names(l_ADAR)]
  r_values_giremi = c(r_values_giremi,cor.test(as.numeric(l_ADAR[1,]), l_AG_edit, method = "pearson")$estimate)
  p_values_giremi = c(p_values_giremi,cor.test(as.numeric(l_ADAR[1,]), l_AG_edit,method =  "pearson")$p.value)
  scatter.smooth(as.numeric(l_ADAR[1,]), l_AG_edit,span = 100,main = paste0("bigger than or equal to ",l_depth[1]," Pearson R = ",r_values_giremi[length(r_values_giremi)]),xlab = "ADAR expression (log2(RPKM))",ylab = "GIREMI_AG")
}
dev.off()
scatter.smooth(selected_depth_ccle[1:52], r_values_giremi,degree = 2,xlab = "Depth Cutoff (bigger than or equal to)", ylab = "Pearson R")


bll = as.matrix(r_values_giremi)
bll=t(bll)
