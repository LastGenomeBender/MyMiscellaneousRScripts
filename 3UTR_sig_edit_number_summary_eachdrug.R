#### sigs summary file
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlation_significants_selected")
files = list.files()
files = files[grepl(x = files,pattern = "AG and TC Variants Only",fixed = T)]
path_to_cor_pvals = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot"
corr_files = list.files(path_to_cor_pvals)
corr_files = corr_files[grepl(x = corr_files,pattern = "AG and TC Variants Only",fixed = T)]
drg_names = substr(files,1,(regexpr(text = files,pattern = "AG and TC Variants Only",fixed = T)-2))
drg_names_cor = substr(corr_files,1,(regexpr(text = corr_files,pattern = "AG and TC Variants Only",fixed = T)-2))
idx = match(drg_names,drg_names_cor)
corr_files_selected = corr_files[idx]
setwd(path_to_cor_pvals)
sumary_stat = data.frame()
for(i in 1:length(corr_files_selected)){
  mat=read.csv(corr_files_selected[i], header = T, row.names = 1,stringsAsFactors = F,check.names = F)
  mat = mat[5:8]
  bool = mat<0.05
  final = colSums(bool)
  final[5] = sum(bool)
  sumary_stat = rbind(sumary_stat, final)
}
rownames(sumary_stat) = drg_names
colnames(sumary_stat) = c("Benjamini-Hochberg Correction Pearson","Bonferroni Correction Pearson","Benjamini-Hochberg Correction Spearman","Bonferroni Correction Spearman","Sum")
sumary_stat = sumary_stat[order(sumary_stat$Sum,decreasing = T),]
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlation_significants_selected")
write.csv(sumary_stat,"100719_Significants_Summary_Statistics_Each_Drug.csv")
