#edit_mat = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
edit_mat = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
edit_mat2 = edit_mat
for (i in 1:(ncol(edit_mat))){
  mat = (as.vector(as.numeric(edit_mat[,i])))
  mat[mat==0] <- NA
  edit_mat[i] = mat
}
load("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/all_lists.RData")
cell_lines = names(all_files)
dots = gregexpr("\\.", cell_lines )


cell_lines_nms = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/cell_line_names.csv",row.names = 1, header = T)
cell_lines_nms = cell_lines_nms[colnames(edit_mat),]
cell_lines_nms = as.character(cell_lines_nms)

#cell_lines_nms = c()
# for(i in 1:length(dots)){
#   dot_pos = dots[[i]]
#   new_nm = substr(cell_lines[i], dot_pos[1]+1, dot_pos[2] -1)
#   cell_lines_nms = c(cell_lines_nms, new_nm)
# }
# cell_lines_nms
# cell_lines[55]
# cell_lines_nms[55]= "SNU-C15"
# write.csv(cell_lines_nms, "cell_line_names")
ccle=read.delim("/media/ahadli/Data/Farid/AOG_lab/transepigenomics/CCLE Published.txt",  header = T, stringsAsFactors = F, check.names = F)


c_cnames_ccle = ccle$`Primary Cell Line Name`
c_cnames_ccle= gsub("[^[:alnum:] ]", "", c_cnames_ccle)
c_cnames_ccle = gsub(" ", "", c_cnames_ccle)
c_cnames_ccle
cell_lines_nms_removed = gsub("[^[:alnum:] ]", "", cell_lines_nms)

cell_lines_nms_removed=gsub(" ", "",cell_lines_nms_removed )
cell_lines_nms_removed
head(ccle)
diff_compounds = levels(factor(ccle$Compound))
compoundwise_ccle = list()
for (i in diff_compounds){
  ccle_s_compound = ccle[ccle$Compound==i,]
  compoundwise_ccle[[i]] = ccle_s_compound
}

for(i in diff_compounds){
  mat = compoundwise_ccle[[i]]
  cnames_ccle = mat$`Primary Cell Line Name`
  cnames_ccle= gsub("[^[:alnum:] ]", "", cnames_ccle)
  cnames_ccle = gsub(" ", "", cnames_ccle)
  bool = c()
  for(j in cnames_ccle){
    tmp_bool = ifelse(test = (sum(toupper(j) == toupper(cell_lines_nms_removed))) > 0, yes = TRUE, no = FALSE)
    bool = c(bool,tmp_bool )
  }
  mat2 = mat[bool,]
  mat2=mat2[order(mat2$ActArea),]
  compoundwise_ccle[[i]] = mat2
}
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/A-to-I GIREMI Colon CCLE")
save(compoundwise_ccle, file = "compoundwise_ccle.RData")
correlations_compound = list()
c1 = 0
for(i in diff_compounds){
  c1 = c1+1
  cat("compound number",c1 )
  mat = compoundwise_ccle[[i]][c(1,8)]
  rnm_ccle = rownames(mat)
  rnm_ccle = mat$`Primary Cell Line Name`
  rnm_ccle= gsub("[^[:alnum:] ]", "", rnm_ccle)
  rnm_ccle = gsub(" ", "", rnm_ccle)
  bool = c()
  for(j in rnm_ccle){
    pos= which(toupper(cell_lines_nms_removed) == toupper(j))
    bool = c(bool,pos)
  }
  edit_temp = edit_mat[bool]
  temp_cell_names = cell_lines_nms[bool]
  temp_cell_names = c(temp_cell_names,c("O","O","O","O"))
  #edit_temp = edit_temp[apply(edit_temp, 1, FUN=function(x) sum((sum(!(x==0)))))>5,]
  edit_temp = edit_temp[apply(edit_temp, 1, FUN=function(x) sum((sum(!is.na(x)))))>5,] 
  sp_pval = c()
  pe_pval = c()
  sp_r = c()
  pe_r =c()
  for(j in 1:nrow(edit_temp)){
    # rp =  cor.test(as.vector(as.numeric(edit_temp[j,])), mat$ActArea, method = "pearson",)$estimate
    # pp = cor.test(as.vector(as.numeric(edit_temp[j,])), mat$ActArea, method = "pearson")$p.value
    # rs = cor.test(as.vector(as.numeric(edit_temp[j,])), mat$ActArea, method = "spearman")$estimate
    # sp = cor.test(as.vector(as.numeric(edit_temp[j,])), mat$ActArea, method = "spearman")$p.value
    rp =  cor.test(as.vector(as.numeric(edit_temp[j,])), mat$ActArea, method = "pearson", use = "na.or.complete")$estimate
    pp = cor.test(as.vector(as.numeric(edit_temp[j,])), mat$ActArea, method = "pearson", use = "na.or.complete")$p.value
    rs = cor.test(as.vector(as.numeric(edit_temp[j,])), mat$ActArea, method = "spearman", use = "na.or.complete")$estimate
    sp = cor.test(as.vector(as.numeric(edit_temp[j,])), mat$ActArea, method = "spearman", use = "na.or.complete")$p.value
    sp_pval = c(sp_pval, sp)
    pe_pval = c(pe_pval, pp)
    sp_r = c(sp_r,rs)
    pe_r =c(pe_r, rp)
  }
  stat = as.data.frame(cbind(pe_r, pe_pval,sp_r,sp_pval))
  rownames(stat) = rownames(edit_temp)
  colnames(stat) = c("Pearson R", "Pearson Pval", "Spearman R", "Spearman Pval")
  stat = stat[order(stat$`Pearson R`,decreasing = F),]
  edit_temp = edit_temp[rownames(stat),]
  stat = cbind(edit_temp, stat)
  activity_area = c(mat$ActArea,0,0,0,0)
  addition = rbind(temp_cell_names,activity_area)
  rownames(addition) = c("Cell Line", "Activity Area")
  colnames(addition) = colnames(stat)
  stat = rbind(addition,stat,stringsAsFactors =F )
  correlations_compound[[i]] =stat
}
#dir.create("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations")
#setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations")
dir.create("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted")
 setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted")

lapply(correlations_compound, function(x) max(x[[c(1:length(x))]]$`Pearson R`) )

for(i in 1:length(correlations_compound)){
  mat = correlations_compound[[i]]
  write.csv(mat, paste0(names(correlations_compound)[i],"_correlation.csv"))
}
getwd()
my_int=correlations_compound[["Topotecan"]][3,]
my_int = my_int[-((ncol(my_int)-3):ncol(my_int))]
my_int2 = correlations_compound[["Topotecan"]][2,]
my_int2 = my_int2[-((ncol(my_int2)-3):ncol(my_int2))]
scatter.smooth(my_int, my_int2, xlab = "Editing Level", ylab = "Activity Area (Topotecan)", main = "RPL7L1(6|42856029|+|AG|)")
nope=read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/Compouds_correlations/Paclitaxel_correlation.csv", row.names = 1, header = T, stringsAsFactors = F)
head(nope)
colnames(aaaa)
scatter.smooth(nope[3,1:(ncol(nope)-4)],nope[2,1:(ncol(nope)-4)],span = 10000000,degree = 1, xlab = "Editing Level",ylab = "Activity Area (Paclitaxel)", main = "RPL7L1                R = -0.74")
text(labels = "R = -0.74",x = 1,y = 1,pos = c(3,4))
aaaa = cbind(nope[3,1:(ncol(nope)-4)],nope[2,1:(ncol(nope)-4)],stringsAsFactors=F)
aaaaaaa = lm(nope[3,1:(ncol(nope)-4)]~nope[2,1:(ncol(nope)-4)])
scatter
install.packages("car")
library(car)
