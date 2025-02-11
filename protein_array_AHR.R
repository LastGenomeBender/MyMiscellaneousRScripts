##### selection of my cells from RPPA array
arr = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
names_as_in_ccle = as.vector(as.character(arr[2, c(3:ncol(arr))]))
rppa = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_RPPA_20181003.csv",header = T, row.names = 1, stringsAsFactors = F, check.names = F)
rppa_selected = rppa[names_as_in_ccle,]
edit = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
edit_rppa_selected = edit[,(colnames(arr)[3:ncol(arr)])]
AHR_edit_rppa_selected = edit_rppa_selected[grepl(pattern ="7|17384354",x =  rownames(edit_rppa_selected),fixed = T),]
AHR_nonedit = colnames(AHR_edit_rppa_selected[,AHR_edit_rppa_selected==0])
AHR_edit = setdiff(colnames(AHR_edit_rppa_selected),y = AHR_nonedit)

library(BiocManager)
install("gage")
library(gage)
rppa_AHR_edit = arr[,AHR_edit]
rppa_AHR_edit = as.vector(as.character(rppa_AHR_edit[2,]))
rppa_AHR_edit = rppa_selected[rppa_AHR_edit,]

rppa_AHR_nonedit = arr[,AHR_nonedit]
rppa_AHR_nonedit = as.vector(as.character(rppa_AHR_nonedit[2,]))
rppa_AHR_nonedit = rppa_selected[rppa_AHR_nonedit,]



boxplot(rppa_AHR_edit$`beta-Catenin`,rppa_AHR_nonedit$`beta-Catenin`,outline = T)
stripchart(list(rppa_AHR_edit$`beta-Catenin`,rppa_AHR_nonedit$`beta-Catenin`), vertical = T, pch=20,method = "jitter",add = T,col ="red")
boxplot(rppa_AHR_edit$`beta-Catenin_pT41_S45`,rppa_AHR_nonedit$`beta-Catenin_pT41_S45`,outline = F)
stripchart(list(rppa_AHR_edit$`beta-Catenin_pT41_S45`,rppa_AHR_nonedit$`beta-Catenin_pT41_S45`), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$`beta-Catenin`, y = rppa_AHR_nonedit$`beta-Catenin`,paired = F)

boxplot(rppa_AHR_edit$Cyclin_D1,rppa_AHR_nonedit$Cyclin_D1,outline = T)
stripchart(list(rppa_AHR_edit$Cyclin_D1,rppa_AHR_nonedit$Cyclin_D1), vertical = T, pch=20,method = "jitter",add = T,col ="red")
boxplot(rppa_AHR_edit$Cyclin_B1,rppa_AHR_nonedit$Cyclin_B1,outline = T)
stripchart(list(rppa_AHR_edit$Cyclin_B1,rppa_AHR_nonedit$Cyclin_B1), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$Cyclin_B1, y = rppa_AHR_nonedit$Cyclin_B1,paired = F)


boxplot(rppa_AHR_edit$GSK3_pS9,rppa_AHR_nonedit$GSK3_pS9,outline = T)
stripchart(list(rppa_AHR_edit$GSK3_pS9,rppa_AHR_nonedit$GSK3_pS9), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$GSK3_pS9, y = rppa_AHR_nonedit$GSK3_pS9,paired = F)

boxplot(rppa_AHR_edit$`GSK3-alpha-beta_pS21_S9`,rppa_AHR_nonedit$`GSK3-alpha-beta_pS21_S9`,outline = T)
stripchart(list(rppa_AHR_edit$`GSK3-alpha-beta_pS21_S9`,rppa_AHR_nonedit$`GSK3-alpha-beta_pS21_S9`), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$`GSK3-alpha-beta_pS21_S9`, y = rppa_AHR_nonedit$`GSK3-alpha-beta_pS21_S9`,paired = F)

boxplot(rppa_AHR_edit$GSK3_pS9,rppa_AHR_nonedit$GSK3_pS9,outline = T)
stripchart(list(rppa_AHR_edit$GSK3_pS9,rppa_AHR_nonedit$GSK3_pS9), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$GSK3_pS9, y = rppa_AHR_nonedit$GSK3_pS9,paired = F)

boxplot(rppa_AHR_edit$Snail_Caution,rppa_AHR_nonedit$Snail_Caution,outline = T)
stripchart(list(rppa_AHR_edit$Snail_Caution,rppa_AHR_nonedit$Snail_Caution), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$Snail_Caution, y = rppa_AHR_nonedit$Snail_Caution,paired = F)

boxplot(rppa_AHR_edit$ADAR1,rppa_AHR_nonedit$ADAR1,outline = T)
stripchart(list(rppa_AHR_edit$Snail_Caution,rppa_AHR_nonedit$Snail_Caution), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$ADAR1, y = rppa_AHR_nonedit$ADAR1,paired = F)

boxplot(rppa_AHR_edit$`NF-kB-p65_pS536_Caution`,rppa_AHR_nonedit$`NF-kB-p65_pS536_Caution`,outline = T)
stripchart(list(rppa_AHR_edit$`NF-kB-p65_pS536_Caution`,rppa_AHR_nonedit$`NF-kB-p65_pS536_Caution`), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$`NF-kB-p65_pS536_Caution`, y = rppa_AHR_nonedit$`NF-kB-p65_pS536_Caution`,paired = F)


boxplot(rppa_AHR_edit$`c-Met_pY1235`,rppa_AHR_nonedit$`c-Met_pY1235`,outline = T)
stripchart(list(rppa_AHR_edit$`NF-kB-p65_pS536_Caution`,rppa_AHR_nonedit$`NF-kB-p65_pS536_Caution`), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$`NF-kB-p65_pS536_Caution`, y = rppa_AHR_nonedit$`NF-kB-p65_pS536_Caution`,paired = F)

boxplot(rppa_AHR_edit$Collagen_VI,rppa_AHR_nonedit$Collagen_VI,outline = T)
stripchart(list(rppa_AHR_edit$Collagen_VI,rppa_AHR_nonedit$Collagen_VI), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$Collagen_VI,rppa_AHR_nonedit$Collagen_VI,paired = F)

########## cutoff as 0.6
AHR_nonedit = colnames(AHR_edit_rppa_selected[,AHR_edit_rppa_selected<0.5])
AHR_edit = setdiff(colnames(AHR_edit_rppa_selected),y = AHR_nonedit)

rppa_AHR_edit = arr[,AHR_edit]
rppa_AHR_edit = as.vector(as.character(rppa_AHR_edit[2,]))
rppa_AHR_edit = rppa_selected[rppa_AHR_edit,]

rppa_AHR_nonedit = arr[,AHR_nonedit]
rppa_AHR_nonedit = as.vector(as.character(rppa_AHR_nonedit[2,]))
rppa_AHR_nonedit = rppa_selected[rppa_AHR_nonedit,]

boxplot(rppa_AHR_edit$`beta-Catenin`,rppa_AHR_nonedit$`beta-Catenin`,outline = T)
stripchart(list(rppa_AHR_edit$`beta-Catenin`,rppa_AHR_nonedit$`beta-Catenin`), vertical = T, pch=20,method = "jitter",add = T,col ="red")
boxplot(rppa_AHR_edit$`beta-Catenin_pT41_S45`,rppa_AHR_nonedit$`beta-Catenin_pT41_S45`,outline = F)
stripchart(list(rppa_AHR_edit$`beta-Catenin_pT41_S45`,rppa_AHR_nonedit$`beta-Catenin_pT41_S45`), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$`beta-Catenin`, y = rppa_AHR_nonedit$`beta-Catenin`,paired = F)

boxplot(rppa_AHR_edit$Cyclin_D1,rppa_AHR_nonedit$Cyclin_D1,outline = T)
stripchart(list(rppa_AHR_edit$Cyclin_D1,rppa_AHR_nonedit$Cyclin_D1), vertical = T, pch=20,method = "jitter",add = T,col ="red")
boxplot(rppa_AHR_edit$Cyclin_B1,rppa_AHR_nonedit$Cyclin_B1,outline = T)
stripchart(list(rppa_AHR_edit$Cyclin_B1,rppa_AHR_nonedit$Cyclin_B1), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$Cyclin_B1, y = rppa_AHR_nonedit$Cyclin_B1,paired = F)


boxplot(rppa_AHR_edit$GSK3_pS9,rppa_AHR_nonedit$GSK3_pS9,outline = T)
stripchart(list(rppa_AHR_edit$GSK3_pS9,rppa_AHR_nonedit$GSK3_pS9), vertical = T, pch=20,method = "jitter",add = T,col ="red")

wilcox.test(x = rppa_AHR_edit$GSK3_pS9, y = rppa_AHR_nonedit$GSK3_pS9,paired = F)

boxplot(rppa_AHR_edit$`GSK3-alpha-beta_pS21_S9`,rppa_AHR_nonedit$`GSK3-alpha-beta_pS21_S9`,outline = T)
stripchart(list(rppa_AHR_edit$`GSK3-alpha-beta_pS21_S9`,rppa_AHR_nonedit$`GSK3-alpha-beta_pS21_S9`), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$`GSK3-alpha-beta_pS21_S9`, y = rppa_AHR_nonedit$`GSK3-alpha-beta_pS21_S9`,paired = F)

boxplot(rppa_AHR_edit$Snail_Caution,rppa_AHR_nonedit$Snail_Caution,outline = T)
stripchart(list(rppa_AHR_edit$Snail_Caution,rppa_AHR_nonedit$Snail_Caution), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$Snail_Caution, y = rppa_AHR_nonedit$Snail_Caution,paired = F)

boxplot(rppa_AHR_edit$ADAR1,rppa_AHR_nonedit$ADAR1,outline = T)
stripchart(list(rppa_AHR_edit$Snail_Caution,rppa_AHR_nonedit$Snail_Caution), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$ADAR1, y = rppa_AHR_nonedit$ADAR1,paired = F)

boxplot(rppa_AHR_edit$`NF-kB-p65_pS536_Caution`,rppa_AHR_nonedit$`NF-kB-p65_pS536_Caution`,outline = T)
stripchart(list(rppa_AHR_edit$`NF-kB-p65_pS536_Caution`,rppa_AHR_nonedit$`NF-kB-p65_pS536_Caution`), vertical = T, pch=20,method = "jitter",add = T,col ="red")
wilcox.test(x = rppa_AHR_edit$`NF-kB-p65_pS536_Caution`, y = rppa_AHR_nonedit$`NF-kB-p65_pS536_Caution`,paired = F)

############### correlations with some staff

indices = order(as.vector(as.numeric(AHR_edit_rppa_selected[1,])),decreasing = F)

AHR_edit_rppa_selected_ordered = AHR_edit_rppa_selected[indices]
rppa_selected_ordered = rppa_selected[indices,]

plot(as.vector(as.numeric(AHR_edit_rppa_selected_ordered[1,])),rppa_selected_ordered$ADAR1)




########### with cutoff == 0,editVsnonedit
vec = c()
diff = c()
for( i in 1:ncol(rppa_AHR_edit)){
  a = wilcox.test(x = as.vector(as.numeric(rppa_AHR_edit[,i])),y =as.vector(as.numeric(rppa_AHR_nonedit[,i])) )$p.value
  vec = c(vec,a)
  diff =c(diff, mean(as.numeric(rppa_AHR_edit[,i])) - mean(as.vector(as.numeric(rppa_AHR_nonedit[,i]))))
  
}

bool = vec<0.05
p.correct = p.adjust(vec,method = "bonferroni")
bool_correct=  p.correct<0.05
sum(bool_correct)
vec = vec[bool]
diff  =diff[bool]
sig_rppa_ahr_edit = rppa_AHR_edit[bool]
sig_rppa_ahr_nonedit = rppa_AHR_nonedit[bool]
setwd("/home/ahadli/Desktop/300719_sunum")
pdf(file = "290719_RPPA_array_differential_protein_exp_AHReditvsAHRnonedit.pdf")
for(i in 1:ncol(sig_rppa_ahr_edit)){
  boxplot(as.vector(as.numeric(sig_rppa_ahr_edit[,i])),as.vector(as.numeric(sig_rppa_ahr_nonedit[,i])),outline = ,names = c("Edited AHR", "Not Edited AHR"))
  stripchart(list(as.vector(as.numeric(sig_rppa_ahr_edit[,i])),as.vector(as.numeric(sig_rppa_ahr_nonedit[,i]))), vertical = T, pch=20,method = "jitter",add = T,col ="red")
  title(paste0(colnames(sig_rppa_ahr_edit)[i], " \np = ",vec[i]))
}
dev.off()
    diff_nonedit_up_bool = diff<0
diff_edit_up_bool = diff>0
 
cat(paste0(colnames(sig_rppa_ahr_edit)[diff_edit_up_bool],collapse = "\n"))
