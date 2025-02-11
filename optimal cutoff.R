install.packages("OptimalCutpoints")
library(OptimalCutpoints)
library(cutpointr)
setwd("D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/log2FC2_RPM_nonprefiltered")
my_AD_mat =  read.csv("t1_AD_mol_RPM_ttest.csv", header=T, row.names = 1)
my_AD_mat = my_AD_mat[1:(length(my_AD_mat)-2)]
my_AD_mat = data.frame(t(my_AD_mat))
my_AD_mat = data.frame(t(my_AD_mat))
mir_cutoff_AD_mol = read.csv("farid_RPM_cutoff_miRNAs_AD_mol.csv", header = T,row.names = 1)
mir_count_AD_mat = my_AD_mat[colnames(mir_cutoff_AD_mol),]
mir_count_AD_mat[11,] = c(rep(1,times=81),rep(0, times=28))
rownames(mir_count_AD_mat)[11] = "condition" 
tr_mir_ad_count = data.frame(t(mir_count_AD_mat))
v_AD_RPM
v_AD =  read.csv("v_AD_mol_RPM_log2.csv", header = T,row.names = 1)
v_AD_norm =  read.csv("v_AD_n_mol_RPM_log2.csv", header = T,row.names = 1)
v_AD_RPM = cbind(v_AD,v_AD_norm)
v_AD_RPM = data.frame(t(v_AD_RPM))
v_AD_RPM = data.frame(t(v_AD_RPM))
mir_ad_v_rpm = v_AD_RPM[colnames(mir_cutoff_AD_mol),]
mir_ad_v_rpm[11,] = c(rep(1,times=92),rep(0, times=18)) 
rownames(mir_ad_v_rpm)[11] = "condition"
colnames(mir_ad_v_rpm)= make.names(names = rep(x="x",times=110),unique = T)
tr_AD_v_rpm = data.frame(t(mir_ad_v_rpm))
t1_RPM_logit = multinom(condition~.,data=tr_mir_ad_count)

cutoff = multi_cutpointr(tr_mir_ad_count, class="condition" , method = oc_youden_normal)

cutoff = optimal.cutpoints(X = hsa.mir.148a+hsa.mir.9.1+hsa.mir.141+hsa.mir.577+hsa.mir.29b.1+hsa.mir.429+hsa.mir.142+hsa.mir.183+hsa.mir.301a+hsa.mir.708~condition, 0, methods =  "Youden", data = tr_mir_ad_count, control = control.cutpoints(generalized.Youden = T,costs.benefits.Youden = T))
AD_predict = predict(t1_RPM_logit,tr_AD_v_rpm, type = 'prob' )
AD_predict = prediction(AD_predict,tr_AD_v_rpm$condition)
eval = performance(AD_predict, "acc")
max_acc = which.max(slot(eval,"y.values")[[1]])
acc = slot(eval,"y.values")[[1]][max_acc]
cutoff_of_max_acc = slot(eval,"x.values")[[1]][max_acc]
ROC_AD_RPM = performance(AD_predict, "tpr","fpr")
print(opt.cut(ROC_AD_RPM, AD_predict))
plot(ROC_AD_RPM, main ="AUC 0.986413")
abline(a=0, b=1)
auc = performance(AD_predict,"auc")
auc


p <- predict(t1_RPM_logit, tr_AD_v_rpm)
tab <- table(p, tr_AD_v_rpm$condition )
tab
