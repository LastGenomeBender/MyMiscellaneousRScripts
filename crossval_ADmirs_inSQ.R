### ROC with ROCR
install.packages("nnet")
install.packages("AUC")
library(AUC)
library(nnet)
library(ROCR)
my_SQ_mat =  cbind(t1_SQ_mol_RPM,t1_SQ_n_mol_RPM)
mir_cutoff_SQ_mol = read.csv(file="farid_RPM_cutoff_miRNAs_AD_mol.csv", row.names = 1,header = T,check.names = F)
mir_count_SQ_mat = my_SQ_mat[colnames(mir_cutoff_SQ_mol),]
mir_count_SQ_mat[11,] = c(rep(1,times=56),rep(0, times=22))
rownames(mir_count_SQ_mat)[11] = "condition" 
tr_mir_SQ_count = data.frame(t(mir_count_SQ_mat))
my_v_SQ_RPM = cbind(v_SQ_mol_RPM,v_SQ_n_mol_RPM)
mir_SQ_v_rpm = my_v_SQ_RPM[colnames(mir_cutoff_SQ_mol),]
mir_SQ_v_rpm[11,] = c(rep(1,times=60),rep(0, times=23)) 
rownames(mir_SQ_v_rpm)[11] = "condition"
colnames(mir_SQ_v_rpm)= make.names(names = rep(x="x",times=83),unique = T)
tr_SQ_v_rpm = data.frame(t(mir_SQ_v_rpm))
t1_RPM_logit = multinom(condition~.,data=tr_mir_SQ_count)
SQ_predict = predict(t1_RPM_logit,tr_SQ_v_rpm, type = 'prob' )
SQ_predict = prediction(SQ_predict,tr_SQ_v_rpm$condition)
eval = performance(SQ_predict, "acc")
max_acc = which.max(slot(eval,"y.values")[[1]])
acc = slot(eval,"y.values")[[1]][max_acc]
cutoff_of_max_acc = slot(eval,"x.values")[[1]][max]
ROC_SQ_RPM = performance(SQ_predict, "tpr","fpr")
plot(ROC_SQ_RPM, main ="AUC 0.9757246")
abline(a=0, b=1)
auc = performance(SQ_predict,"auc")
auc
