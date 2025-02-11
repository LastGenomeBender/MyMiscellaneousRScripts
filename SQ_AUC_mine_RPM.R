### ROC with ROCR
install.packages("nnet")
install.packages("AUC")
library(AUC)
library(nnet)
library(ROCR)
my_SQ_mat =  t1_SQ__mol_ttest[1:(length(t1_SQ__mol_ttest)-2)]
mir_cutoff_SQ_mol
mir_count_SQ_mat = my_SQ_mat[colnames(mir_cutoff_SQ_mol),]
mir_count_SQ_mat[(dim(mir_count_SQ_mat)[1]+1),] = c(rep(1,times=56),rep(0, times=22))
rownames(mir_count_SQ_mat)[dim(mir_count_SQ_mat)[1]] = "condition" 
tr_mir_SQ_count = data.frame(t(mir_count_SQ_mat))
v_SQ_RPM
mir_SQ_v_rpm = v_SQ_RPM[colnames(mir_cutoff_SQ_mol),]
mir_SQ_v_rpm[(dim(mir_SQ_v_rpm)[1]+1),] = c(rep(1,times=60),rep(0, times=23)) 
rownames(mir_SQ_v_rpm)[dim(mir_SQ_v_rpm)[1]] = "condition"
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
plot(ROC_SQ_RPM, main ="AUC=0.9782609")
abline(a=0, b=1)
auc = performance(SQ_predict,"auc")

