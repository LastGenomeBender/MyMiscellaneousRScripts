### ROC with ROCR
install.packages("nnet")
install.packages("AUC")
library(AUC)
library(nnet)
library(ROCR)
my_AD_mat =  cbind(t1_AD_mol_RPM,t1_AD_n_mol_RPM)
mir_cutoff_AD_mol = read.csv(file="mamed-separator-miRNAs-for-ADvsNorm-rpkm.csv", header=T, row.names = 1,stringsAsFactors = F)
mir_count_AD_mat = my_AD_mat[as.character(mir_cutoff_AD_mol[1,]),]
mir_count_AD_mat[12,] = c(rep(1,times=81),rep(0, times=28))
rownames(mir_count_AD_mat)[12] = "condition" 
tr_mir_ad_count = as.data.frame(t(mir_count_AD_mat))
my_v_AD_RPM = cbind(v_AD_mol_RPM,v_AD_n_mol_RPM)
mir_ad_v_rpm = v_AD_RPM[as.character(mir_cutoff_AD_mol[1,]),]
mir_ad_v_rpm[12,] = c(rep(1,times=92),rep(0, times=18)) 
rownames(mir_ad_v_rpm)[12] = "condition"
colnames(mir_ad_v_rpm)= make.names(names = rep(x="x",times=110),unique = T)
tr_AD_v_rpm = as.data.frame(t(mir_ad_v_rpm))
t1_RPM_logit = multinom(condition~.,data=tr_mir_ad_count)
AD_predict = predict(t1_RPM_logit,tr_AD_v_rpm, type = 'prob' )
AD_predict = prediction(AD_predict,tr_AD_v_rpm$condition)
eval = performance(AD_predict, "acc")
max_acc = which.max(slot(eval,"y.values")[[1]])
acc = slot(eval,"y.values")[[1]][max_acc]
cutoff_of_max_acc = slot(eval,"x.values")[[1]][max]
ROC_AD_RPM = performance(AD_predict, "tpr","fpr")
plot(ROC_AD_RPM, main ="AUC 0.9785628")
abline(a=0, b=1)
auc = performance(AD_predict,"auc")
auc
