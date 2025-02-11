#### ROC for t1 AD RPM
install.packages("ROCR")
install.packages("ROSE")
library(ROSE)
library(ROCR)
my_AD_mat =  t1_AD__mol_ttest[1:(length(t1_AD__mol_ttest)-2)]
mir_cutoff_AD_mol
mir_count_AD_mat = my_AD_mat[colnames(mir_cutoff_AD_mol),]
mir_count_AD_mat[11,] = c(rep(1,times=81),rep(0, times=28))
rownames(mir_count_AD_mat)[11] = "condition" 
tr_mir_ad_count = data.frame(t(mir_count_AD_mat))
v_AD_RPM
mir_ad_v_rpm = v_AD_RPM[colnames(mir_cutoff_AD_mol),]
mir_ad_v_rpm[11,] = c(rep(1,times=92),rep(0, times=18)) 
rownames(mir_ad_v_rpm)[11] = "condition"
colnames(mir_ad_v_rpm)= make.names(names = rep(x="x",times=110),unique = T)
tr_AD_v_rpm = data.frame(t(mir_ad_v_rpm))
t1_RPM_logit = glm(condition~hsa.mir.148a+hsa.mir.9.1+hsa.mir.141+hsa.mir.577+hsa.mir.29b.1+hsa.mir.429+hsa.mir.142+hsa.mir.183+hsa.mir.301a+hsa.mir.708,data=tr_mir_ad_count,family=binomial(link="logit"))
AD_predict = predict(t1_RPM_logit, tr_AD_v_rpm )
ROC = roc.curve(tr_AD_v_rpm$condition,AD_predict)
