setwd("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project")

t1_brca = read.csv("T1_NSCLC_TCGA.csv", row.names = 1)
t2_brca = read.csv("T2_NSCLC_TCGA.csv", row.names = 1)
V_brca = read.csv("V_NSCLC_TCGA.csv", row.names = 1)

t1_norm = read.csv("T1_norm_TCGA.csv", row.names = 1)
t2_norm = read.csv("T2_norm_TCGA.csv", row.names = 1)
V_norm = read.csv("V_norm_TCGA.csv", row.names = 1)

#################### for t1s
logfc = data.frame()
ad_RPM_ttest = data.frame()
for(i in 1:dim(t1_brca)[1]){
  test = t.test(t1_brca[i,], t1_norm[i,])
  log = mean(data.matrix(t1_brca[i,])) - mean(data.matrix(t1_norm[i,]))
  logfc = rbind(logfc, log)
  ad_RPM_ttest = rbind(ad_RPM_ttest, test$p.value)
}
colnames(ad_RPM_ttest) = "not_corrected_pvalue"
colnames(logfc) = "Log2FC"
BH = p.adjust(p=ad_RPM_ttest$not_corrected_pvalue, method = "BH")
benferonni = p.adjust(p=ad_RPM_ttest$not_corrected_pvalue, method = "bonferroni")
ad_RPM_ttest$benjaminihochberg_correction = BH
ad_RPM_ttest$bonferroni_corection = benferonni
canc=rep("cancer",length(t1_brca))
nor = rep("normal",length(t1_norm))
colnames(t1_brca)=paste0(canc,"/",colnames(t1_brca))
colnames(t1_norm)=paste0(nor,"/",colnames(t1_norm))
t1_brca_ttest =  cbind(t1_brca, t1_norm, ad_RPM_ttest ,logfc)
write.csv(t1_brca_ttest, "T1_NSCLC_ttest.csv")
t1_brca_ttest = t1_brca_ttest[t1_brca_ttest$"Log2FC">2,]
t1_brca_ttest = t1_brca_ttest[order(t1_brca_ttest$"Log2FC"),]
t1_brca_ttest_BH = t1_brca_ttest[t1_brca_ttest$benjaminihochberg_correction<0.05,]
write.csv(t1_brca_ttest_BH, "T1_NSCLC_BH_significant.csv")
t1_brca_ttest_bon = t1_brca_ttest[t1_brca_ttest$bonferroni_corection < 0.05,]
write.csv(t1_brca_ttest_bon, "T1_NSCLC_Bonferroni_significant.csv")

###### t2
logfc = data.frame()
ad_RPM_ttest = data.frame()
for(i in 1:dim(t2_brca)[1]){
  test = t.test(t2_brca[i,], t2_norm[i,])
  log = mean(data.matrix(t2_brca[i,])) - mean(data.matrix(t2_norm[i,]))
  logfc = rbind(logfc, log)
  ad_RPM_ttest = rbind(ad_RPM_ttest, test$p.value)
}
colnames(ad_RPM_ttest) = "not_corrected_pvalue"
colnames(logfc) = "Log2FC"
BH = p.adjust(p=ad_RPM_ttest$not_corrected_pvalue, method = "BH")
benferonni = p.adjust(p=ad_RPM_ttest$not_corrected_pvalue, method = "bonferroni")

ad_RPM_ttest$benjaminihochberg_correction = BH
ad_RPM_ttest$bonferroni_corection = benferonni
canc=rep("cancer",length(t2_brca))
nor = rep("normal",length(t2_norm))
colnames(t2_brca)=paste0(canc,"/",colnames(t2_brca))
colnames(t2_norm)=paste0(nor,"/",colnames(t2_norm))
t2_brca_ttest =  cbind(t2_brca, t2_norm, ad_RPM_ttest ,logfc)
write.csv(t2_brca_ttest, "t2_NSCLC_ttest.csv")
t2_brca_ttest = t2_brca_ttest[t2_brca_ttest$"Log2FC">2,]
t2_brca_ttest = t2_brca_ttest[order(t2_brca_ttest$"Log2FC"),]
t2_brca_ttest_BH = t2_brca_ttest[t2_brca_ttest$benjaminihochberg_correction<0.05,]
write.csv(t2_brca_ttest_BH, "t2_NSCLC_BH_significant.csv")
t2_brca_ttest_bon = t2_brca_ttest[t2_brca_ttest$bonferroni_corection   <0.05,]
write.csv(t2_brca_ttest_bon, "t2_NSCLC_Bonferroni_significant.csv")

################ v
logfc = data.frame()
ad_RPM_ttest = data.frame()
for(i in 1:dim(V_brca)[1]){
  test = t.test(V_brca[i,], V_norm[i,])
  log = mean(data.matrix(V_brca[i,])) - mean(data.matrix(V_norm[i,]))
  logfc = rbind(logfc, log)
  ad_RPM_ttest = rbind(ad_RPM_ttest, test$p.value)
}
colnames(ad_RPM_ttest) = "not_corrected_pvalue"
colnames(logfc) = "Log2FC"
BH = p.adjust(p=ad_RPM_ttest$not_corrected_pvalue, method = "BH")
benferonni = p.adjust(p=ad_RPM_ttest$not_corrected_pvalue, method = "bonferroni")
ad_RPM_ttest$benjaminihochberg_correction = BH
ad_RPM_ttest$bonferroni_corection = benferonni
canc=rep("cancer",length(V_brca))
nor = rep("normal",length(V_norm))
colnames(V_brca)=paste0(canc,"/",colnames(V_brca))
colnames(V_norm)=paste0(nor,"/",colnames(V_norm))
V_brca_ttest =  cbind(V_brca, V_norm, ad_RPM_ttest ,logfc)
write.csv(V_brca_ttest, "V_NSCLC_ttest.csv")
V_brca_ttest = V_brca_ttest[V_brca_ttest$"Log2FC">2,]
V_brca_ttest = V_brca_ttest[order(V_brca_ttest$"Log2FC"),]
V_brca_ttest_BH = V_brca_ttest[V_brca_ttest$benjaminihochberg_correction<0.05,]
write.csv(V_brca_ttest_BH, "V_NSCLC_BH_significant.csv")
V_brca_ttest_bon = V_brca_ttest[V_brca_ttest$bonferroni_corection <0.05,]
write.csv(V_brca_ttest_bon, "V_NSCLC_Bonferroni_significant.csv")
