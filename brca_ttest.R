t1_brca = read.csv("T1_BRCA_TCGA.csv", row.names = 1)
t2_brca = read.csv("T2_BRCA_TCGA.csv", row.names = 1)
V_brca = read.csv("V_BRCA_TCGA.csv", row.names = 1)

t1_norm = read.csv("T1_norm_TCGA.csv", row.names = 1)
t2_norm = read.csv("T2_norm_TCGA.csv", row.names = 1)
V_norm = read.csv("V_norm_TCGA.csv", row.names = 1)

####################
logfc = data.frame()
ad_RPM_ttest = data.frame()
for(i in 1:dim(t1_brca)[1]){
  test = t.test(t1_brca[i,], t1_norm[i,])
  log = mean(data.matrix(t1_brca[i,])) - mean(data.matrix(t1_norm[i,]))
  logfc = rbind(logfc, log)
  ad_RPM_ttest = rbind(ad_RPM_ttest, test$p.value)
}
colnames(ad_RPM_ttest) = "p value"
colnames(logfc) = "Log2FC"
t1_brca_ttest=  cbind(t1_brca, t1_norm, ad_RPM_ttest ,logfc)
t1_brca_ttest = t1_brca_ttest[t1_brca_ttest$"p value"<0.05,]
t1_brca_ttest = t1_brca_ttest[t1_brca_ttest$"Log2FC">2,]
t1_brca_ttest = t1_brca_ttest[order(t1_brca_ttest$"p value"),]

write.csv(file = "T1_BRCA_ttest.csv", t1_brca_ttest)
