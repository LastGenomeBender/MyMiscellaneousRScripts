setwd("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project/GSE64591/")
table = read.delim2("GSE64591_normalized.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F)
table$ID_REF[duplicated(table$ID_REF)]
rownames(table) = make.unique(table$ID_REF)
table = table[-1]
table = data.matrix(table)
table = data.frame(table)
table = log2(table)
norm = table[1:100]
cancer = table[101:(length(table))]
### sep cancer into two groups 
rnd=sample(1:2,length(cancer), replace = T)
cancer_test = cancer[rnd==2]        # cancer_test = cancer
cancer_val = cancer[rnd==1]
######## sep norm
rnd2 = sample(1:2,length(norm), replace = T) 
norm_test = norm[rnd==2]    #cacner_test
norm_val = norm[rnd2 == 1]

write.csv(cancer_test, "cancer_testcsv.csv")
write.csv(cancer_val, "cancer_val.csv")
write.csv(norm_test, "norm_test.csv")
write.csv(norm_val,"norm_val.csv")
#### ttest for the test set
# pval = c()
# BH_pval = c()
# logFC=c()
# for(i in 1:nrow(cancer_test)){
#   test = t.test(cancer_test[i,], norm_test[i,])
#   fc = mean(data.matrix(cancer_test[i,]))- mean(data.matrix(norm_test[i,]))
#   pval = c(pval, test$p.value)
#   logFC = c(logFC, fc)
# }
# BH_pval = p.adjust(pval, method = "BH")
# ttest_matrix = cbind(cancer_test, norm_test,pval, BH_pval, logFC)
# sig_matrix = ttest_matrix[ttest_matrix$BH_pval<0.05,]
# sig_matrix_ordered = sig_matrix[order(sig_matrix$logFC,decreasing = T),]
# sig_matrix_ordered = sig_matrix_ordered[sig_matrix_ordered$logFC>0,]
# write.csv(ttest_matrix,"ttest_matrix.csv")
# write.csv(sig_matrix_ordered, "significant_ttest.csv")
#################### ttest for all
pval = c()
BH_pval = c()
logFC=c()
for(i in 1:nrow(cancer)){
  test = t.test(cancer[i,], norm[i,])
  fc = mean(data.matrix(cancer[i,]))- mean(data.matrix(norm[i,]))
  pval = c(pval, test$p.value)
  logFC = c(logFC, fc)
}
BH_pval = p.adjust(pval, method = "BH")
ttest_matrix = cbind(cancer, norm,pval, BH_pval, logFC)
sig_matrix = ttest_matrix[ttest_matrix$BH_pval<0.05,]
sig_matrix_ordered = sig_matrix[order(sig_matrix$logFC,decreasing = T),]
############################# annotation
matrix = sig_matrix_ordered
matrix$newAnnot = substr(rownames(matrix), 1, (nchar(rownames(matrix))-7))
library(miRNAmeConverter)
nc = MiRNANameConverter() # Create MiRNANameConverter object
for(i in 1:nrow(matrix)){
  if(paste0(translateMiRNAName(nc,matrix$newAnnot[i])$v21.0, "a") == "a"){
    matrix$newAnnot[i] = "a"
  }
  else{
    matrix$newAnnot[i] = translateMiRNAName(nc,matrix$newAnnot[i])$v21.0
  }
}

sig_matrix_ordered = matrix
sig_matrix_ordered_1.5 = sig_matrix_ordered[sig_matrix_ordered$logFC>1.5,]
sig_matrix_ordered_1 = sig_matrix_ordered[sig_matrix_ordered$logFC>1,]
write.csv(ttest_matrix,"ttest_matrix_all.csv")
write.csv(sig_matrix_ordered_1.5, "significant_ttest_all_1.5.csv")
write.csv(sig_matrix_ordered_1, "significant_ttest_all_1.csv")

#### lasst time
lasso_data1=cbind(cancer_test,norm_test)
lasso_data1 = lasso_data1[rownames(sig_matrix_ordered),]
lasso_data2 = t(cbind(lasso_data1))
ad1 = rep(x = 1, times=ncol(cancer))
ad2 = rep(x= 0, times =ncol(norm))
response = c(ad1, ad2)
lasso_data2 = cbind(lasso_data2,response) 
lasso_data2 = as.data.frame(lasso_data2)
colnames(lasso_data2)[ncol(lasso_data2)] = "resp"
x <- model.matrix(resp~.,lasso_data2)
#conT2ert class to numerical T2ariable

cv.out <- cv.glmnet(x,as.numeric(t(lasso_data2[ncol(lasso_data2)])),alpha=1,family="binomial",type.measure = "mse" )
#plot result
plot(cv.out)
lambda_min <- cv.out$lambda.min
cv.out$lambda
cv.out
#best T2alue of lambda
lambda_1se <- cv.out$lambda.1se
#regression coefficients
coef(cv.out,s=lambda_1se)
class(cv.out$glmnet.fit)
######### ROC
# t3_brca =read.csv("Guven_T2_Lung_TCGA.csv", row.names=1)
# t3_norm = read.csv("T2_norm_TCGA.csv", row.names=1)
# t3_data = cbind(t3_brca, t3_norm)
# significant = read.csv("T1_BRCA_ttest.csv",row.names = 1)
# t3_data = t3_data[rownames(significant),]
# t3_data = as.data.frame(t(t3_data))
# 
# ad1 = rep(x = 1, times=ncol(t3_brca))
# ad2 = rep(x= 0, times =ncol(t3_norm))
# response = c(ad1, ad2)
# t3_data = cbind(t3_data,response) 
# colnames(t3_data)[ncol(t3_data)]="resp"
# x_test <- model.matrix(resp~.,t3_data)

lasso_prob <- predict(cv.out,newx = x,s=lambda_1se,type="response")

ROC = roc.curve(lasso_data2$resp,lasso_prob)
ROC$auc
title("                                             
      
      AUC is 0.8444444")


########### validation lasso
lasso_data1=cbind(cancer_val,norm_val)
lasso_data1 = lasso_data1[rownames(sig_matrix_ordered),]
lasso_data2 = t(cbind(lasso_data1))
ad1 = rep(x = 1, times=ncol(cancer_val))
ad2 = rep(x= 0, times =ncol(norm_val))
response = c(ad1, ad2)
lasso_data2 = cbind(lasso_data2,response) 
lasso_data2 = as.data.frame(lasso_data2)
colnames(lasso_data2)[ncol(lasso_data2)] = "resp"
x <- model.matrix(resp~.,lasso_data2)
lasso_prob <- predict(cv.out,newx = x,s=lambda_1se,type="response")
ROC = roc.curve(lasso_data2$resp,lasso_prob)
ROC$auc
title("                                             
      
    AUC is 0.6357708")
