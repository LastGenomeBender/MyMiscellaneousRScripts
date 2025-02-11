library(glmnet)
gse6 = read.csv("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project/commons_between_febbits_gse61741.csv", header = T,row.names = 1, check.names = F,stringsAsFactors = F)
gse6 = gse6[1:(length(gse6)-4)]
cancer = gse6[1:73]
normal = gse6[74:length(gse6)]
a = sample(1:2, length(cancer), replace = T)
c_t = cancer[a==1]
c_v = cancer[a==2]
a = sample(1:2, length(cancer), replace = T)
n_t = normal[a==1]
n_v = normal[a==2]
setwd("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project/LASSO_Febits")
write.csv(c_t, "cancer_testcsv.csv")
write.csv(c_v, "cancer_val.csv")
write.csv(n_t, "norm_test.csv")
write.csv(n_v,"norm_val.csv")
###### make lasso model
lasso_data1=cbind(c_t,n_t)
lasso_data2 = t(cbind(lasso_data1))
ad1 = rep(x = 1, times=ncol(c_t))
ad2 = rep(x= 0, times =ncol(n_t))
response = c(ad1, ad2)
lasso_data2 = cbind(lasso_data2,response) 
lasso_data2 = as.data.frame(lasso_data2)
colnames(lasso_data2)[ncol(lasso_data2)] = "resp"
x <- model.matrix(resp~.,lasso_data2)
#conT2ert class to numerical T2ariable

cv.out <- cv.glmnet(x,as.numeric(t(lasso_data2[ncol(lasso_data2)])),alpha=1,family="binomial",type.measure = "auc" )
#plot result
plot(cv.out)
lambda_min <- cv.out$lambda.min
cv.out$lambda
cv.out
#best T2alue of lambda
lambda_1se <- cv.out$lambda.1se
#regression coefficients
write.csv(data.matrix(coef(cv.out,s=lambda_1se)), "overall_coeffs.csv")

lasso_prob <- predict(cv.out,newx = x,s=lambda_1se,type="response")

ROC = roc.curve(lasso_data2$resp,lasso_prob)
ROC$auc
title("                                             
      
      AUC is  0.8612811")


########### validation lasso
lasso_data1=cbind(c_v,n_v)
lasso_data2 = t(cbind(lasso_data1))
ad1 = rep(x = 1, times=ncol(c_v))
ad2 = rep(x= 0, times =ncol(n_v))
response = c(ad1, ad2)
lasso_data2 = cbind(lasso_data2,response) 
lasso_data2 = as.data.frame(lasso_data2)
colnames(lasso_data2)[ncol(lasso_data2)] = "resp"
x <- model.matrix(resp~.,lasso_data2)
lasso_prob <- predict(cv.out,newx = x,s=lambda_1se,type="response")
ROC = roc.curve(lasso_data2$resp,lasso_prob)
ROC$auc
title("                                             
      
    AUC is 0.8135135")


###### 
gse2 = read.csv("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project/commons_between_febbits_gse24709.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
gse2 = gse2[1:(length(gse2)-4)]
c_t = gse2[1:23]
n_t = gse2[24:length(gse2)]
lass_mirs=read.csv("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project/LASSO_Febits/GSE61741/overall_coeffs.csv", header = T, stringsAsFactors = F, check.names = F)
lass_mirs = lass_mirs[lass_mirs$`1`>0,]
lass_mirs[,1] = substr(lass_mirs[,1], 2, (nchar(lass_mirs[,1])-1))
c_t = c_t[lass_mirs[,1],]
n_t = n_t[lass_mirs[,1],]
############ model
lasso_data1=cbind(c_t,n_t)
lasso_data2 = t(cbind(lasso_data1))
ad1 = rep(x = 1, times=ncol(c_t))
ad2 = rep(x= 0, times =ncol(n_t))
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
write.csv(data.matrix(coef(cv.out,s=lambda_1se)), "overall_coeffs.csv")

lasso_prob <- predict(cv.out,newx = x,s=lambda_1se,type="response")

ROC = roc.curve(lasso_data2$resp,lasso_prob)
ROC$auc
title("                                             
      
      AUC is  0.9290618")
