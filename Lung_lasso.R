install.packages("glmnet")
library(glmnet)
lasso_br = read.csv("Guven_T1_Lung_TCGA.csv",row.names=1)
lasso_norm = read.csv("T1_norm_TCGA.csv",row.names=1)
lasso_data = cbind(lasso_br,lasso_norm)
significant = read.csv("Guven_T1_Lung_ttest_fc3.csv",row.names = 1)
lasso_data = lasso_data[rownames(significant),]
head(lasso_data)
lasso_data2 = t(lasso_data)
ad1 = rep(x = 1, times=ncol(lasso_br))
ad2 = rep(x= 0, times =ncol(lasso_norm))
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
#best T2alue of lambda
lambda_1se <- cv.out$lambda.1se
#regression coefficients
coef(cv.out,s=lambda_1se)
class(cv.out$glmnet.fit)
######### ROC
t3_brca =read.csv("Guven_T2_Lung_TCGA.csv", row.names=1)
t3_norm = read.csv("T2_norm_TCGA.csv", row.names=1)
t3_data = cbind(t3_brca, t3_norm)
significant = read.csv("T1_BRCA_ttest.csv",row.names = 1)
t3_data = t3_data[rownames(significant),]
t3_data = as.data.frame(t(t3_data))

ad1 = rep(x = 1, times=ncol(t3_brca))
ad2 = rep(x= 0, times =ncol(t3_norm))
response = c(ad1, ad2)
t3_data = cbind(t3_data,response) 
colnames(t3_data)[ncol(t3_data)]="resp"
x_test <- model.matrix(resp~.,t3_data)

lasso_prob <- predict(cv.out,newx = x_test,s=lambda_1se,type="response")

ROC = roc.curve(t3_data$resp,lasso_prob)
ROC$auc
title("                                             
      
    AUC is 0.9877549")




#######
install.packages("InformationValue") 
library(InformationValue)
a = youdensIndex(t3_data$resp,lasso_prob)


#########
V_brca =read.csv("Guven_V_Lung_TCGA.csv", row.names=1)
V_norm = read.csv("V_norm_TCGA.csv", row.names=1)
V_data = cbind(V_brca, V_norm)
significant = read.csv("T1_BRCA_ttest.csv",row.names = 1)
V_data = V_data[rownames(significant),]
V_data = as.data.frame(t(V_data))

ad1 = rep(x = 1, times=ncol(V_brca))
ad2 = rep(x= 0, times =ncol(V_norm))
response = c(ad1, ad2)
V_data = cbind(V_data,response) 
colnames(V_data)[ncol(V_data)]="resp"
x_test <- model.matrix(resp~.,V_data)
lasso_prob_V3 <- predict(cv.out,newx = x_test,s=lambda_1se,type="response")
ROC_V3 = roc.curve(V_data$resp,lasso_prob_V3)
ROC_V3$auc
title("                                       
      
    AUC is 0.9850198")
glm_predict <- rep(0,nrow(V_data))
glm_predict[lasso_prob_V3>0.8053612] <- 1
table(pred=glm_predict,true=V_data$resp)
mean(glm_predict==V_data$resp)
