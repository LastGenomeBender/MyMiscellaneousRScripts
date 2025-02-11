setwd("/media/ahadli/Data/Farid/AOG_lab/COAD/Secil_microarray_020718")
x = read.csv("Untitled 1.csv",header = T,row.names=1)
a = rep(F,times = 20)
for (i in seq (1,20,by=2)){
  a[i]=T
}
preop = x[,a]
postop = x[,!a]
ttest = data.frame()
df = data.frame()
for (i in 1:nrow(x)){
  t=t.test(as.numeric(preop[i,]),as.numeric(postop[i,]),paired = T,mu = 0)
  d = mean(as.numeric(preop[i,]))- mean(as.numeric(postop[i,]))
  df =rbind(df,d)
  ttest = rbind(ttest,t$p.value)
}
stat = cbind(ttest,df)
rownames(stat)=rownames(x)
sel_stat= stat[stat[1]<0.01 & stat[2]>2,]
#sel_stat2 = sel_stat 
#adjustp = p.adjust(p = as.matrix(sel_stat2[1]),method = "BH", n = nrow(sel_stat2))
#sel_stat2[1]=adjustp
#sel_stat2 = sel_stat2[sel_stat2[,1]<0.05,]
##annot
anno=read.delim("Annot.adf.csv",stringsAsFactors = F,header = T,row.names = 1)
anno_selected= anno[rownames(sel_stat),]
bool= duplicated(anno_selected[,1])
sel_stat = sel_stat[!bool,]
anno_selected=anno_selected[!bool,]
preop = preop[rownames(sel_stat),]
postop = postop[rownames(sel_stat),]
rownames(sel_stat)=as.vector(t(anno_selected[1]))
rownames(preop)=as.vector(t(anno_selected[1]))
rownames(postop)=as.vector(t(anno_selected[1]))
sel_stat = sel_stat[rownames(sel_stat)!="1",]
write.csv(sel_stat,"signif_mirs.csv")
############# Lassotime

library(glmnet)
lasso_br = preop[rownames(sel_stat),]
lasso_norm = postop[rownames(sel_stat),]
lasso_data = cbind(lasso_br,lasso_norm)
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