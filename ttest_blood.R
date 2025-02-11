setwd("D:/Farid/AOG lab/Breast Cancer/GSE83270_normalized_matrix_wholoe_blood.txt")
x=read.table("GSE83270_normalized_matrix.txt", row.names = 1, stringsAsFactors = F, header = T, sep = "\t")
x=x[(rowSums(x=="N/A")==0),]
for(i in 2:ncol(x)){
  x[,i]=as.numeric(x[,i])
}
mirnames = x[1]
x= x[-1]
head(x)
x= log2(x)
logfc=data.frame()
ad_RPM_ttest = data.frame()
for(i in 1:dim(x)[1]){
  test = t.test(x[i,1:6], x[i,7:12])
  log = mean(data.matrix(x[i,1:6])) - mean(data.matrix(x[i,7:12]))
  logfc = rbind(logfc, log)
  ad_RPM_ttest = rbind(ad_RPM_ttest, test$p.value)
}
x= cbind(x,ad_RPM_ttest,logfc)
head(x)
bool = x[13]<0.05 & x[14]>2
x = x[bool,]
mirnames = mirnames[bool] 
rownames(x)= mirnames
write.csv(x,file ="blood_significant.csv" )
