mamed=read.csv("significant-miRNAs-from-GSE16512_mamed.csv",  header = T, stringsAsFactors = F)
mamed= mamed[(rowSums(!is.na(mamed[2])))>0,]
mamed[,1]=paste0("hsa-",mamed[,1])
a = c()
for(i in 1:nrow(mamed)){
  for(j in 1:nrow(sel_stat)){
    if(mamed[i,1]==rownames(sel_stat)[j]){
      a= c(a,rownames(sel_stat)[j])
    }
  }
}
for_preop=sel_stat[a,]
rownames(mamed)=mamed[,1]
for_sec = mamed[a,]
total = cbind(for_preop, for_sec)
tatal= total[,-3]
colnames(tatal)[1:2]=c("p_value_first","LogFC_first")
write.csv(tatal, "First_vs_Second_dataset_common.csv")
