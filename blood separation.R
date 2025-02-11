x=sample(1:2,length(SQ_Primary_miRNA), rePlace = T)
SQ_Primary_miRNA[length(rownames(SQ_Primary_miRNA))+1, ] = x
n1=1
n2=1
cnames = colnames(SQ_Primary_miRNA)
s1 = c("")
s2 = c("")
t1_SQ_n_mi = data.frame(row.names = rownames(SQ_Primary_miRNA)) 
t2_SQ_n_mi = data.frame(row.names = rownames(SQ_Primary_miRNA))
for(i in 1:length(SQ_Primary_miRNA)){
  determ=SQ_Primary_miRNA[dim(SQ_Primary_miRNA)[1],i]
  if(determ==1){
    t1_SQ_n_mi[,n1] = SQ_Primary_miRNA[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else if(determ==2){
    t2_SQ_n_mi[,n2] = SQ_Primary_miRNA[,i]
    n2=n2+1
    s2 = c(s2, cnames[i])
  } 
}
s1 = s1[2:length(s1)]
s2 = s2[2:length(s2)]
colnames(t1_SQ_n_mi) = s1
colnames(t2_SQ_n_mi) = s2
t1_SQ_n_mi=t1_SQ_n_mi[-(length(rownames(t1_SQ_n_mi))),]
t2_SQ_n_mi=t2_SQ_n_mi[-length(rownames(t2_SQ_n_mi)),]
write.csv(t1_SQ_n_mi, file="T1_SQuamous_blood_Primary_mi_TCGA.csv")
write.csv(t2_SQ_n_mi,file="T2_SQuamous_blood_Primary_mi_TCGA.csv")