## ACC molecular seperation
common_counts_AC_molecular = common_counts_AC_molecular[1:1881,]
x=sample(1:2,length(common_counts_AC_molecular), replace = T)
common_counts_AC_molecular[length(rownames(common_counts_AC_molecular))+1, ] = x
n1=1
n3=1
cnames = colnames(common_counts_AC_molecular)
s1 = c("")
s3 = c("")
t1_AD_mol = data.frame(row.names = rownames(common_counts_AC_molecular)) 
v_AD_mol = data.frame(row.names = rownames(common_counts_AC_molecular))
for(i in 1:length(common_counts_AC_molecular)){
  determ=common_counts_AC_molecular[dim(common_counts_AC_molecular)[1],i]
  if(determ==1){
    t1_AD_mol[,n1] = common_counts_AC_molecular[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else{
    v_AD_mol[,n3]= common_counts_AC_molecular[,i]
    n3=n3+1
    s3= c(s3, cnames[i])
  }
}
s1 = s1[2:length(s1)]
s3 = s3[2:length(s3)]
colnames(t1_AD_mol) = s1
colnames(v_AD_mol) =s3
t1_AD_mol=t1_AD_mol[-(length(rownames(t1_AD_mol))),]
v_AD_mol=v_AD_mol[-(length(rownames(v_AD_mol))),]
write.csv(t1_AD_mol, file="T1_Adeno_molecular_mi_TCGA.csv")
write.csv(v_AD_mol,file="V_Adeno_molecular_mi_TCGA.csv")

## SCC molecular seperation
x=sample(1:2,length(common_counts_SCC_molecular), replace = T)
common_counts_SCC_molecular = common_counts_SCC_molecular[1:1881,]
common_counts_SCC_molecular[length(rownames(common_counts_SCC_molecular))+1, ] = x
n1=1
n3=1
cnames = colnames(common_counts_SCC_molecular)
s1 = c("")
s3 = c("")
t1_SQ_mol = data.frame(row.names = rownames(common_counts_SCC_molecular)) 
v_SQ_mol = data.frame(row.names = rownames(common_counts_SCC_molecular))
for(i in 1:(length(common_counts_SCC_molecular)-1)){
  determ=common_counts_SCC_molecular[dim(common_counts_SCC_molecular)[1],i]
  if(determ==1){
    t1_SQ_mol[,n1] = common_counts_SCC_molecular[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else{
    v_SQ_mol[,n3]= common_counts_SCC_molecular[,i]
    n3=n3+1
    s3= c(s3, cnames[i])
  }
}
s1 = s1[2:length(s1)]
s3 = s3[2:length(s3)]
colnames(t1_SQ_mol) = s1
colnames(v_SQ_mol) =s3
t1_SQ_mol=t1_SQ_mol[-(length(rownames(t1_SQ_mol))),]
v_SQ_mol=v_SQ_mol[-(length(rownames(v_SQ_mol))),]
write.csv(t1_SQ_mol, file="T1_Squamous_molecular_mi_TCGA.csv")
write.csv(v_SQ_mol,file="V_Squamous_molecular_mi_TCGA.csv")

###Adeno normal seperation
t1_AD_n_mol = data.frame(row.names = rownames(AC_molecular))
t1_AD_n_mol = data.frame(t1_AD_n_mi, t2_AD_n_mi[,1:(round(length(t2_AD_n_mi)/2))])
v_AD_n_mol = data.frame(row.names = rownames(AC_molecular))
v_AD_n_mol = data.frame(t2_AD_n_mi[ ,(round(length(t2_AD_n_mi)/2)+1):length(t2_AD_n_mi)], v_AD_n_mi )

### squamous normal seperation
t1_SQ_n_mol = data.frame(row.names = rownames(SCC_molecular))
t1_SQ_n_mol = data.frame(t1_SQ_n_mi, t2_SQ_n_mi[,1:(round(length(t2_SQ_n_mi)/2))])
v_SQ_n_mol = data.frame(row.names = rownames(SCC_molecular))
v_SQ_n_mol = data.frame(t2_SQ_n_mi[ ,(round(length(t2_SQ_n_mi)/2)+1):length(t2_SQ_n_mi)], v_SQ_n_mi )



write.csv(t1_SQ_n_mol, file="T1_Squamous_normal_molecular_mi_TCGA.csv")
write.csv(t1_AD_n_mol, file="T1_Adeno_normal_molecular_mi_TCGA.csv")
write.csv(v_AD_n_mol, file="V_Adeno_normal_molecular_mi_TCGA.csv")
write.csv(v_SQ_n_mol, file="V_Squamous_normal_molecular_mi_TCGA.csv")


