setwd("/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/Lasso lung")

install.packages("stringdist")
install.packages("evd")
library(stringdist)
library(RecordLinkage)
AD_normal_miRNA_RPM=RPMmerge(mypath = "/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/LUAD normal miRNA/miRNA")
AD_primary_miRNA_RPM = RPMmerge(mypath = "/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/LUAD_primarytumor/miRNA")
SQ_normal_miRNA_RPM = RPMmerge(mypath = "/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/LUSC_mormal_miRNA/l/miRNA")
SQ_Primary_miRNA_RPM = RPMmerge(mypath = "/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/LUASC_primary/miRNA")

######### elemination of duplicates
AD_normal_miRNA_RPM = AD_normal_miRNA_RPM[ ,(!duplicated(colnames(AD_normal_miRNA_RPM)))]
AD_primary_miRNA_RPM = AD_primary_miRNA_RPM[, (!duplicated(colnames(AD_primary_miRNA_RPM)))]
SQ_normal_miRNA_RPM = SQ_normal_miRNA_RPM[, (!duplicated(colnames(SQ_normal_miRNA_RPM)))]
SQ_Primary_miRNA_RPM = SQ_Primary_miRNA_RPM[, (!duplicated(colnames(SQ_Primary_miRNA_RPM)))]

AD_normal_miRNA_RPM = log2(AD_normal_miRNA_RPM)
AD_primary_miRNA_RPM = log2(AD_primary_miRNA_RPM)
SQ_normal_miRNA_RPM = log2(SQ_normal_miRNA_RPM)
SQ_Primary_miRNA_RPM = log2(SQ_Primary_miRNA_RPM)

write.csv(AD_normal_miRNA_RPM,file="Guven_AD_normal_miRNA_RPM_log2.csv")
write.csv(AD_primary_miRNA_RPM,file="Guven_AD_primary_miRNA_RPM_log2.csv")
write.csv(SQ_normal_miRNA_RPM, file="Guven_SQ_normal_miRNA_RPM_log2.csv")
write.csv(SQ_Primary_miRNA_RPM, file="Guven_SQ_primary_miRNA_RPM_log2.csv")
#######
BRCA = cbind(AD_primary_miRNA_RPM, SQ_Primary_miRNA_RPM)
norm= cbind(AD_normal_miRNA_RPM, SQ_normal_miRNA_RPM)




x=sample(1:3,length(BRCA), replace = T)
BRCA[length(rownames(BRCA))+1, ] = x
n1=1
n2=1
n3=1
cnames = colnames(BRCA)
s1 = c("")
s2 = c("")
s3 = c("")
t1_AD_n_mi = data.frame(row.names = rownames(BRCA)) 
t2_AD_n_mi = data.frame(row.names = rownames(BRCA))
v_AD_n_mi = data.frame(row.names = rownames(BRCA))
for(i in 1:length(BRCA)){
  determ=BRCA[dim(BRCA)[1],i]
  if(determ==1){
    t1_AD_n_mi[,n1] = BRCA[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else if(determ==2){
    t2_AD_n_mi[,n2] = BRCA[,i]
    n2=n2+1
    s2 = c(s2, cnames[i])
  } 
  else{
    v_AD_n_mi[,n3]= BRCA[,i]
    n3=n3+1
    s3= c(s3, cnames[i])
  }
}
s1 = s1[2:length(s1)]
s2 = s2[2:length(s2)]
s3 = s3[2:length(s3)]
colnames(t1_AD_n_mi) = s1
colnames(t2_AD_n_mi) = s2
colnames(v_AD_n_mi) =s3
t1_AD_n_mi=t1_AD_n_mi[-(length(rownames(t1_AD_n_mi))),]
t2_AD_n_mi=t2_AD_n_mi[-length(rownames(t2_AD_n_mi)),]
v_AD_n_mi=v_AD_n_mi[-(length(rownames(v_AD_n_mi))),]

write.csv(t1_AD_n_mi, file="Guven_T1_Lung_TCGA.csv")
write.csv(t2_AD_n_mi,file="Guven_T2_Lung_TCGA.csv")
write.csv(v_AD_n_mi,file="Guven_V_Lung_TCGA.csv")

##### normal separation
x=sample(1:3,length(norm), replace = T)
norm[length(rownames(norm))+1, ] = x
n1=1
n2=1
n3=1
cnames = colnames(norm)
s1 = c("")
s2 = c("")
s3 = c("")
t1_AD_n_mi = data.frame(row.names = rownames(norm)) 
t2_AD_n_mi = data.frame(row.names = rownames(norm))
v_AD_n_mi = data.frame(row.names = rownames(norm))
for(i in 1:length(norm)){
  determ=norm[dim(norm)[1],i]
  if(determ==1){
    t1_AD_n_mi[,n1] = norm[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else if(determ==2){
    t2_AD_n_mi[,n2] = norm[,i]
    n2=n2+1
    s2 = c(s2, cnames[i])
  } 
  else{
    v_AD_n_mi[,n3]= norm[,i]
    n3=n3+1
    s3= c(s3, cnames[i])
  }
}
s1 = s1[2:length(s1)]
s2 = s2[2:length(s2)]
s3 = s3[2:length(s3)]
colnames(t1_AD_n_mi) = s1
colnames(t2_AD_n_mi) = s2
colnames(v_AD_n_mi) =s3
t1_AD_n_mi=t1_AD_n_mi[-(length(rownames(t1_AD_n_mi))),]
t2_AD_n_mi=t2_AD_n_mi[-length(rownames(t2_AD_n_mi)),]
v_AD_n_mi=v_AD_n_mi[-(length(rownames(v_AD_n_mi))),]
write.csv(t1_AD_n_mi, file="T1_norm_TCGA.csv")
write.csv(t2_AD_n_mi,file="T2_norm_TCGA.csv")
write.csv(v_AD_n_mi,file="V_norm_TCGA.csv")
