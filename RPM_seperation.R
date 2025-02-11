setwd("D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/log2FC2_RPM_nonprefiltered")

install.packages("stringdist")
install.packages("evd")
library(stringdist)
library(RecordLinkage)
AD_normal_miRNA_RPM=RPMmerge(mypath = "D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/LUAD normal miRNA/miRNA")
AD_primary_miRNA_RPM = RPMmerge(mypath = "D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/LUAD_primarytumor/miRNA")
SQ_normal_miRNA_RPM = RPMmerge(mypath = "D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/LUSC_mormal_miRNA/l/miRNA")
SQ_Primary_miRNA_RPM = RPMmerge(mypath = "D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/LUASC_primary/miRNA")

######### elemination of duplicates
AD_normal_miRNA_RPM = AD_normal_miRNA_RPM[ ,(!duplicated(colnames(AD_normal_miRNA_RPM)))]
AD_primary_miRNA_RPM = AD_primary_miRNA_RPM[, (!duplicated(colnames(AD_primary_miRNA_RPM)))]
SQ_normal_miRNA_RPM = SQ_normal_miRNA_RPM[, (!duplicated(colnames(SQ_normal_miRNA_RPM)))]
SQ_Primary_miRNA_RPM = SQ_Primary_miRNA_RPM[, (!duplicated(colnames(SQ_Primary_miRNA_RPM)))]

AD_normal_miRNA_RPM = log2(AD_normal_miRNA_RPM)
AD_primary_miRNA_RPM = log2(AD_primary_miRNA_RPM)
SQ_normal_miRNA_RPM = log2(SQ_normal_miRNA_RPM)
SQ_Primary_miRNA_RPM = log2(SQ_Primary_miRNA_RPM)

write.csv(AD_normal_miRNA_RPM,file="AD_normal_miRNA_RPM_log2.csv")
write.csv(AD_primary_miRNA_RPM,file="AD_primary_miRNA_RPM_log2.csv")
write.csv(SQ_normal_miRNA_RPM, file="SQ_normal_miRNA_RPM_log2.csv")
write.csv(SQ_Primary_miRNA_RPM, file="SQ_primary_miRNA_RPM_log2.csv")
tum_all = data.frame(AD_primary_miRNA_RPM, SQ_Primary_miRNA_RPM, check.names = F)




###### Separation of adeno primary into groups
colnames(t1_AD_mol)
zata = substr(colnames(tum_all),1,nchar(colnames(t1_AD_mol)))
t1_AD_mol_RPM = data.frame(row.names = rownames(t1_AD_mol))
for (i in 1 : length(t1_AD_mol)){
  for(j  in 1: length(zata)){
    if(stringdist(colnames(t1_AD_mol)[i],zata[j])==0){
      t1_AD_mol_RPM = cbind(t1_AD_mol_RPM, tum_all[j])
    }
  }
}

v_AD_mol
v_AD_mol_RPM = data.frame(row.names = rownames(v_AD_mol))
for (i in 1 : length(v_AD_mol)){
  for(j  in 1: length(zata)){
    if(colnames(v_AD_mol)[i]==zata[j]){
      v_AD_mol_RPM = cbind(v_AD_mol_RPM, tum_all[j])
    }
  }
}
v_AD_mol_RPM=v_AD_mol_RPM[,!duplicated(substr(colnames(v_AD_mol_RPM),1,nchar(colnames(t1_AD_mol))))]

t1_SQ_mol
t1_SQ_mol_RPM = data.frame(row.names = rownames(t1_SQ_mol))
for (i in 1 : length(t1_SQ_mol)){
  for(j  in 1: length(zata)){
    if(stringdist(colnames(t1_SQ_mol)[i],zata[j])==0){
      t1_SQ_mol_RPM = cbind(t1_SQ_mol_RPM, tum_all[j])
    }
  }
}

v_SQ_mol
v_SQ_mol_RPM = data.frame(row.names = rownames(v_SQ_mol))
for (i in 1 : length(v_SQ_mol)){
  for(j  in 1: length(zata)){
    if(stringdist(colnames(v_SQ_mol)[i],zata[j])==0){
      v_SQ_mol_RPM = cbind(v_SQ_mol_RPM, tum_all[j])
    }
  }
}

### normal seperation
t1_AD_n_mol
colnames(t1_AD_n_mol) = gsub("[.]","-",colnames(t1_AD_n_mol))
t1_AD_n_mol_RPM = data.frame(row.names= rownames(v_AD_mol))
for(i in 1:length(t1_AD_n_mol)){
  for(j in 1:length(AD_normal_miRNA_RPM)){
    if(colnames(t1_AD_n_mol)[i]==colnames(AD_normal_miRNA_RPM)[j]){
      t1_AD_n_mol_RPM = cbind(t1_AD_n_mol_RPM, AD_normal_miRNA_RPM[j] )
    }
  }
}

v_AD_n_mol
colnames(v_AD_n_mol) = gsub("[.]","-",colnames(v_AD_n_mol))
v_AD_n_mol_RPM = data.frame(row.names= rownames(v_AD_mol))
for(i in 1:length(v_AD_n_mol)){
  for(j in 1:length(AD_normal_miRNA_RPM)){
    if(colnames(v_AD_n_mol)[i]==colnames(AD_normal_miRNA_RPM)[j]){
      v_AD_n_mol_RPM = cbind(v_AD_n_mol_RPM, AD_normal_miRNA_RPM[j] )
    }
  }
}

t1_SQ_n_mol
colnames(t1_SQ_n_mol) = gsub("[.]","-",colnames(t1_SQ_n_mol))
t1_SQ_n_mol_RPM = data.frame(row.names= rownames(v_SQ_mol))
for(i in 1:length(t1_SQ_n_mol)){
  for(j in 1:length(SQ_normal_miRNA_RPM)){
    if(colnames(t1_SQ_n_mol)[i]==colnames(SQ_normal_miRNA_RPM)[j]){
      t1_SQ_n_mol_RPM = cbind(t1_SQ_n_mol_RPM, SQ_normal_miRNA_RPM[j] )
    }
  }
}

v_SQ_n_mol
colnames(v_SQ_n_mol) = gsub("[.]","-",colnames(v_SQ_n_mol))
v_SQ_n_mol_RPM = data.frame(row.names= rownames(v_SQ_mol))
for(i in 1:length(v_SQ_n_mol)){
  for(j in 1:length(SQ_normal_miRNA_RPM)){
    if(colnames(v_SQ_n_mol)[i]==colnames(SQ_normal_miRNA_RPM)[j]){
      v_SQ_n_mol_RPM = cbind(v_SQ_n_mol_RPM, SQ_normal_miRNA_RPM[j] )
    }
  }
}


write.csv(t1_AD_mol_RPM, file = "t1_AD_mol_RPM_log2.csv")
write.csv(v_AD_mol_RPM, file="v_AD_mol_RPM_log2.csv")
write.csv(t1_AD_n_mol_RPM, file="t1_AD_n_mol_RPM_log2.csv")
write.csv(v_AD_n_mol_RPM, file = "v_AD_n_mol_RPM_log2.csv")
write.csv(t1_SQ_mol_RPM, file = "t1_SQ_mol_RPM_log2.csv")
write.csv(v_SQ_mol_RPM, file="v_SQ_mol_RPM_log2.csv")
write.csv(t1_SQ_n_mol_RPM, file="t1_SQ_n_mol_RPM_log2.csv")
write.csv(v_SQ_n_mol_RPM, file = "v_SQ_n_mol_RPM_log2.csv")

### ttest of AD mol
logfc = data.frame()
ad_RPM_ttest = data.frame()
for(i in 1:dim(t1_AD_mol_RPM)[1]){
  test = t.test(t1_AD_mol_RPM[i,], t1_AD_n_mol_RPM[i,])
  log = mean(data.matrix(t1_AD_mol_RPM[i,])) - mean(data.matrix(t1_AD_n_mol_RPM[i,]))
  logfc = rbind(logfc, log)
  ad_RPM_ttest = rbind(ad_RPM_ttest, test$p.value)
}
colnames(ad_RPM_ttest) = "p value"
colnames(logfc) = "Log2FC"
t1_AD__mol_ttest=  cbind(t1_AD_mol_RPM, t1_AD_n_mol_RPM, ad_RPM_ttest ,logfc)
t1_AD__mol_ttest = t1_AD__mol_ttest[t1_AD__mol_ttest$"p value"<0.05,]
t1_AD__mol_ttest = t1_AD__mol_ttest[t1_AD__mol_ttest$"Log2FC">2,]
t1_AD__mol_ttest = t1_AD__mol_ttest[order(t1_AD__mol_ttest$"p value"),]

### ttest of SQ mol
logfc = data.frame()
SQ_RPM_ttest = data.frame()
for(i in 1:dim(t1_SQ_mol_RPM)[1]){
  test = t.test(t1_SQ_mol_RPM[i,], t1_SQ_n_mol_RPM[i,])
  log = mean(data.matrix(t1_SQ_mol_RPM[i,])) - mean(data.matrix(t1_SQ_n_mol_RPM[i,]))
  logfc = rbind(logfc, log)
  SQ_RPM_ttest = rbind(SQ_RPM_ttest, test$p.value)
}
colnames(SQ_RPM_ttest) = "p value"
colnames(logfc) = "Log2FC"
t1_SQ__mol_ttest=  cbind(t1_SQ_mol_RPM, t1_SQ_n_mol_RPM, SQ_RPM_ttest ,logfc)
t1_SQ__mol_ttest = t1_SQ__mol_ttest[t1_SQ__mol_ttest$"p value"<0.05,]
t1_SQ__mol_ttest = t1_SQ__mol_ttest[t1_SQ__mol_ttest$"Log2FC">2,]
t1_SQ__mol_ttest = t1_SQ__mol_ttest[order(t1_AD__mol_ttest$"p value"),]

write.csv(t1_SQ__mol_ttest,file="t1_SQ_mol_RPM_ttest.csv")
write.csv(t1_AD__mol_ttest,file="t1_AD_mol_RPM_ttest.csv")
