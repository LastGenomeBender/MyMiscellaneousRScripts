setwd("D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/log2FC2_RPM_nonprefiltered")
### find count matrix for ADmol genes
install.packages("rlist")
library(rlist)
t1_ad_sig_mol=read.delim("T1_NormalvsTumor_AD_Molecular_DESeq2.csv", sep=",", row.names = 1)
t1_mol_sig_ad_counts = t1_AD__mol_ttest[1:(length(t1_AD__mol_ttest)-2)]
colnames = c(rep(x = "AD", times =81), rep(x="norm",times=28))
colnames(t1_mol_sig_ad_counts)=colnames
for(i in 1:length(rownames(t1_ad_sig_mol))){
  for(j in 1:length(rownames(log2norm_AD_mol))){
    if(rownames(t1_ad_sig_mol)[i]==rownames(log2norm_AD_mol)[j]){
      t1_mol_sig_ad_counts = rbind(t1_mol_sig_ad_counts, log2norm_AD_mol[j,])
    }
  }
}

### standardize ADmol
rwstd_AD_mol = t1_mol_sig_ad_counts
for(i in 1:length(rownames(rwstd_AD_mol))){
  stdev_row = sd(data.matrix(rwstd_AD_mol[i,]))
  mean_row = mean(data.matrix((rwstd_AD_mol[i,])))
  for(j in 1:length(rwstd_AD_mol)){
    std_value = (rwstd_AD_mol[i,j] - mean_row)/stdev_row
    rwstd_AD_mol[i,j] = std_value
  }
}

###### find miRs in ADmol
miR_list = list()
logic_frame = data.matrix(rwstd_AD_mol>0.75)
logic_of_tumor = logic_frame[,1:81]
logic_of_norm = logic_frame[,82:109]
rsumm = rowSums(logic_of_norm) 
rsum = rsumm==0
logic_frame = logic_frame[rsum,]
while(sum(rowSums(logic_frame)) > 0){
  rowsums = rowSums(logic_frame)
  max = max(rowsums)
  locations = which(data.matrix(rowsums) %in% c(max))
  names = rownames(logic_frame)[locations[1]]
  miR_list= list.append(miR_list, names)
  logic_frame = logic_frame[,!logic_frame[locations[1],]]
}
colnames(logic_frame)
##### find cutoff values for AD miRNAs
cutoff_AD_molecular=data.frame(row.names= "cutoff value")
for(i in 1:length(miR_list)){
  row = log2norm_AD_mol[miR_list[[i]],]
  cutoff = max(row)
  cutoff_AD_molecular[1,i]= cutoff
  colnames(cutoff_AD_molecular)[i] = miR_list[[i]]
}
write.csv(cutoff_AD_molecular ,file="cutoff_miRNAs_AD_mol.csv")

#### find count matrix for SQmol genes
t1_SQ_sig_mol=read.delim("T1_NormalvsTumor_SQ_Molecular_DESeq2.csv", sep=",", row.names = 1)
t1_mol_sig_SQ_counts = t1_SQ__mol_ttest[1:(length(t1_SQ__mol_ttest)-2)]

colnames(t1_mol_sig_SQ_counts) = c(rep(x="SQ",times=56),rep(x="norm",times=22))
for(i in 1:length(rownames(t1_SQ_sig_mol))){
  for(j in 1:length(rownames(log2norm_SQ_mol))){
    if(rownames(t1_SQ_sig_mol)[i]==rownames(log2norm_SQ_mol)[j]){
      t1_mol_sig_SQ_counts = rbind(t1_mol_sig_SQ_counts, log2norm_SQ_mol[j,])
    }
  }
}


### standardize SQmol
rwstd_SQ_mol = t1_mol_sig_SQ_counts
for(i in 1:length(rownames(rwstd_SQ_mol))){
  stdev_row = sd(data.matrix(rwstd_SQ_mol[i,]))
  mean_row = mean(data.matrix((rwstd_SQ_mol[i,])))
  for(j in 1:length(rwstd_SQ_mol)){
    std_value = (rwstd_SQ_mol[i,j] - mean_row)/stdev_row
    rwstd_SQ_mol[i,j] = std_value
  }
}



### find miRNAs in 
miR_list = list()
logic_frame = data.matrix(rwstd_SQ_mol>0.75)
logic_of_tumor = logic_frame[,1:56]
logic_of_norm = logic_frame[,57:78]
rsumm = rowSums(logic_of_norm) 
rsum = rsumm==0
logic_frame = logic_frame[rsum,]
while(sum(rowSums(logic_frame)) > 0){
  rowsums = rowSums(logic_frame)
  max = max(rowsums)
  locations = which(data.matrix(rowsums) %in% c(max))
  names = rownames(logic_frame)[locations[1]]
  miR_list= list.append(miR_list, names)
  logic_frame = logic_frame[,!logic_frame[locations[1],]]
}
colnames(logic_frame)







###  find cutoff values for SQ miRNAs
cutoff_SQ_molecular=data.frame(row.names= "cutoff value")
for(i in 1:length(miR_list)){
  row = log2norm_SQ_mol[miR_list[[i]],]
  cutoff = max(row)
  cutoff_SQ_molecular[1,i]= cutoff
  colnames(cutoff_SQ_molecular)[i] = miR_list[[i]]
}
write.csv(cutoff_SQ_molecular ,file="cutoff_miRNAs_SQ_mol.csv")

