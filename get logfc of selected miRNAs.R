AD_molecular = data.frame(t1_AD_mol, v_AD_mol)
SQ_molecular = data.frame(t1_SQ_mol, v_SQ_mol)
write.csv(AD_molecular, file="farid_AD_molecular_raw.csv")
write.csv(SQ_molecular, file="farid_SQ_molecular_raw.csv")

ttest_AD_mol = read.csv(file="T1_NormalvsTumor_AD_Molecular_DESeq2.csv", row.names = 1, check.names = F)
logFC = data.frame(row.names = "logfc")
for(i in 1:length(mir_cutoff_AD_mol)){
  for(j in 1:dim(ttest_AD_mol)[1]){
    if(colnames(mir_cutoff_AD_mol)[i]==rownames(ttest_AD_mol)[j]){
      logFC = cbind(logFC, ttest_AD_mol[j,2])
    }
  }
}

ttest_SQ_mol = read.csv(file="T1_NormalvsTumor_SQ_Molecular_DESeq2.csv", row.names = 1, check.names = F)
logFC_SQ = data.frame(row.names = "logfc")
for(i in 1:length(mir_cutoff_SQ_mol)){
  for(j in 1:dim(ttest_SQ_mol)[1]){
    if(colnames(mir_cutoff_SQ_mol)[i]==rownames(ttest_SQ_mol)[j]){
      logFC_SQ = cbind(logFC_SQ, ttest_SQ_mol[j,2])
    }
  }
}
