##### correlation matrix of significant miRNAs

### AD

ad_sig = read.csv("t1_AD_mol_RPM_only_sig_mirnascsv.csv", header = T, row.names = 1)
correlation_mat = cor(data.matrix(t(ad_sig)))
mt_mirnas = read.csv("farid_RPM_cutoff_miRNAs_AD_mol.csv", header = T, row.names = 1, check.names = F)
mt_sig = ad_sig[colnames(mt_mirnas),]
mor_mat = cor(data.matrix(t(mt_sig)))
for (i in 2:length(mor_mat)){
  mor_mat[c(1:(i-1)),i] = 0
}
rwstd_AD_mol = ad_sig
for(i in 1:length(rownames(rwstd_AD_mol))){
  stdev_row = sd(data.matrix(rwstd_AD_mol[i,]))
  mean_row = mean(data.matrix((rwstd_AD_mol[i,])))
  for(j in 1:length(rwstd_AD_mol)){
    std_value = (rwstd_AD_mol[i,j] - mean_row)/stdev_row
    rwstd_AD_mol[i,j] = std_value
  }
}
heatmap(data.matrix(t(rwstd_AD_mol)))
