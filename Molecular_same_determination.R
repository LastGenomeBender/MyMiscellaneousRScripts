AC_molecular = read.csv(file = "AC_molecular.csv",check.names = F)
SCC_molecular = read.csv(file = "SCC_molecular.csv", check.names = F)
mixed_molecular = read.csv(file = "mixed_molecular.csv",check.names = F)
combinedrnames = colnames(AD_primary_miRNA)
combinedrnames = c(combinedrnames, colnames(SQ_Primary_miRNA))
for(i in 1:length(combinedrnames)){
  combinedrnames[i] = substr(combinedrnames[i], 1, (nchar(combinedrnames[i])-1))
}
combined_mi_counts = data.frame(AD_primary_miRNA, SQ_Primary_miRNA)
colnames(combined_mi_counts) = combinedrnames
verified_mixed_molecular = c("")
verified_SCC_molecular = c("")
verified_AC_molecular = c("")
for(i in 1:length(AC_molecular)){
  for(j in 1:length(combinedrnames)){
    if (colnames(AC_molecular)[i]==combinedrnames[j]){
      verified_AC_molecular= c(verified_AC_molecular, colnames(AC_molecular)[i])
    }
  }
}
for(i in 1:length(SCC_molecular)){
  for(j in 1:length(combinedrnames)){
    if (colnames(SCC_molecular)[i]==combinedrnames[j]){
      verified_SCC_molecular= c(verified_SCC_molecular, colnames(SCC_molecular)[i])
    }
  }
}
for(i in 1:length(mixed_molecular)){
  for(j in 1:length(combinedrnames)){
    if (colnames(mixed_molecular)[i]==combinedrnames[j]){
      verified_mixed_molecular= c(verified_mixed_molecular, colnames(mixed_molecular)[i])
    }
  }
}
verified_AC_molecular = verified_AC_molecular[-1]
verified_SCC_molecular = verified_SCC_molecular[-1]
verified_mixed_molecular = verified_mixed_molecular[-1]
verified_mixed_molecular = verified_AC_molecular[!duplicated(verified_mixed_molecular)]
verified_AC_molecular = verified_AC_molecular[!duplicated(verified_AC_molecular)]
verified_SCC_molecular = verified_SCC_molecular[!duplicated(verified_SCC_molecular)]
common_counts_AC_molecular = data.frame(row.names = rownames(SQ_Primary_miRNA))
common_counts_SCC_molecular = data.frame(row.names = rownames(SQ_Primary_miRNA))
common_counts_mixed_molecular = data.frame(row.names = rownames(SQ_Primary_miRNA))
n = 1
cnames= c("")
for(i in 1:length(AC_molecular)){
  for(j in 1: length(combined_mi_counts)){
    if(colnames(AC_molecular)[i]==colnames(combined_mi_counts)[j]){
      common_counts_AC_molecular = data.frame(common_counts_AC_molecular, combined_mi_counts[,j])
      cnames = c(cnames, colnames(combined_mi_counts)[j])
    } 
  }
}
cnames = cnames[-1]
colnames(common_counts_AC_molecular)=cnames
cnames= c("")
for(i in 1:length(SCC_molecular)){
  for(j in 1: length(combined_mi_counts)){
    if(colnames(SCC_molecular)[i]==colnames(combined_mi_counts)[j]){
      common_counts_SCC_molecular = data.frame(common_counts_SCC_molecular, combined_mi_counts[,j])
      cnames = c(cnames, colnames(combined_mi_counts)[j])
    } 
  }
}
cnames = cnames[-1]
colnames(common_counts_SCC_molecular)= cnames
common_counts_AC_molecular = common_counts_AC_molecular[ ,!duplicated(colnames(common_counts_AC_molecular))]
common_counts_SCC_molecular = common_counts_SCC_molecular[ ,!duplicated(colnames(common_counts_SCC_molecular))]
