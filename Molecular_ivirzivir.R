AC_molecular = read.csv(file = "AC_molecular.csv", stringsAsFactors = F,check.names = F)
SCC_molecular = read.csv(file = "SCC_molecular.csv", check.names = F)
mixed_molecular = read.csv(file = "mixed_molecular.csv", check.names = F)
combinedrnames = colnames(AD_primary_miRNA)
combinedrnames = c(combinedrnames, colnames(SQ_Primary_miRNA))
for(i in 1:length(combinedrnames)){
  combinedrnames[i] = substr(combinedrnames[i], 1, (nchar(combinedrnames[i])-1))
}
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
