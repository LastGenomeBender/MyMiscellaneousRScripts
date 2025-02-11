library(stringr)
###### annoting the miRNA IDs.
setwd("D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/BLOOD")
adeno_blood = read.csv("adeno_blood.csv", header = T, stringsAsFactors = F)
adeno_blood = adeno_blood[!duplicated(adeno_blood[1]),]
normal_blood = read.csv("normal_blood.csv", header = T, stringsAsFactors = F)
normal_blood = normal_blood[!duplicated(normal_blood[1]),]
squamous_blood = read.csv("squamous_blood.csv", header = T, stringsAsFactors = F)
squamous_blood = squamous_blood[!duplicated(squamous_blood[1]),]

setwd("D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/log2FC2_RPM_nonprefiltered")

sig_AD = read.csv("t1_AD_mol_RPM_ttest.csv", header = T, row.names = 1)
sig_SQ = read.csv("t1_SQ_mol_RPM_ttest.csv", header = T, row.names = 1)
setwd("D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/BLOOD")
ad_blood_selected_counts = data.frame()
rownames_short = as.character(adeno_blood[,1])
rad_blood = ""
for (i in 1: length(rownames(sig_AD))){
  for (j in 1: length(rownames_short)){
    if (rownames(sig_AD)[i]==tolower(substring(rownames_short[j],1,nchar(rownames(sig_AD)[i])))){
      ad_blood_selected_counts = rbind (ad_blood_selected_counts, adeno_blood[j, 2:length(adeno_blood)])
      rad_blood = c(rad_blood, rownames_short[j])
    }
  }
}
rad_blood = rad_blood[-1]
for(i in 1:length(rad_blood)){
  if(!is.na(str_locate(rad_blood[i],"/")[1])){
    rad_blood[i] = substring(rad_blood[i],1,(str_locate(rad_blood[i],"/")-1)[1])
  }
}
duplicate= duplicated(rad_blood)
rad_blood = rad_blood[!duplicate]
ad_blood_selected_counts = ad_blood_selected_counts[!duplicate,]
rownames(ad_blood_selected_counts) = rad_blood

sq_blood_selected_counts = data.frame()
rownames_short = as.character(squamous_blood[,1])
rsq_blood = ""
for (i in 1: length(rownames(sig_SQ))){
  for (j in 1: length(rownames_short)){
    if (rownames(sig_SQ)[i]==tolower(substring(rownames_short[j],1,nchar(rownames(sig_SQ)[i])))){
      sq_blood_selected_counts = rbind (sq_blood_selected_counts, squamous_blood[j, 2:length(squamous_blood)])
      rsq_blood = c(rsq_blood, rownames_short[j])
    }
  }
}
rsq_blood = rsq_blood[-1]
for(i in 1:length(rsq_blood)){
  if(!is.na(str_locate(rsq_blood[i],"/")[1])){
    rsq_blood[i] = substring(rsq_blood[i],1,(str_locate(rsq_blood[i],"/")-1)[1])
  }
}
duplicate= duplicated(rsq_blood)
rsq_blood = rsq_blood[!duplicate]
sq_blood_selected_counts = sq_blood_selected_counts[!duplicate,]
rownames(sq_blood_selected_counts) = rsq_blood


n_ad_blood_selected_counts = data.frame()
rownames_short = as.character(normal_blood[,1])
rn_n_ad_blood = ""
for (i in 1: length(rownames(sig_AD))){
  for (j in 1: length(rownames_short)){
    if (rownames(sig_AD)[i]==tolower(substring(rownames_short[j],1,nchar(rownames(sig_AD)[i])))){
      n_ad_blood_selected_counts = rbind (n_ad_blood_selected_counts, normal_blood[j, 2:length(normal_blood)])
      rn_n_ad_blood = c(rn_n_ad_blood, rownames_short[j])
    }
  }
}
rn_n_ad_blood = rn_n_ad_blood[-1]
for(i in 1:length(rn_n_ad_blood)){
  if(!is.na(str_locate(rn_n_ad_blood[i],"/")[1])){
    rn_n_ad_blood[i] = substring(rn_n_ad_blood[i],1,(str_locate(rn_n_ad_blood[i],"/")-1)[1])
  }
}
duplicate= duplicated(rn_n_ad_blood)
rn_n_ad_blood = rn_n_ad_blood[!duplicate]
n_ad_blood_selected_counts = n_ad_blood_selected_counts[!duplicate,]
rownames(n_ad_blood_selected_counts) = rn_n_ad_blood



n_sq_blood_selected_counts = data.frame()
rownames_short = as.character(normal_blood[,1])
rn_n_sq_blood = ""
for (i in 1: length(rownames(sig_SQ))){
  for (j in 1: length(rownames_short)){
    if (rownames(sig_SQ)[i]==tolower(substring(rownames_short[j],1,nchar(rownames(sig_SQ)[i])))){
      n_sq_blood_selected_counts = rbind (n_sq_blood_selected_counts, normal_blood[j, 2:length(normal_blood)])
      rn_n_sq_blood = c(rn_n_sq_blood, rownames_short[j])
    }
  }
}
rn_n_sq_blood = rn_n_sq_blood[-1]
for(i in 1:length(rn_n_sq_blood)){
  if(!is.na(str_locate(rn_n_sq_blood[i],"/")[1])){
    rn_n_sq_blood[i] = substring(rn_n_sq_blood[i],1,(str_locate(rn_n_sq_blood[i],"/")-1)[1])
  }
}
duplicate= duplicated(rn_n_sq_blood)
rn_n_sq_blood = rn_n_sq_blood[!duplicate]
n_sq_blood_selected_counts = n_sq_blood_selected_counts[!duplicate,]
rownames(n_sq_blood_selected_counts) = rn_n_sq_blood

write.csv(n_sq_blood_selected_counts,"n_sq_blood_selected_counts.csv")
write.csv(n_ad_blood_selected_counts,"n_ad_blood_selected_counts.csv")
write.csv(ad_blood_selected_counts,"ad_blood_selected_counts.csv")
write.csv(sq_blood_selected_counts,"sq_blood_selected_counts.csv")

