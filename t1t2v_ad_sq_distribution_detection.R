setwd('/media/ahadli/Data/Farid/AOG_lab/Alperin isleri/miRNA_analysis/october_new_project')
ad_cancer = read.csv("AD_primary_miRNA_RPM_log2.csv",row.names = 1, header = T)
sq_cancer = read.csv("SQ_primary_miRNA_RPM_log2.csv",row.names = 1, header = T)
ad_normal  = read.csv("AD_normal_miRNA_RPM_log2.csv",row.names = 1, header = T)
sq_normal = read.csv("SQ_normal_miRNA_RPM_log2.csv",row.names = 1, header = T)
t1 =  read.csv("T1_NSCLC_ttest.csv",row.names = 1, header = T)
t2 =  read.csv("t2_NSCLC_ttest.csv",row.names = 1, header = T)
v = read.csv("V_NSCLC_ttest.csv",row.names = 1, header = T)
t1 = t1[1:(length(t1)-4)]
t2 = t2[1:(length(t2)-4)]
v = v[1:(length(v)-4)]
t1_new_nms = substr(colnames(t1), 15, nchar(colnames(t1)))
t2_new_nms = substr(colnames(t2), 8, nchar(colnames(t2)))
v_new_nms = substr(colnames(v), 8, nchar(colnames(v)))
c_t1_ad = 0
for(i in 1: length(t1_new_nms)){
  for( j in 1:length(ad_cancer)){
    if(colnames(ad_cancer)[j]==t1_new_nms[i]){
      c_t1_ad = c_t1_ad +1
    }
  }
}
c_t2_ad = 0
for(i in 1: length(t2_new_nms)){
  for( j in 1:length(ad_cancer)){
    if(colnames(ad_cancer)[j]==t2_new_nms[i]){
      c_t2_ad = c_t2_ad +1
    }
  }
}
c_v_ad = 0
for(i in 1: length(v_new_nms)){
  for( j in 1:length(ad_cancer)){
    if(colnames(ad_cancer)[j]==v_new_nms[i]){
      c_v_ad = c_v_ad +1
    }
  }
}

c_t1_sq = 0
for(i in 1: length(t1_new_nms)){
  for( j in 1:length(sq_cancer)){
    if(colnames(sq_cancer)[j]==t1_new_nms[i]){
      c_t1_sq = c_t1_sq +1
    }
  }
}
c_t2_sq = 0
for(i in 1: length(t2_new_nms)){
  for( j in 1:length(sq_cancer)){
    if(colnames(sq_cancer)[j]==t2_new_nms[i]){
      c_t2_sq = c_t2_sq +1
    }
  }
}
c_v_sq = 0
for(i in 1: length(v_new_nms)){
  for( j in 1:length(sq_cancer)){
    if(colnames(sq_cancer)[j]==v_new_nms[i]){
      c_v_sq = c_v_sq +1
    }
  }
}

############# normal
n_t1_ad = 0
for(i in 1: length(t1_new_nms)){
  for( j in 1:length(ad_normal)){
    if(colnames(ad_normal)[j]==t1_new_nms[i]){
      n_t1_ad = n_t1_ad +1
    }
  }
}
n_t2_ad = 0
for(i in 1: length(t2_new_nms)){
  for( j in 1:length(ad_normal)){
    if(colnames(ad_normal)[j]==t2_new_nms[i]){
      n_t2_ad = n_t2_ad +1
    }
  }
}
n_v_ad = 0
for(i in 1: length(v_new_nms)){
  for( j in 1:length(ad_normal)){
    if(colnames(ad_normal)[j]==v_new_nms[i]){
      n_v_ad = n_v_ad +1
    }
  }
}

n_t1_sq = 0
for(i in 1: length(t1_new_nms)){
  for( j in 1:length(sq_normal)){
    if(colnames(sq_normal)[j]==t1_new_nms[i]){
      n_t1_sq = n_t1_sq +1
    }
  }
}
n_t2_sq = 0
for(i in 1: length(t2_new_nms)){
  for( j in 1:length(sq_normal)){
    if(colnames(sq_normal)[j]==t2_new_nms[i]){
      n_t2_sq = n_t2_sq +1
    }
  }
}
n_v_sq = 0
for(i in 1: length(v_new_nms)){
  for( j in 1:length(sq_normal)){
    if(colnames(sq_normal)[j]==v_new_nms[i]){
      n_v_sq = n_v_sq +1
    }
  }
}
