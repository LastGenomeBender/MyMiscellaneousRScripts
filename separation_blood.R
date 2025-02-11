AD_normal_miRNA= read.csv("n_ad_blood_selected_counts.csv", header =T,row.names = 1, stringsAsFactors = F)
AD_primary_miRNA = read.csv("ad_blood_selected_counts.csv", header =T,row.names = 1, stringsAsFactors = F)
SQ_normal_miRNA =  read.csv("n_sq_blood_selected_counts.csv", header =T,row.names = 1, stringsAsFactors = F)
SQ_Primary_miRNA = read.csv("sq_blood_selected_counts.csv", header =T,row.names = 1, stringsAsFactors = F)
######### elemination of duplicates
AD_normal_miRNA = AD_normal_miRNA[ ,(!duplicated(colnames(AD_normal_miRNA)))]
AD_primary_miRNA = AD_primary_miRNA[, (!duplicated(colnames(AD_primary_miRNA)))]
SQ_normal_miRNA = SQ_normal_miRNA[, (!duplicated(colnames(SQ_normal_miRNA)))]
SQ_Primary_miRNA = SQ_Primary_miRNA[, (!duplicated(colnames(SQ_Primary_miRNA)))]

#######################
# seperation into T1 T2 and V
# sepertion of AD_normal
x=sample(1:2,length(AD_normal_miRNA), replace = T)
AD_normal_miRNA[length(rownames(AD_normal_miRNA))+1, ] = x
n1=1
n2=1
cnames = colnames(AD_normal_miRNA)
s1 = c("")
s2 = c("")
t1_AD_n_mi = data.frame(row.names = rownames(AD_normal_miRNA)) 
t2_AD_n_mi = data.frame(row.names = rownames(AD_normal_miRNA))
for(i in 1:length(AD_normal_miRNA)){
  determ=AD_normal_miRNA[dim(AD_normal_miRNA)[1],i]
  if(determ==1){
    t1_AD_n_mi[,n1] = AD_normal_miRNA[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else if(determ==2){
    t2_AD_n_mi[,n2] = AD_normal_miRNA[,i]
    n2=n2+1
    s2 = c(s2, cnames[i])
  } 
}
s1 = s1[2:length(s1)]
s2 = s2[2:length(s2)]
colnames(t1_AD_n_mi) = s1
colnames(t2_AD_n_mi) = s2
t1_AD_n_mi=t1_AD_n_mi[-(length(rownames(t1_AD_n_mi))),]
t2_AD_n_mi=t2_AD_n_mi[-length(rownames(t2_AD_n_mi)),]
write.csv(t1_AD_n_mi, file="T1_Adeno_blood_normal_mi_TCGA.csv")
write.csv(t2_AD_n_mi,file="T2_Adeno_blood_normal_mi_TCGA.csv")

##### SQ normal
x=sample(1:2,length(SQ_normal_miRNA), replace = T)
SQ_normal_miRNA[length(rownames(SQ_normal_miRNA))+1, ] = x
n1=1
n2=1
cnames = colnames(SQ_normal_miRNA)
s1 = c("")
s2 = c("")
t1_SQ_n_mi = data.frame(row.names = rownames(SQ_normal_miRNA)) 
t2_SQ_n_mi = data.frame(row.names = rownames(SQ_normal_miRNA))
for(i in 1:length(SQ_normal_miRNA)){
  determ=SQ_normal_miRNA[dim(SQ_normal_miRNA)[1],i]
  if(determ==1){
    t1_SQ_n_mi[,n1] = SQ_normal_miRNA[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else if(determ==2){
    t2_SQ_n_mi[,n2] = SQ_normal_miRNA[,i]
    n2=n2+1
    s2 = c(s2, cnames[i])
  } 
}
s1 = s1[2:length(s1)]
s2 = s2[2:length(s2)]
colnames(t1_SQ_n_mi) = s1
colnames(t2_SQ_n_mi) = s2
t1_SQ_n_mi=t1_SQ_n_mi[-(length(rownames(t1_SQ_n_mi))),]
t2_SQ_n_mi=t2_SQ_n_mi[-length(rownames(t2_SQ_n_mi)),]
write.csv(t1_SQ_n_mi, file="T1_squamous_blood_normal_mi_TCGA.csv")
write.csv(t2_SQ_n_mi,file="T2_squamous_blood_normal_mi_TCGA.csv")

####### AD primary
x=sample(1:2,length(AD_primary_miRNA), replace = T)
AD_primary_miRNA[length(rownames(AD_primary_miRNA))+1, ] = x
n1=1
n2=1
cnames = colnames(AD_primary_miRNA)
s1 = c("")
s2 = c("")
t1_AD_n_mi = data.frame(row.names = rownames(AD_primary_miRNA)) 
t2_AD_n_mi = data.frame(row.names = rownames(AD_primary_miRNA))
for(i in 1:length(AD_primary_miRNA)){
  determ=AD_primary_miRNA[dim(AD_primary_miRNA)[1],i]
  if(determ==1){
    t1_AD_n_mi[,n1] = AD_primary_miRNA[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else if(determ==2){
    t2_AD_n_mi[,n2] = AD_primary_miRNA[,i]
    n2=n2+1
    s2 = c(s2, cnames[i])
  } 
}
s1 = s1[2:length(s1)]
s2 = s2[2:length(s2)]
colnames(t1_AD_n_mi) = s1
colnames(t2_AD_n_mi) = s2
t1_AD_n_mi=t1_AD_n_mi[-(length(rownames(t1_AD_n_mi))),]
t2_AD_n_mi=t2_AD_n_mi[-length(rownames(t2_AD_n_mi)),]
write.csv(t1_AD_n_mi, file="T1_Adeno_blood_primary_mi_TCGA.csv")
write.csv(t2_AD_n_mi,file="T2_Adeno_blood_primary_mi_TCGA.csv")

#### squamous primary
x=sample(1:2,length(SQ_Primary_miRNA), replace = T)
SQ_Primary_miRNA[length(rownames(SQ_Primary_miRNA))+1, ] = x
n1=1
n2=1
cnames = colnames(SQ_Primary_miRNA)
s1 = c("")
s2 = c("")
t1_SQ_n_mi = data.frame(row.names = rownames(SQ_Primary_miRNA)) 
t2_SQ_n_mi = data.frame(row.names = rownames(SQ_Primary_miRNA))
for(i in 1:length(SQ_Primary_miRNA)){
  determ=SQ_Primary_miRNA[dim(SQ_Primary_miRNA)[1],i]
  if(determ==1){
    t1_SQ_n_mi[,n1] = SQ_Primary_miRNA[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else if(determ==2){
    t2_SQ_n_mi[,n2] = SQ_Primary_miRNA[,i]
    n2=n2+1
    s2 = c(s2, cnames[i])
  } 
}
s1 = s1[2:length(s1)]
s2 = s2[2:length(s2)]
colnames(t1_SQ_n_mi) = s1
colnames(t2_SQ_n_mi) = s2
t1_SQ_n_mi=t1_SQ_n_mi[-(length(rownames(t1_SQ_n_mi))),]
t2_SQ_n_mi=t2_SQ_n_mi[-length(rownames(t2_SQ_n_mi)),]
write.csv(t1_SQ_n_mi, file="T1_SQuamous_blood_Primary_mi_TCGA.csv")
write.csv(t2_SQ_n_mi,file="T2_SQuamous_blood_Primary_mi_TCGA.csv")
