install.packages("standardize")
library(standardize)
t1_squamous=read.delim("T1_Squamous_TCGA.csv", sep=",", row.names = 1)
t1_sq_100sig = read.delim("T1_ADvsSCC_100mostSig_TCGA.csv", sep=",", row.names = 1)
t2_sq_100siq = read.delim("T2_ADvsSCC_100mostSig_TCGA.csv", sep=",", row.names = 1)
common_genes = c("")
for(i in 1:100){
  for(j in 1:100){
    if(rownames(t1_sq_100sig)[i]==rownames(t2_sq_100siq)[j]){
      common_genes = c(common_genes, rownames(t1_sq_100sig)[i])
    }
  }
}
common_genes= common_genes[-1]
SQandADmedians=data.frame(row.names = common_genes)
SQmedians = c(0)
for (k in 1:71){
  for (l in 1:length(rownames(t1_squamous))){
    if(common_genes[k]==rownames(t1_squamous)[l]){
      med = median(data.matrix(t1_squamous[l,]))
      SQmedians = c(SQmedians, med)
    }
  }
}
SQmedians = SQmedians[-1]
SQandADmedians[,1] = SQmedians


####################
t1_adeno=read.delim("T1_Adeno_TCGA.csv", sep=",", row.names = 1)
ADmedians = c(0)
for (k in 1:71){
  for (l in 1:length(rownames(t1_adeno))){
    if(common_genes[k]==rownames(t1_adeno)[l]){
      med = median(data.matrix(t1_adeno[l,]))
      ADmedians = c(ADmedians, med)
    }
  }
}
ADmedians = ADmedians[-1]
SQandADmedians[,2] = ADmedians
colnames(SQandADmedians) = c("SQmedian","ADmedian")

############################################
V_squamous = read.delim("V_Squamous_TCGA.csv", sep=",", row.names = 1)
V_adeno = read.delim("V_Adeno_TCGA.csv", sep=",", row.names = 1)
V_combined = data.frame(V_squamous, V_adeno)
V_combined_71 = data.frame()
for (k in 1:71){
  for (l in 1:length(rownames(V_combined))){
    if(common_genes[k]==rownames(V_combined)[l]){
      V_combined_71 = rbind(V_combined_71,V_combined[l,])
    }
  }
}

###########################################
correlation_rnames = c("SQ(r)","ACC(r)")
correlation = data.frame(row.names = correlation_rnames)
for (k in 1:length(V_combined_71)){
    SQ_r= cor(data.matrix(V_combined_71[,k]),data.matrix(SQandADmedians[,1]))
    AD_r= cor(data.matrix(V_combined_71[,k]),data.matrix(SQandADmedians[,2]))
    correlation[,k] = c(SQ_r,AD_r)
}
colnames(correlation)= colnames(V_combined_71)
SQ = rep("SQ", times=192)
AD = rep("AD",times=192)
ActualPheno = c(SQ,AD)
correlation = rbind(correlation, ActualPheno)
rownames(correlation)
rownames(correlation)[3]="ActualPheno"
PredictedPheno = c("")
for(i in 1:384){
  if(correlation[1,i]>correlation[2,i]){
    PredictedPheno = c(PredictedPheno, "SQ") 
  }
  else{
    PredictedPheno = c(PredictedPheno, "AD")
  }
}
PredictedPheno= PredictedPheno[-1]
correlation = rbind(correlation, PredictedPheno)
dim(correlation)
rownames(correlation)[4] = "PredictedPheno"
rownames(correlation)
colnames(correlation)
write.csv(correlation, file="SQrandADr.csv")
SQcorAD=c();
SQcorSQ=c();
for(i in 1:384){
  if(correlation[3,i]=="SQ"){
    SQcorAD = c(SQcorAD, correlation[2,i])    
    SQcorSQ = c(SQcorSQ, correlation[1,i])
  }
}
ADcorAD=c();
ADcorSQ=c();
for(i in 1:384){
  if(correlation[3,i]=="AD"){
    ADcorAD = c(ADcorAD, correlation[2,i])    
    ADcorSQ = c(ADcorSQ, correlation[1,i])
  }
}
plot(SQcorAD,SQcorSQ, type="p",  pch=19, col="red", xlab="AC(r)",ylab="SCC(r)", ylim=c(-0.3,1.0),xlim=c(-0.3,1.0))
par(new=T)
plot(ADcorAD,ADcorSQ, type="p", pch=19,col="blue", xlab='', ylab='', axes=F,ylim=c(-0.3,1.0),xlim=c(-0.3,1.0))
legend(-0.1,-0.1, legend = c("Squamous Cell Carcinoma","Adenocarcinoma"),col=c("red","blue"),pch=19)
par(new=F)
oylesine= matrix(data=0, nrow=2,ncol=2)
PrevsAC =data.frame(data=oylesine,row.names=c("ACacttual", "SCCactual"))
for(i in 1:384){
  if(correlation[3,i]=="AD"&& correlation[4,i]=="AD"){
    PrevsAC[1,1] =PrevsAC[1,1]+1     
  }
  else if (correlation[3,i]=="AD"&& correlation[4,i]=="SQ"){
    PrevsAC[1,2] = PrevsAC[1,2]+1
  }
  else if (correlation[3,i]=="SQ"&&  correlation[4,i]=="AD"){
    PrevsAC[2,1]= PrevsAC[2,1]+1
  }
  else{
    PrevsAC[2,2]= PrevsAC[2,2]+1
  }
    
}
colnames(PrevsAC)=c("ACpredicted","SCCpredicted")
ol= matrix(data=0, ncol=2,nrow=2)
SpecificitySensitivity=data.frame(row.names = c("AC","SCC"), data = ol)
      SpecificitySensitivity[1,1] = PrevsAC[1,1]/(PrevsAC[1,1]+PrevsAC[1,2]) *100
      SpecificitySensitivity[2,1] = PrevsAC[2,2]/(PrevsAC[2,2]+PrevsAC[2,1]) *100
      SpecificitySensitivity[1,2] = PrevsAC[2,2]/(PrevsAC[2,2]+PrevsAC[2,1]) *100
      SpecificitySensitivity[2,2] = PrevsAC[1,1]/(PrevsAC[1,1]+PrevsAC[1,2]) *100
colnames(SpecificitySensitivity) = c("Sensitivity", "Specifity")
heatmap=data.frame()
heatmap_rnames = common_genes
uy = "_AC"
uy1= "_SCC"
for(i in 1:71){
  for(j in 1:100){
  if(heatmap_rnames[i]==rownames(t1_sq_100sig)[j]){
    if(t1_sq_100sig[j,1]>0){
      heatmap_rnames[i] = paste(heatmap_rnames[i], uy)
    }
    else{
      heatmap_rnames[i] = paste(heatmap_rnames[i] , uy1)
    }
}
  }
}
for(i in 1:71){
  for(i in 1:100)
}
rownames(V_combined_71)
heatmap= V_combined_71
rownames(heatmap) = heatmap_rnames
rownames(heatmap)
for (i in 1:384){
  if(i<=192){
    colnames(heatmap)[i]= paste(colnames(heatmap)[i], "_SCC")
  }
  else{
    colnames(heatmap)[i]= paste(colnames(heatmap)[i], "_AC")
  }
}
heatmap
write.table(heatmap, file="V_71_heatmap.txt",sep="\t",col.names = NA,row.names = T)
heatmap_rev= data.frame(row.names = rownames(heatmap))
for(i in 384:1){
    heatmap_rev = cbind(heatmap_rev, heatmap[,i])
}
colnames(heatmap)[1]==colnames(heatmap_rev)[384]
colnames(heatmap_rev)
heatmap_rev
write.table(heatmap_rev, file="V_71_heatmap_unstd.txt",sep="\t",col.names = NA,row.names = T)
heatmap_std= data.frame(row.names = rownames(heatmap))
heatmap_col <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(heatmap_col)
heatmap.2(data.matrix(heatmap_rev),col=rev(morecols(50)),trace="none",scale="row")
heatmap[1,1]==heatmap_rev[1,384]
rstd_heatmap = heatmap_rev
stdev_heatmap = sd(data.matrix(rstd_heatmap))
stdev_heatmap
mean_heatmap = mean(data.matrix((rstd_heatmap)))
mean_heatmap
for(i in 1:71){
  for(j in 1:384){
    std_value = (rstd_heatmap[i,j] - mean_heatmap)/stdev_heatmap
    rstd_heatmap[i,j] = std_value
  }
}
morecols <- colorRampPalette(heatmap_col)
heatmap.2(data.matrix(rstd_heatmap),col=rev(morecols(50)),trace="none",scale="row")
##### by row standardization
rwstd_heatmap = heatmap_rev
for(i in 1:71){
  stdev_row = sd(data.matrix(rwstd_heatmap[i,]))
  mean_row = mean(data.matrix((rwstd_heatmap[i,])))
  for(j in 1:384){
    std_value = (rwstd_heatmap[i,j] - mean_row)/stdev_row
    rwstd_heatmap[i,j] = std_value
  }
}
morecols <- colorRampPalette(heatmap_col)
heatmap.2(data.matrix(rwstd_heatmap),col=rev(my_palette),trace="none",scale="row")
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
write.table(rwstd_heatmap,file="RowStadndardize_Validation_TCGA.txt",col.names = NA,row.names = T,sep = "\t")
rrwstd_heatmap = heatmap
for(i in 1:71){
  stdev_row = sd(data.matrix(rrwstd_heatmap[i,]))
  mean_row = mean(data.matrix((rrwstd_heatmap[i,])))
  for(j in 1:384){
    std_value = (rrwstd_heatmap[i,j] - mean_row)/stdev_row
    rrwstd_heatmap[i,j] = std_value
  }
}
#####################################################
Correlation_cluster= correlation[c(1,2), ]
heatmap(data.matrix(Correlation_cluster))
