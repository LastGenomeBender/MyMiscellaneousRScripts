library("limma")
library("edgeR")
library(gplots)
library("org.Mm.eg.db")
library("RColorBrewer")
library("Glimma")
#set the working directory with setwd()

seqdata <- read.delim("LUSCHiSeqV2", row.names = "sample")
x=sample(1:3,553, replace = T)
seqdata[20531,] = x
n1=1
n2=1
n3=1
cnames = colnames(seqdata)
s1 = c("")
s2 = c("")
s3 = c("")
t1_squamous_rep2 = data.frame(row.names = rownames(seqdata)) 
t2_squamous_rep2 = data.frame(row.names = rownames(seqdata))
v_squamous_rep2 = data.frame(row.names = rownames(seqdata))
for(i in 1:553){
  determ=seqdata[20531,i]
  if(determ==1){
    t1_squamous_rep2[,n1] = seqdata[,i] 
    n1= n1+1
    s1 = c(s1, cnames[i])
  } 
  else if(determ==2){
    t2_squamous_rep2[,n2] = seqdata[,i]
    n2=n2+1
    s2 = c(s2, cnames[i])
  } 
  else{
    v_squamous_rep2[,n3]= seqdata[,i]
    n3=n3+1
    s3= c(s3, cnames[i])
  }
}
s1 = s1[2:length(s1)]
s2 = s2[2:length(s2)]
s3 = s3[2:length(s3)]
colnames(t1_squamous_rep2) = s1
colnames(t2_squamous_rep2) = s2
colnames(v_squamous_rep2) =s3
t1_squamous_rep2=t1_squamous_rep2[-20531,]
t2_squamous_rep2=t2_squamous_rep2[-20531,]
v_squamous_rep2=v_squamous_rep2[-20531,]
write.csv(t1_squamous_rep2, file="T1_Squamous_TCGA_rep2.csv")
write.csv(t2_squamous_rep2,file="T2_Squamous_TCGA_rep2.csv")
write.csv(v_squamous_rep2,file="V_Squamous_TCGA_rep2.csv")
