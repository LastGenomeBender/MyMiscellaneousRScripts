setwd("/media/ahadli/Data/Farid/AOG_lab/COAD/Secil_microarray_020718")
x = read.csv("Untitled 1.csv",header = T,row.names=1)
anno=read.delim("Annot.adf.csv",stringsAsFactors = F,header = T,row.names = 1)
sel_count = x[rownames(sel_stat),]
anno_selected= anno[rownames(sel_stat),]
bool= duplicated(anno_selected[,1])
sel_count = sel_count[!bool,]
anno_selected=anno_selected[!bool,]
rownames(sel_count)=as.vector(t(anno_selected[1]))
sel_count = sel_count[rownames(sel_count)!="1",]
samplewise_diff = data.frame(row.names = rownames(sig))
for(i in seq(1,20,2)){
  diff = sel_count[i]- sel_count[i+1]
  samplewise_diff = cbind(samplewise_diff,diff)
}
write.csv(samplewise_diff,"preop-postop_samplewise.csv")

mwise = matrix()
for (i in 1:nrow(samplewise_diff)){
  mwise = c(mwise, sum(samplewise_diff[i,]>0))
}
swise = matrix()
for ( i in 1:ncol(samplewise_diff)){
  swise = c(swise, sum(samplewise_diff[,i]>0))
}
swise = swise[-1]
mwise = mwise[-1]
samplewise_diff = cbind(samplewise_diff, mwise)
samplewise_diff = rbind(samplewise_diff, swise)
write.csv(samplewise_diff,"preop-postop_samplewise.csv")
