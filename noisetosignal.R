setwd("/media/ahadli/Data/Farid/AOG_lab/COAD/Secil_microarray_020718")
x = read.csv("GSE55139-withnonpaired-logfc.csv",header = T,row.names=1, stringsAsFactors = F)
x_new = x[c(-1,-length(x))]
x_nnew = log2(x_new+0.00000001)
ediff = c()
for(i in 1:nrow(x_new)){
  mm = max(x_new[i,]) - min(x_new[i,])
  diff = c(diff, mm)
}
alll = sum(diff<1)
