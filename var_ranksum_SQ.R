var_SQ_mol = t1_SQ__mol_ttest[1:(length(t1_SQ__mol_ttest)-2)]
var_tumor = data.frame()
var_normal = data.frame()
for(i in 1:dim(var_SQ_mol)[1]){
  var_tumor =  rbind(var_tumor,apply(t1_SQ__mol_ttest[i,1:60],1,var))
  var_normal = rbind(var_normal,apply(t1_SQ__mol_ttest[i,61:23],1,var))
}
rownames(var_tumor)=rownames(var_SQ_mol)
rownames(var_normal)=rownames(var_SQ_mol)
var_tumor[,2] = 0
var_normal[,2]= 0
var_tumor = var_tumor[order(var_tumor[1]),]
var_normal = var_normal[order(var_normal[1]),]
summed_rank = data.frame()
for(i in 1: dim(var_tumor)[1]){
  for (j in 1:dim(var_normal)[1]){
    if(rownames(var_tumor)[i]==rownames(var_normal)[j]){
      summed_rank[i,1] = i+j 
    }
  }
}
rownames(summed_rank)= rownames(var_tumor)
summed_rank[,2] = 0
summed_rank = summed_rank[order(summed_rank[1]),]
summed_rank = t(summed_rank)
summed_rank = summed_rank[1,]
summed_rank
write.csv(summed_rank, file="miRs from lowest var ranksum_SQ.csv")
