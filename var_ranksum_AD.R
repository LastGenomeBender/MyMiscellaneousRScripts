var_AD_mol = t1_AD__mol_ttest[1:(length(t1_AD__mol_ttest)-2)]
var_tumor = data.frame()
var_normal = data.frame()
for(i in 1:dim(var_AD_mol)[1]){
  var_tumor =  rbind(var_tumor,apply(t1_AD__mol_ttest[i,1:81],1,var))
  var_normal = rbind(var_normal,apply(t1_AD__mol_ttest[i,82:109],1,var))
  }
rownames(var_tumor)=rownames(var_AD_mol)
rownames(var_normal)=rownames(var_AD_mol)
var_tumor[,2] = 0
var_normal[,2]= 0
var_tumor = var_tumor[order(var_tumor$X1.09504445462298),]
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
write.csv(summed_rank, file="miRs from lowest var ranksum_AD.csv")
