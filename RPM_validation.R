# import the mirlist
setwd("D:/Farid/AOG lab/Alperin isleri/miRNA_analysis")
v_AD_RPM=data.frame(v_AD_mol_RPM, v_AD_n_mol_RPM)
colnames(v_AD_RPM) = c(rep(x="AD",times=92), rep(x="norm", times=18))

mir_cutoff_AD_mol=read.csv(file="farid_cutoff_miRNAs_AD_mol.csv", stringsAsFactors = F, check.names = F, row.names = 1)
colnames(mir_cutoff_AD_mol) = rownames(t1_AD__mol_ttest)[5:14]
mir_cutoff_AD_mol[1,10]=0
test_counts = data.frame()
n=0
for(i in 1: length(mir_cutoff_AD_mol)){
  for(j in 1: dim(v_AD_RPM)[1]){
    if(colnames(mir_cutoff_AD_mol)[i]==rownames(v_AD_RPM)[j]){
      test_counts = rbind(test_counts, v_AD_RPM[j,])
      n=n+1
    }
  }
}

cutoffs = data.frame()
counts_normal = test_counts[,(length(v_AD_mol)+1):length(test_counts)]
for(i in 1: dim(test_counts)[1]){
  single_cutoff = max(counts_normal[i,])
  cutoffs = rbind(cutoffs, single_cutoff)
}
rownames(cutoffs)= rownames(test_counts)
logic_validation = data.frame(row.names = rownames(test_counts))
for(i in 1:dim(test_counts)[2]){
  logic = test_counts[,i] > cutoffs
  logic_validation = cbind(logic_validation, logic)
}
colnames(logic_validation) = colnames(test_counts) 
logic_validation = logic_validation[,(colSums(logic_validation)==0)]
colnames(logic_validation)

# import the mirlist
setwd("D:/Farid/AOG lab/Alperin isleri/miRNA_analysis")
v_SQ_RPM=data.frame(v_SQ_mol_RPM, v_SQ_n_mol_RPM)
colnames(v_SQ_RPM) = c(rep(x="SQ",times=60), rep(x="norm", times=23))

mir_cutoff_SQ_mol=read.csv(file="cutoff_miRNAs_SQ_mol.csv", stringsAsFactors = F, check.names = F, row.names = 1)
test_counts = data.frame()
n=0
for(i in 1: length(mir_cutoff_SQ_mol)){
  for(j in 1: dim(v_SQ_RPM)[1]){
    if(colnames(mir_cutoff_SQ_mol)[i]==rownames(v_SQ_RPM)[j]){
      test_counts = rbind(test_counts, v_SQ_RPM[j,])
      n=n+1
    }
  }
}

cutoffs = data.frame()
counts_normal = test_counts[,(length(v_SQ_mol)+1):length(test_counts)]
for(i in 1: dim(test_counts)[1]){
  single_cutoff = max(counts_normal[i,])
  cutoffs = rbind(cutoffs, single_cutoff)
}
rownames(cutoffs)= rownames(test_counts)
logic_validation = data.frame(row.names = rownames(test_counts))
for(i in 1:dim(test_counts)[2]){
  logic = test_counts[,i] > cutoffs
  logic_validation = cbind(logic_validation, logic)
}
colnames(logic_validation) = colnames(test_counts) 
logic_validation = logic_validation[,(colSums(logic_validation)==0)]
colnames(logic_validation)

