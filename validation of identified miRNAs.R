#### validation of Adeno miRNAs
v_AD_mol
v_AD_n_mol
v_AD_mol_cnames = rep("ADmol", times = length(v_AD_mol)) 
v_AD_n_mol_cnames  = rep("normal", times = length(v_AD_n_mol))
validation_AD_mol_cnames = c(v_AD_mol_cnames, v_AD_n_mol_cnames)
unicln = make.names(validation_AD_mol_cnames, unique = T)
validation_AD_molecular = data.frame(row.names = rownames(v_AD_mol))
validation_AD_molecular = cbind(v_AD_mol, v_AD_n_mol) 
boolean = validation_AD_molecular==0
boolean = rowSums(boolean)
boolean = boolean<=18
validation_AD_molecular= validation_AD_molecular[boolean,]
colnames(validation_AD_molecular) = unicln
coldata = data.frame(row.names = unicln)
coldata[,1] = as.factor(validation_AD_mol_cnames)
colnames(coldata) = "condition"
ddsHTSeq <- DESeqDataSetFromMatrix( countData = validation_AD_molecular, colData= coldata,
                                    design= ~ condition)
ddsHTSeq
dds <- DESeq(ddsHTSeq)
log2norm_AD_mol_validation  = counts(dds,normalized=T)
log2norm_AD_mol_validation = log2(log2norm_AD_mol_validation+1)
log2norm_AD_mol_validation = data.frame(log2norm_AD_mol_validation)
write.csv(log2norm_AD_mol_validation, file = "V_normalized_AD_molecular.csv")

# import the mirlist

mir_cutoff_AD_mol=read.csv(file="cutoff_miRNAs_AD_mol.csv", stringsAsFactors = F, check.names = F, row.names = 1)

test_counts = data.frame()
n=0
for(i in 1: length(mir_cutoff_AD_mol)){
  for(j in 1: dim(log2norm_AD_mol_validation)[1]){
    if(colnames(mir_cutoff_AD_mol)[i]==rownames(log2norm_AD_mol_validation)[j]){
      test_counts = rbind(test_counts, log2norm_AD_mol_validation[j,])
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


########################### validation of SQ miRNAs
v_SQ_mol
v_SQ_n_mol
v_SQ_mol_cnames = rep("SQmol", times = length(v_SQ_mol)) 
v_SQ_n_mol_cnames  = rep("normal", times = length(v_SQ_n_mol))
validation_SQ_mol_cnames = c(v_SQ_mol_cnames, v_SQ_n_mol_cnames)
unicln = make.names(validation_SQ_mol_cnames, unique = T)
validation_SQ_molecular = data.frame(row.names = rownames(v_SQ_mol))
validation_SQ_molecular = cbind(v_SQ_mol, v_SQ_n_mol) 
boolean = validation_SQ_molecular==0
boolean = rowSums(boolean)
boolean = boolean<=length(v_SQ_n_mol)
validation_SQ_molecular= validation_SQ_molecular[boolean,]
colnames(validation_SQ_molecular) = unicln
coldata = data.frame(row.names = unicln)
coldata[,1] = as.factor(validation_SQ_mol_cnames)
colnames(coldata) = "condition"
ddsHTSeq <- DESeqDataSetFromMatrix( countData = validation_SQ_molecular, colData= coldata,
                                    design= ~ condition)
ddsHTSeq
dds <- DESeq(ddsHTSeq)
log2norm_SQ_mol_validation  = counts(dds,normalized=T)
log2norm_SQ_mol_validation = log2(log2norm_SQ_mol_validation+1)
log2norm_SQ_mol_validation = data.frame(log2norm_SQ_mol_validation)
write.csv(log2norm_SQ_mol_validation, file = "V_normalized_SQ_molecular.csv")

# import the mirlist

mir_cutoff_SQ_mol=read.csv(file="cutoff_miRNAs_SQ_mol.csv", stringsAsFactors = F, check.names = F, row.names = 1)

test_counts = data.frame()
n=0
for(i in 1: length(mir_cutoff_SQ_mol)){
  for(j in 1: dim(log2norm_SQ_mol_validation)[1]){
    if(colnames(mir_cutoff_SQ_mol)[i]==rownames(log2norm_SQ_mol_validation)[j]){
      test_counts = rbind(test_counts, log2norm_SQ_mol_validation[j,])
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
