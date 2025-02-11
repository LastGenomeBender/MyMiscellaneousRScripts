#### adeno molecular DGE without prefiltering using DESeq2 and min LogFC of 2
setwd("D:/Farid/AOG lab/Alperin isleri/miRNA_analysis/log2FC2_DESeq2_prefiltered")
t1_AD_primary_cnames= rep("ADmol", times=length(t1_AD_mol))
t1_AD_normal_cnames =  rep("ADnorm", times=length(t1_AD_n_mol))
t1.ADprivsADnorm = data.frame(t1_AD_mol,t1_AD_n_mol)
colnames(t1.ADprivsADnorm)= c(t1_AD_primary_cnames,t1_AD_normal_cnames)

library("DESeq2")
oricln = colnames(t1.ADprivsADnorm)
cln = make.names(oricln, unique = T)
colnames(t1.ADprivsADnorm)= cln
coldata = data.frame(row.names = cln)
coldata[,1] = as.factor(oricln)
colnames(coldata) = "condition"
ddsHTSeq <- DESeqDataSetFromMatrix( countData = t1.ADprivsADnorm, colData= coldata,
                                    design= ~ condition)
ddsHTSeq
dds <- DESeq(ddsHTSeq)
log2norm_AD_mol  = counts(dds,normalized=T)
log2norm_AD_mol = log2(log2norm_AD_mol+1)
log2norm_AD_mol = data.frame(log2norm_AD_mol)
res <- results(dds, contrast = c("condition","ADmol","ADnorm"))
res
res = res[!is.na(res$padj),]
res = res[order(res$pvalue),]
res2= res[res$log2FoldChange>2,]
res2[,7]= res2$pvalue * dim(res)[1]
colnames(res2)[7] = "Benferonni"
res2
write.csv(res2, file = "T1_NormalvsTumor_AD_Molecular_DESeq2.csv")
write.csv(log2norm_AD_mol, file = "T1_AD_mol_normalizedlog2_count.csv")

#### squamous molecular DGE without prefiltering using DESeq2 and min LogFC of 2

t1_SQ_primary_cnames = rep("SQmol", times=length(t1_SQ_mol))
t1_SQ_normal_cnames =  rep("SQnorm", times=length(t1_SQ_n_mol))
t1.SQprivsSQnorm = data.frame(t1_SQ_mol,t1_SQ_n_mol)
colnames(t1.SQprivsSQnorm)= c(t1_SQ_primary_cnames,t1_SQ_normal_cnames)
library("DESeq2")
oricln = colnames(t1.SQprivsSQnorm)
cln = make.names(oricln, unique = T)

colnames(t1.SQprivsSQnorm)= cln
coldata = data.frame(row.names = cln)
coldata[,1] = as.factor(oricln)
colnames(coldata) = "condition"
ddsHTSeq <- DESeqDataSetFromMatrix( countData = t1.SQprivsSQnorm, colData= coldata,
                                    design= ~ condition)
ddsHTSeq
dds <- DESeq(ddsHTSeq)
log2norm_SQ_mol = counts(dds, normalized = T)
log2norm_SQ_mol = log2(log2norm_SQ_mol+1)
log2norm_SQ_mol = data.frame(log2norm_SQ_mol)
dds
res <- results(dds,contrast = c("condition","SQmol","SQnorm"))
res
res = res[!is.na(res$padj),]
res = res[order(res$pvalue),]
res2= res[res$log2FoldChange>2,]
res2[,7]= res2$pvalue * dim(res)[1]
colnames(res2)[7] = "Benferonni"
res2
write.csv(res2, file = "T1_NormalvsTumor_SQ_Molecular_DESeq2.csv")
write.csv(log2norm_SQ_mol, file = "T1_SQ_mol_normalizedlog2_count.csv")
