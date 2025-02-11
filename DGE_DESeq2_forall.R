source("https://bioconductor.org/biocLite.R")
biocLite("pasilla")
biocLite("DESeq2")
install.packages("tibble")
library(tibble)
t1_AD_primary_cnames= rep("ADpri", times=length(t1_AD_p_mi))
t1_AD_normal_cnames =  rep("ADnorm", times=length(t1_AD_n_mi))
t1.ADprivsADnorm = data.frame(t1_AD_p_mi,t1_AD_n_mi)
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
res <- results(dds)
res
res = res[!is.na(res$pvalue),]
res = res[order(res$pvalue),]
res2= res[1:100,]
res2[,7]= res2$pvalue * dim(res)[1]
colnames(res2)[7] = "Benferonni"
write.csv(res2, file = "T1_NormalvsTumor_AD_DESeq2.csv")
colnames(t1.ADprivsADnorm) = oricln
res_logfc= res[order(-res$log2FoldChange),]
res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
res_logfc2= res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
write.csv(res_logfc2, file = "T1_TumorvsNormal_AD_LogFCsorted_DESeq2.csv")

#### t2 adeno toumor vs normal
t2_AD_primary_cnames= rep("ADpri", times=length(t2_AD_p_mi))
t2_AD_normal_cnames =  rep("ADnorm", times=length(t2_AD_n_mi))
t2.ADprivsADnorm = data.frame(t2_AD_p_mi,t2_AD_n_mi)
colnames(t2.ADprivsADnorm)= c(t2_AD_primary_cnames,t2_AD_normal_cnames)
library("DESeq2")
oricln = colnames(t2.ADprivsADnorm)
cln = make.names(oricln, unique = T)
colnames(t2.ADprivsADnorm)= cln
coldata = data.frame(row.names = cln)
coldata[,1] = as.factor(oricln)
colnames(coldata) = "condition"
ddsHTSeq <- DESeqDataSetFromMatrix( countData = t2.ADprivsADnorm, colData= coldata,
                                    design= ~ condition)
ddsHTSeq
dds <- DESeq(ddsHTSeq)
res <- results(dds)
res
res = res[!is.na(res$pvalue),]
res = res[order(res$pvalue),]
res2= res[1:100,]
res2[,7]= res2$pvalue * dim(res)[1]
colnames(res2)[7] = "Benferonni"
write.csv(res2, file = "t2_TumorvsNormal_AD_DESeq2.csv")
colnames(t2.ADprivsADnorm) = oricln
res_logfc= res[order(-res$log2FoldChange),]
res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
res_logfc2= res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
write.csv(res_logfc2, file = "t2_TumorvsNormal_AD_LogFCsorted_DESeq2.csv")

###### T1 squamous
t1_SQ_primary_cnames= rep("SQpri", times=length(t1_SQ_p_mi))
t1_SQ_normal_cnames =  rep("SQnorm", times=length(t1_SQ_n_mi))
t1.SQprivsSQnorm = data.frame(t1_SQ_p_mi,t1_SQ_n_mi)
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
res <- results(dds)
res
res = res[!is.na(res$pvalue),]
res = res[order(res$pvalue),]
res2= res[1:100,]
res2[,7]= res2$pvalue * dim(res)[1]
colnames(res2)[7] = "Benferonni"
write.csv(res2, file = "T1_NormalvsTumor_SQ_DESeq2.csv")
colnames(t1.SQprivsSQnorm) = oricln
res_logfc= res[order(-res$log2FoldChange),]
res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
res_logfc2= res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
write.csv(res_logfc2, file = "T1_TumorvsNormal_SQ_LogFCsorted_DESeq2.csv")

###### T2 squamous
t2_SQ_primary_cnames= rep("SQpri", times=length(t2_SQ_p_mi))
t2_SQ_normal_cnames =  rep("SQnorm", times=length(t2_SQ_n_mi))
t2.SQprivsSQnorm = data.frame(t2_SQ_p_mi,t2_SQ_n_mi)
colnames(t2.SQprivsSQnorm)= c(t2_SQ_primary_cnames,t2_SQ_normal_cnames)
library("DESeq2")
oricln = colnames(t2.SQprivsSQnorm)
cln = make.names(oricln, unique = T)
colnames(t2.SQprivsSQnorm)= cln
coldata = data.frame(row.names = cln)
coldata[,1] = as.factor(oricln)
colnames(coldata) = "condition"
ddsHTSeq <- DESeqDataSetFromMatrix( countData = t2.SQprivsSQnorm, colData= coldata,
                                    design= ~ condition)
ddsHTSeq
dds <- DESeq(ddsHTSeq)
res <- results(dds)
res
res = res[!is.na(res$pvalue),]
res = res[order(res$pvalue),]
res2= res[1:100,]
res2[,7]= res2$pvalue * dim(res)[1]
colnames(res2)[7] = "Benferonni"
write.csv(res2, file = "t2_NormalvsTumor_SQ_DESeq2.csv")
colnames(t2.SQprivsSQnorm) = oricln
res_logfc= res[order(-res$log2FoldChange),]
res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
res_logfc2= res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
write.csv(res_logfc2, file = "t2_TumorvsNormal_SQ_LogFCsorted_DESeq2.csv")

