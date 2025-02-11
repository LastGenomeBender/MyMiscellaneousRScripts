################ expresion difference between edited and nonedited AHR
## read the rpkm files of that bad 33
arr = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
meta = arr[c(1,2)]
arr = arr[c(-1,-2),c(-1,-2)]
arr = as.data.frame(sapply(arr, as.numeric))
arr = log2(arr+1)
edit = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
edit_selected = edit[,(colnames(arr)[1:ncol(arr)])]
AHR_edit_rppa_selected = edit_selected[grepl(pattern ="7|17384354",x =  rownames(edit_selected),fixed = T),]
  AHR_nonedit = colnames(AHR_edit_rppa_selected[,AHR_edit_rppa_selected==0])
  AHR_edit = setdiff(colnames(AHR_edit_rppa_selected),y = AHR_nonedit)
sum(meta$V2=="ULBP2")

arr_edit_tot = arr[c(AHR_edit, AHR_nonedit)]
arr_edit_
############## compare PAI-1
arr_pai=arr_edit_tot[meta$V2=="ULBP2",]
arr_edit_ulbp =  t(arr_pai[AHR_edit])
arr_nonedit_ulbp =  t(arr_pai[AHR_nonedit])
boxplot(arr_edit_ulbp)
boxplot(arr_nonedit_ulbp)

  meta
arr_edit_tot_rnms = meta$V2[c(-1,-2)]
arr_edit_tot_rnms = arr_edit_tot_rnms[!duplicated(arr_edit_tot_rnms)]
arr_edit_tot_unique = data.frame()
for(i in 1:length(arr_edit_tot_rnms)){
  temp = arr_edit_tot[meta$V2 == arr_edit_tot_rnms[i],]
   if(dim(temp)[1]>1) {temp = colMeans((temp))} else{} 
  arr_edit_tot_unique = rbind(arr_edit_tot_unique, temp) 
}

ensgs = vapply()

rownames(arr_edit_tot_unique) = arr_edit_tot_rnms
write.csv(arr_edit_tot_unique)
arr_edit = arr[AHR_edit]
write.csv(AHR)
setwd("/home/ahadli/Desktop")

write.table(arr_edit_tot,"AHR_edir_vs_nonedit_rpkm.txt",sep = "\t" )
arr_nonedit = arr[AHR_nonedit]
vec.ttest = c()
log_fc_vec = c()
for( i in 1:nrow(arr_edit)){
  a=t.test(x = as.vector(as.numeric(arr_nonedit[i,])),y=as.vector(as.numeric(arr_edit[i,])))$p.value
  logFC_nonedit_edit = mean(as.vector(as.numeric(arr_nonedit[i,])))- mean(as.vector(as.numeric(arr_edit[i,])))
  vec.ttest = c(vec.ttest,a )
  log_fc_vec = c(log_fc_vec,logFC_nonedit_edit)
}

lvec.ttest_adjust = p.adjust(vec.ttest,method = "BH")
bool = vec.ttest<0.05
which(bool==T)
bool_corrected = vec.ttest_adjust<0.05
which(bool_corrected==T)

arr_nonedit = arr_nonedit[bool,]
arr_edit = arr_edit[bool,]
log_fc_vec = log_fc_vec[bool]
vec.ttest = vec.ttest[bool]
meta = meta[c(-1,-2),]
meta = meta[bool,]
meta 
final = cbind(arr_nonedit,arr_edit,log_fc_vec,vec.ttest)
final = final[order(final$vec.ttest,decreasing = F),]
library(BiocManager)
install("gage")
library(gage)

final_gage=final[c(-35,-34)]
get.gset
aaa=gage(exprs = final_gage,gsets = ,samp = 1:10,compare = ,ref = 11:23)

kegg.gsets(species = hsa,id.type = )


##### try deseq2

library(BiocManager)
library(gage)
install("pasilla")
install("DESeq2")
install.packages("tibble")
arr_rpkm = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
arr_readcount = read.delim("/home/ahadli/Desktop/CCLE_RNAseq_genes_counts_20180929.gct",header = T, sep="\t",row.names = 1, stringsAsFactors = F, check.names = F)
my_cells_as_ccle = as.vector(as.character(arr_rpkm[2,c(2:ncol(arr_rpkm))]))
gnames = arr_readcount$Description
arr_readcount_selected = arr_readcount[my_cells_as_ccle]
counts = arr_readcount_selected[c(1:nrow(arr_readcount_selected)),c(2:ncol(arr_readcount_selected))]
colnames(counts) = colnames(AHR_edit_rppa_selected)
counts = counts[c(AHR_edit,AHR_nonedit)]
data(kegg.gs)

class(counts[1,1])
cnts=counts
sel.rn=rowSums(cnts) != 0
gnames=gnames[sel.rn]
cnts=cnts[sel.rn,]
##joint workflow with DEseq/edgeR/limma/Cufflinks forks here
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
cnts.norm=log2(cnts.norm+8)
# t2_SQ_primary_cnames= rep("SQpri", times=length(t2_SQ_p_mi))
# t2_SQ_normal_cnames =  rep("SQnorm", times=length(t2_SQ_n_mi))
# t2.SQprivsSQnorm = data.frame(t2_SQ_p_mi,t2_SQ_n_mi)
# colnames(t2.SQprivsSQnorm)= c(t2_SQ_primary_cnames,t2_SQ_normal_cnames)
# library("DESeq2")
# oricln = colnames(t2.SQprivsSQnorm)
# cln = make.names(oricln, unique = T)
# colnames(t2.SQprivsSQnorm)= cln
# coldata = data.frame(row.names = cln)
# coldata[,1] = as.factor(oricln)
# colnames(coldata) = "condition"
# ddsHTSeq <- DESeqDataSetFromMatrix( countData = t2.SQprivsSQnorm, colData= coldata,
#                                     design= ~ condition)
# ddsHTSeq
# dds <- DESeq(ddsHTSeq)
# res <- results(dds)
# res
# res = res[!is.na(res$pvalue),]
# res = res[order(res$pvalue),]
# res2= res[1:100,]
# res2[,7]= res2$pvalue * dim(res)[1]
# colnames(res2)[7] = "Benferonni"
# write.csv(res2, file = "t2_NormalvsTumor_SQ_DESeq2.csv")
# colnames(t2.SQprivsSQnorm) = oricln
# res_logfc= res[order(-res$log2FoldChange),]
# res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
# res_logfc2= res_logfc[c(1:50, (length(rownames(res_logfc))-49):length(rownames(res_logfc))),]
 ##step 3: gage
   ##joint workflow with DEseq/edgeR/limma/Cufflinks merges around here
   library(gage)
data("egSymb")
 ref.idx=1:23
 samp.idx=24:33
 kegg.gs = kegg.gsets("hsa")
 aaa=sym2eg(gnames)
 aaa_bool = !is.na(aaa)
 cnts.norm = cnts.norm[aaa_bool,]
 rownames(cnts.norm) = aaa[aaa_bool]
 sum(is.na(aaa))
 data("kegg.gs")
 cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = ref.idx,
                       samp = samp.idx, compare ="unpaired")
 cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])
  sel <- cnts.kegg.p$greater[, "q.val"] < 0.25 &
   !is.na(cnts.kegg.p$greater[,"q.val"])
  path.ids <- rownames(cnts.kegg.p$greater)[sel]
  sel.l <- cnts.kegg.p$less[, "q.val"] < 0.25 &
    !is.na(cnts.kegg.p$less[,"q.val"])
  path.ids.l <- rownames(cnts.kegg.p$less)[sel.l]
  path.ids2 <- substr(c(path.ids.l), 1, 8)
  path.ids1 = substr(c(path.ids), 1, 8)
  install("pathview")
  library(pathview)
  setwd(dir = "/home/ahadli/Desktop/Less_in_nonedit")
  pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = cnts.d, pathway.id = pid,species = "hsa"))
  setwd("/home/ahadli/Desktop/greater_in_nonedit")
  pv.out.list <- sapply(path.ids1, function(pid) pathview(gene.data = cnts.d, pathway.id = pid,species = "hsa"))
  
    getwd()

  cnts.kegg.p
str(cnts.kegg.p, strict.width='wrap')
str(cnts.kegg.p)

write.csv(cnts.kegg.p$greater,"greater.csv")
