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
library(gage)
library("pasilla")
library("DESeq2")
library("tibble")
library(gageData)
library(BiocManager)
install("gageData")
library(gageData)
arr_rpkm = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
arr_readcount = read.delim("/home/ahadli/Desktop/CCLE_RNAseq_genes_counts_20180929.gct",header = T, sep="\t",row.names = 1, stringsAsFactors = F, check.names = F)
my_cells_as_ccle = as.vector(as.character(arr_rpkm[2,c(2:ncol(arr_rpkm))]))
gnames = arr_readcount$Description
arr_readcount_selected = arr_readcount[my_cells_as_ccle]
counts = arr_readcount_selected[c(1:nrow(arr_readcount_selected)),c(2:ncol(arr_readcount_selected))]
colnames(counts) = colnames(AHR_edit_rppa_selected)
counts = counts[c(AHR_edit,AHR_nonedit)]
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]

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
data("egSymb")
####### Edits are the Reference
ref.idx=1:23
samp.idx=24:33
kegg.gs = kegg.sets.hs
aaa=sym2eg(gnames)
aaa_bool = !is.na(aaa)
cnts.norm = cnts.norm[aaa_bool,]
rownames(cnts.norm) = aaa[aaa_bool]
sum(is.na(aaa))
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
setwd(dir = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/01102019_RNA_seq_GAGE_readcounts/more in nonedited AHR")
pv.out.list <- sapply(path.ids1, function(pid) pathview(gene.data = cnts.d, pathway.id = pid,species = "hsa"))
write.csv(cnts.kegg.p$greater,"greater.csv")
setwd(dir = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/01102019_RNA_seq_GAGE_readcounts/less in nonedited AHR")
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = cnts.d, pathway.id = pid,species = "hsa"))
write.csv(cnts.kegg.p$less,"less.csv")
getwd()

cnts.kegg.p
str(cnts.kegg.p, strict.width='wrap')
str(cnts.kegg.p)

write.csv(cnts.kegg.p$greater,"greater.csv")