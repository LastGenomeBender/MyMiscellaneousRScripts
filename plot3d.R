##### my miRNAs
mymirs = "hsa-miR-34b-3p|hsa-miR-34c-3p"
mymirs = unlist(strsplit(mymirs, split = "|", fixed = T))
arr = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
names_as_in_ccle = as.vector(as.character(arr[2, c(3:ncol(arr))]))
####### read the mir mat  
mir_mat = read.delim2("/home/ahadli/Desktop/CCLE_miRNA_20181103.gct",sep = "\t")
sum(duplicated(mir_mat$Description))
rownames(mir_mat) = mir_mat$Description
mir_mat = mir_mat[c(-1,-2)]
mir_mat = mir_mat[,names_as_in_ccle]
mir_mat2 = mir_mat[mymirs,]
mir_mat2[1,] = mir_mat3[1,]
rownames(mir_mat2)[1] = "hsa-miR-34b"
mir_mat3 = mir_mat[grepl(x=rownames(mir_mat), pattern = "hsa-miR-34b",fixed = T),]
grepl(x=rownames(mir_mat), pattern = "hsa-miR-34b",fixed = T)
sum(grepl(x=rownames(mir_mat), pattern = "hsa-miR-34b",fixed = T))


indices = order(as.vector(as.numeric(AHR_edit_rppa_selected[1,])),decreasing = F)

AHR_edit_rppa_selected_ordered = AHR_edit_rppa_selected[indices]

rppa_selected_ordered = rppa_selected[indices,]
mir_mat2 = mir_mat2[indices]

plot(as.vector(as.numeric(AHR_edit_rppa_selected_ordered[1,])),mir_mat2[2,])
install.packages("plot3D")
install.packages("plot3Drgl")
library(plot3D)
library(plot3Drgl)
plotrgl(lighting = T)
scatter3D(x=as.numeric(AHR_edit_rppa_selected_ordered[1,]), y=as.vector(as.numeric(mir_mat2[2,])),z=rppa_selected_ordered$Snail_Caution,xlab = "edit",ylab="miRNA",zlab="snail",plot = T)
scatter2D(x=rppa_selected_ordered$Snail_Caution, y=as.vector(as.numeric(mir_mat2[2,])))
