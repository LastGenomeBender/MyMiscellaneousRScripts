##### read the RNA editing file
editings = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/significant_editing.csv", header = T, row.names = 1, check.names = F, stringsAsFactors = F)
rpkm_data = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv",header = T, row.names = 1, check.names = F, stringsAsFactors = F) 

#### there are 33 cells with the expression value while 36 cells with editings. lets intersect them
genes = strsplit(x = rownames(editings),split = "|",fixed = T)
genes2 = unlist(lapply(genes, function(X) X[5])) 
new_mat = data.frame()
rownames(rpkm_data) = rpkm_data$V1
genenames = rpkm_data$V2
genenames = genenames[c(-1,-2)]
my_cell_names = as.vector(as.character(rpkm_data[1,]))
my_cell_names = my_cell_names[c(-1,-2)]
ccle_cell_names = as.vector(as.character(rpkm_data[2,]))
ccle_cell_names = ccle_cell_names[c(-1,-2)]
rpkm_data = rpkm_data[c(-1,-2)]
rpkm_data = rpkm_data[c(-1,-2),]
names(genenames)  = rownames(rpkm_data)

for(i in 1:ncol(rpkm_data)){
  rpkm_data[,i] = as.numeric(rpkm_data[,i])
}
mean_RPKM_each_gene = data.frame()
for(i in 1:length(genes2)){
  gene = genes2[i]
  bool = genenames == gene
  temp = rpkm_data[bool,]
  if(dim(temp)[1]>1){
    ###take mean rpkm as it is
    out = colMeans(temp)
  }
  else{
    out = temp
  }
  mean_RPKM_each_gene = rbind(mean_RPKM_each_gene, out)
}
rownames(mean_RPKM_each_gene) = rownames(editings)
diff

  genes2

bool2 = editings == 0
sum(bool)
for( i in 1:nrow(editings))
