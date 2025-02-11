load("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/all_lists.RData")
library(ggplot2)
######## histogram of editing types
e_types = colnames(all_count_etypes[[2]])
overall_edit = data.frame()
for(i in 1:length(all_count_etypes)){
  mat = as.data.frame(t(all_count_etypes[[i]]))
  mat = as.data.frame(t(mat[e_types,]))
  overall_edit = rbind(overall_edit, mat)
}
colnames(overall_edit) = e_types
overall_edit[is.na(overall_edit)] = 0
a = getwd()
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE")
pdf(file = "overall_general_stats_snpless_SNVs.pdf",paper = "a4")
setwd(a)
boxplot(overall_edit, outline = F, xlab= "SNV Type", ylab="# of SNVs")
stripchart(overall_edit, vertical = T, pch=20,method = "jitter",add = T,col ="red")
title("SNP-less SNVs")

########################### overall distribution to the AG sites
seq_type_dist = data.frame()
for(i in 1:length(all_exon_intron_dist_each_etype)){
  mat = all_exon_intron_dist_each_etype[[i]]
  mat = mat[rownames(mat)=="AG",]
  seq_type_dist = rbind(seq_type_dist,mat)
}
colnames(seq_type_dist) = colnames(all_exon_intron_dist_each_etype[[1]])
boxplot(seq_type_dist, xlab= "Sequence Type", ylab = "# of Editings")
title(main = "A-to-G SNV")
stripchart(seq_type_dist, vertical = T, pch=20,method = "jitter",add = T,col ="red")
explanat_AG_seqtype_mat = data.frame(row.names = colnames(seq_type_dist))
AG_avg_each_seqtype = colSums(seq_type_dist)/nrow(seq_type_dist)
AG_max_each_seqtype = apply(seq_type_dist, 2,max)
AG_min_each_seqtype = apply(seq_type_dist, 2,min)
explanat_AG_seqtype_mat=cbind(explanat_AG_seqtype_mat,AG_min_each_seqtype,AG_max_each_seqtype,AG_avg_each_seqtype)
colnames(explanat_AG_seqtype_mat) = c("Minimum", "Maximum", "Average")
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE")
write.csv(explanat_AG_seqtype_mat, "AtoG_sequence_distribution_stats.csv")
########## END


########################### overall distribution to the TC sites
seq_type_dist = data.frame()
for(i in 1:length(all_exon_intron_dist_each_etype)){
  mat = all_exon_intron_dist_each_etype[[i]]
  mat = mat[rownames(mat)=="TC",]
  seq_type_dist = rbind(seq_type_dist,mat)
}
colnames(seq_type_dist) = colnames(all_exon_intron_dist_each_etype[[1]])
boxplot(seq_type_dist, xlab= "Sequence Type", ylab = "# of Editings")
title(main = "T-to-C Editings")
stripchart(seq_type_dist, vertical = T, pch=20,method = "jitter",add = T,col ="red")
explanat_AG_seqtype_mat = data.frame(row.names = colnames(seq_type_dist))
AG_avg_each_seqtype = colSums(seq_type_dist)/nrow(seq_type_dist)
AG_max_each_seqtype = apply(seq_type_dist, 2,max)
AG_min_each_seqtype = apply(seq_type_dist, 2,min)
explanat_AG_seqtype_mat=cbind(explanat_AG_seqtype_mat,AG_min_each_seqtype,AG_max_each_seqtype,AG_avg_each_seqtype)
colnames(explanat_AG_seqtype_mat) = c("Minimum", "Maximum", "Average")
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI")
write.csv(explanat_AG_seqtype_mat, "TtoC_sequence_distribution_stats.csv")

########################## distribution of seq type for each Etype
seq_type_dist2 = all_exon_intron_dist_each_etype[[2]]
for(i in c(1,3:length(all_exon_intron_dist_each_etype))){
  mat = all_exon_intron_dist_each_etype[[i]]
  mat2 = mat[rownames(seq_type_dist2),]
  mat2[is.na(mat2)] = 0
  seq_type_dist2 = seq_type_dist2 + mat2
}
####################
etype_seq_dist = list()
for (i in 1:length(e_types)){
  mat = data.frame()
  for(j in 1:length(all_exon_intron_dist_each_etype)){
    mat2 = all_exon_intron_dist_each_etype[[j]]
    mat2 = mat2[e_types[i],]
    mat = rbind(mat, mat2)
  }
  mat[is.na(mat)] = 0
  colnames(mat) = colnames(all_exon_intron_dist_each_etype[[j]])
  rownames(mat) = names(all_exon_intron_dist_each_etype)
  etype_seq_dist[[e_types[i]]] = mat
}
for(i in 1:length(etype_seq_dist)){
  mat = etype_seq_dist[[i]]
  boxplot(mat, xlab= "Sequence Type", ylab = "# of Editings", outline = F,title = paste0(rownames(seq_type_dist2)[i]," Editing"))
  stripchart(mat, vertical = T, pch=20,method = "jitter",add = T,col ="red")
  title(paste(names(etype_seq_dist),"Editing Sequence Type Distribution"))
}

############# Chromosomewise editings
chr_names = c(seq(1:22), "X","Y")
edit_each_chr = list()
for (i in 1:length(chr_names)){
  mat = data.frame()
  for(j in 1:length(chrwise_etype_dist)){
    mat2 = chrwise_etype_dist[[j]]
    mat2 = as.data.frame(t(mat2[chr_names[i],]))
    mat2 = as.data.frame(t(mat2[e_types,]))
    mat = rbind(mat, mat2)
  }
  mat[is.na(mat)] = 0
  colnames(mat) = e_types
  rownames(mat) = names(chrwise_etype_dist)
  edit_each_chr[[chr_names[i]]] = mat
}

for(i in 1:length(edit_each_chr)){
  mat = edit_each_chr[[i]]
  boxplot(mat, xlab= "Editing Type", ylab = "# of Editings", outline = F)
  stripchart(mat, vertical = T, pch=20,method = "jitter",add = T,col ="red") 
  title(paste0("Chromosome ", names(edit_each_chr)[i]))
}
#########chrwise editing A-to-G
chrwise_dist = data.frame(row.names = chr_names)
for(i in 1:length(chrwise_etype_dist)){
  mat = as.data.frame(chrwise_etype_dist[[i]])
  mat = mat[chr_names,]
  mat[is.na(mat)] = 0
  mat = mat
  chrwise_dist = cbind((chrwise_dist),(mat))
}
chrwise_dist = chrwise_dist$AG
View(chrwise_dist)
colnames(chrwise_dist) = names(chrwise_etype_dist)

#####################END
chrwise_dist = data.frame(row.names = chr_names)
for(i in 1:length(chrwise_etype_dist)){
  mat = as.data.frame(chrwise_etype_dist[[i]])
  mat = mat[chr_names,]
  mat[is.na(mat)] = 0
  mat = mat
  ifelse(test = is.null(mat), yes = mat <- rep(0,times=24),no = mat<-mat)
  chrwise_dist = cbind((chrwise_dist),(mat))
}
View(chrwise_dist)
colnames(chrwise_dist) = names(chrwise_etype_dist)


boxplot(t(chrwise_dist), xlab= "Chromosome", ylab = "# of SNVs",outline = F)
stripchart(data.frame(t(chrwise_dist)), vertical = T, pch=20,method = "jitter",add = T,col ="red")
title("Chromosome Distribution of A-to-G SNVs")
write.csv(chrwise_dist,"A-to-G_chrwise_dist.csv")
###################  Gene wise decision and ranksum test
unique_genes = c()
for(i in 1:length(hedited_genes)){
  mat = hedited_genes[[i]]
  rank = rank(mat$max_genes, ties.method = "max")
  rank = nrow(mat) + 1 - rank
  gns = rownames(mat)
  unique_genes = c(unique_genes, gns)
  mat$rank = rank
  hedited_genes[[i]] = mat
}

unique_genes = unique_genes[!duplicated(unique_genes)]

all_ranks = data.frame(row.names = unique_genes)
for(i in 1:length(hedited_genes)){
  mat = hedited_genes[[i]]
  mat = mat[unique_genes,]
  mat[is.na(mat)] = 480
  all_ranks = cbind(all_ranks, mat$rank)
}
colnames(all_ranks) = names(hedited_genes)
all_ranks$rankSum = rowSums(all_ranks)
all_ranks=all_ranks[order(all_ranks$rankSum),]
all_ranks$rankSum


rownames(all_ranks)
max(all_ranks[1,c(1:(ncol(all_ranks)-1))])
all_ranks[1,c(1:(ncol(all_ranks)-1))]
all_ranks$G30645.SW48.3_final_funcanno_usethis.csv  
all_ranks2 = all_ranks[-ncol(all_ranks)]
all_ranks2 = all_ranks2[!(colnames(all_ranks2)=="G30645.SW48.3_final_funcanno_usethis.csv")]
max_20 =as.data.frame(t(all_ranks2[1:20,]))
boxplot(max_20)
stripchart(max_20, vertical = T, pch=20,method = "jitter",add = T,col ="red")
max_50 = as.data.frame(t(all_ranks2[1:50,]))
boxplot(max_50)
stripchart(max_50, vertical = T, pch=20,method = "jitter",add = T,col ="red")

colnames(max_20)
cat(paste0(), sep = "\n")
all_ranks[1:20, ncol(all_ranks)] 
write.csv(all_ranks, "ranks_genes_number_of_editing sites.csv")

############### editing num of max 200
max200nms = rownames(all_ranks)[1:200]
edit_max200 = data.frame(row.names = max200nms) 
for(i in 1:length(hedited_genes)){
  mat = hedited_genes[[i]]
  mat = mat[max200nms,]$max_genes
  mat[is.na(mat)]=0
  edit_max200=cbind(edit_max200, mat)
}
colnames(edit_max200) = names(hedited_genes)
frplot = as.data.frame(t(edit_max200))
frplot20 = frplot[1:20]
boxplot(frplot[-1], outline = F)
stripchart(frplot[-1], vertical = T, pch=20,method = "jitter",add = T,col ="red")

boxplot(frplot20[-1], outline = F)
stripchart(frplot20[-1], vertical = T, pch=20,method = "jitter",add = T,col ="red")
##################### mean editing site ranking ranking
editing_genes_all = data.frame()
for(i in 1:length(hedited_genes)){
  mat = hedited_genes[[i]]
  mat = mat[unique_genes,]$max_genes
  mat[is.na(mat)]=0
  editing_genes_all=rbind(editing_genes_all, mat)
}
editing_genes_all = as.data.frame(t(editing_genes_all))
rownames(editing_genes_all) = unique_genes
colnames(editing_genes_all) = names(hedited_genes)
editing_genes_all$Total_edit_site = rowSums(editing_genes_all)
editing_genes_all2 = editing_genes_all[order(editing_genes_all$Total_edit_site, decreasing = T), ]
plot(editing_genes_all2$Total_edit_site[-1], ylab = "Sum of the all edited site in all cell lines", xlab = "Rank")
write.csv(editing_genes_all2,"sum_edited_sites_eachGene_ranked.csv")
boxplot(as.data.frame(t(editing_genes_all2[-1, -ncol(editing_genes_all2)]))[1:200],title = "Top 200 genes (Tot Edit Sites)")
boxplot(as.data.frame(t(editing_genes_all2[-1, -ncol(editing_genes_all2)]))[1:20],title = "Top 200 genes (Tot Edit Sites)",outline = F)
stripchart(as.data.frame(t(editing_genes_all2[-1, -ncol(editing_genes_all2)]))[1:20], vertical = T, pch=20,method = "jitter",add = T,col ="red")


########### AtoG mean Edit site Ranking
editing_genes_all = data.frame()
for(i in 1:length(hedited_genes_AG)){
  if(is.null(hedited_genes_AG[[i]])){
    mat = rep(0, times =  ncol(editing_genes_all))
  }
  else{
    mat = hedited_genes_AG[[i]]
    mat = mat[match(unique_genes,rownames(mat)),]
  }
  mat[is.na(mat)]=0
  editing_genes_all=rbind(editing_genes_all, mat)
}
editing_genes_all = as.data.frame(t(editing_genes_all))
rownames(editing_genes_all) = unique_genes
colnames(editing_genes_all) = names(hedited_genes)
editing_genes_all$Mean_edit_site = rowSums(editing_genes_all)/ncol(editing_genes_all)
editing_genes_all2 = editing_genes_all[order(editing_genes_all$Mean_edit_site, decreasing = T), ]
plot(editing_genes_all2$Mean_edit_site[-1], ylab = " Mean # of the AG SNV sites in all cell lines", xlab = "Rank")
title("Mean A-to-G SNV #s of Genes")
qeqes = editing_genes_all2$Mean_edit_site[-1]
hist(qeqes,breaks = 12, plot = T,right = F,xlab = "Mean # of SNV Sites", ylab = "# of genes",main = "Histogram of Mean # of SNV sites (AG)")
write.csv(editing_genes_all2,"AG_mean_SNV_sites_eachGene_ranked.csv")
boxplot(as.data.frame(t(editing_genes_all2[-1, -ncol(editing_genes_all2)]))[1:200],title = "Top 200 genes (AG Edit Sites)")
title("Top 200 genes (AG SNV Sites)")
boxplot(as.data.frame(t(editing_genes_all2[-1, -ncol(editing_genes_all2)]))[1:20],title = "Top 20 genes (AG Edit Sites)",varwidth = T,outline = F,las=2)
title("Top 20 genes (AG SNV Sites)")
stripchart(as.data.frame(t(editing_genes_all2[-1, -ncol(editing_genes_all2)]))[1:20], vertical = T, pch=20,method = "jitter",add = T,col ="red")
max(editing_genes_all2$Mean_edit_site[-1])



################TC mean edit sites all genes
editing_genes_all = data.frame()
for(i in 1:length(hedited_genes_TC)){
  if(is.na(hedited_genes_TC[[i]])){
    mat = rep(0, times = ncol(editing_genes_all))
  }
  else{
    mat = hedited_genes_TC[[i]]
    mat = mat[unique_genes,]
  }
  mat[is.na(mat)]=0
  editing_genes_all=rbind(editing_genes_all, mat)
}
editing_genes_all = as.data.frame(t(editing_genes_all))
rownames(editing_genes_all) = unique_genes
colnames(editing_genes_all) = names(hedited_genes)
editing_genes_all$Mean_edit_site = rowSums(editing_genes_all)/ncol(editing_genes_all)
editing_genes_all2 = editing_genes_all[order(editing_genes_all$Mean_edit_site, decreasing = T), ]
plot(editing_genes_all2$Mean_edit_site[-1], ylab = "Mean of the all edited site in all cell lines", xlab = "Rank")
title("Mean T-to-C Editing of Genes")
qeqes = editing_genes_all2$Mean_edit_site[-1]
hist(qeqes,breaks = 3, plot = F,right = F,xlab = "Mean Edited Sites", ylab = "# of genes",main = "Histogram of Mean Edited Sites (TC Editing)")
write.csv(editing_genes_all2,"TC_mean_edited_sites_eachGene_ranked.csv")
boxplot(as.data.frame(t(editing_genes_all2[-1, -ncol(editing_genes_all2)]))[1:200],title = "Top 200 genes (TC Edit Sites)")
title("Top 200 genes (TC Edit Sites)")
boxplot(as.data.frame(t(editing_genes_all2[-1, -ncol(editing_genes_all2)]))[1:20],title = "Top 20 genes (TC Edit Sites)",varwidth = T,outline = F)
title("Top 20 genes (TC Edit Sites)")
stripchart(as.data.frame(t(editing_genes_all2[-1, -ncol(editing_genes_all2)]))[1:200], vertical = T, pch=20,method = "jitter",add = T,col ="red")

############# Top 20 genes and their editing sites
edit_mat = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
edit_mat_AG = edit_mat[grepl("\\|AG\\|", rownames(edit_mat)),]
edit_mat_TC = edit_mat[grepl("\\|TC\\|", rownames(edit_mat)),]
top20_AG = read.csv("AG_mean_SNV_sites_eachGene_ranked.csv", row.names = 1,header = T, stringsAsFactors = F, check.names = F)
top20_AG = top20_AG[-1,]
top20_AG = top20_AG[1:(ncol(top20_AG)-7)]
top20_TC = read.csv("TC_mean_edited_sites_eachGene_ranked.csv", row.names = 1,header = T, stringsAsFactors = F, check.names = F)
top20_TC = top20_TC[-1,]
top20_TC = top20_TC[1:(ncol(top20_TC)-6)]
top20_exp_eachsite_AG = data.frame()
for(i in 1:20){
  geneofint = paste0("|",rownames(top20_AG)[i],"||||")
  edit_int = edit_mat_AG[grepl(geneofint,rownames(edit_mat_AG),fixed = T),]
  edit_int$mean_edit = rowMeans(edit_int)
  final_mean = c(edit_int$mean_edit, rep(0,times=1000-nrow(edit_int))) 
  top20_exp_eachsite_AG = rbind(top20_exp_eachsite_AG,final_mean)
}
rownames(top20_exp_eachsite_AG) = rownames(top20_AG)[1:20]
write.csv(top20_exp_eachsite_AG,"AG_top20genes_expression_of_eachsite.csv")

geneofint = paste0("|GNL3||||")
edit_int = edit_mat_AG[grepl(geneofint,rownames(edit_mat_AG),fixed = T),]
edit_int = edit_int>0
sum(edit_int)

########
top20_exp_eachsite_TC = data.frame()
for(i in 1:20){
  geneofint = paste0("\\|",rownames(top20_TC)[i],"\\|")
  edit_int = edit_mat_TC[grepl(geneofint,rownames(edit_mat_TC)),]
  edit_int$mean_edit = rowMeans(edit_int)
  final_mean = c(edit_int$mean_edit, rep(0,times=100-nrow(edit_int))) 
  top20_exp_eachsite_TC = rbind(top20_exp_eachsite_TC,final_mean)
}
rownames(top20_exp_eachsite_TC) = rownames(top20_TC)[1:20]
write.csv(top20_exp_eachsite_TC,"TC_top20genes_expression_of_eachsite.csv")

###################### editing level of cell lines.

edit_mat = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
edit_mat = edit_mat[grepl("\\|AG\\|", rownames(edit_mat)),] #"\\|AG\\|" for AG
edit_mat_bool = edit_mat>0
edit_mat_bool2 = colSums(edit_mat_bool)
edit_mat_bool2
mean(edit_mat_bool2)
min(edit_mat_bool2)
max(edit_mat_bool2)
sum(edit_mat_bool2)
class(edit_mat_bool2)
cell_lines_nms = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/cell_line_names.csv",row.names = 1, header = T)
cllll_nms = as.character(cell_lines_nms$x)
names(edit_mat_bool2) = cell_lines_nms
barplot(edit_mat_bool2, las=2)
title("A-to-G Editings")
dev.off()
