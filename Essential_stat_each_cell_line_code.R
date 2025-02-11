all_count_etypes = list()
all_exon_intron_dist_each_etype = list()
chrwise_etype_dist = list()
hedited_genes = list()
hedited_genes_AG = list()
hedited_genes_TC = list()
for(i in 1:length(all_files)) {
  setwd(paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/filtered files/",names(all_files)[i]))
  mat = all_files[[i]]
  ##### count of each editing type
  tot_etype = all_files[[i]]$RNAE_t
  e_types = levels(factor(tot_etype))
  count_etypes = c()
  for(j in 1:length(e_types)) {
    count_etypes = c(count_etypes, sum(e_types[j]==tot_etype))
  }
  type_dist = data.frame(row.names = "1")
  count_etypes=rbind(type_dist, count_etypes)
  colnames(count_etypes) = e_types
  all_count_etypes[[names(all_files)[i]]] = count_etypes
  write.csv(count_etypes,paste0("Edit_type_distribution_",names(all_files)[i]))
  # ###### END
  #
  # ##### Exon, int, ve.s dist of each editing type
  editing_location_by_etype = data.frame()
  seq_ontology = c("synonymous_variant","missense_variant","stop_codon_gain","3_prime_UTR","5_prime_UTR","intron_variant")
  for(j in 1:length(e_types)){
    bool_temp = tot_etype == e_types[j]
    mat_temp = mat[bool_temp,]
    temp_count = c()
    for(k in 1:length(seq_ontology)){
      temp_count = c(temp_count, sum(grep(seq_ontology[k],mat_temp$func_anno)>0))
    }
    editing_location_by_etype = rbind(editing_location_by_etype, temp_count)
  }
  editing_location_by_etype = as.data.frame(editing_location_by_etype)
  colnames(editing_location_by_etype) = c("Synonymous","Missense","Nonsense","3'UTR","5'UTR","Intron")
  rownames(editing_location_by_etype) = e_types
  all_exon_intron_dist_each_etype[[names(all_files)[i]]] = editing_location_by_etype
  write.csv(editing_location_by_etype,paste0("Exon_vs_distribution_of_each_Edit_Type",names(all_files)[i]))
  # ######## END
  #
  # ####### Chromosome-wise editing distribution
  chr = mat$chr
  chr_lev =levels(factor(chr))
  chrwise_editingType_dist = data.frame()
  for(j in 1:length(chr_lev)){
    bool_temp = chr == chr_lev[j]
    mat_temp = mat[bool_temp,]
    temp_count = c()
    for(k in 1:length(e_types)){
      temp_count = c(temp_count, sum(mat_temp$RNAE_t == e_types[k]))
    }
    chrwise_editingType_dist = rbind(chrwise_editingType_dist,temp_count)
  }
  chrwise_editingType_dist = as.data.frame(chrwise_editingType_dist)
  rownames(chrwise_editingType_dist) = chr_lev
  colnames(chrwise_editingType_dist) = e_types
  chrwise_etype_dist[[names(all_files)[i]]] = chrwise_editingType_dist
  write.csv(chrwise_editingType_dist, paste0("Chromosomewise_distribution_of_each_Edit_Type",names(all_files)[i]))
  # ####### END
  
  ####### Genes with highest number of edited sites
  tot_genes = mat$gene
  genes = levels(factor(tot_genes))
  tot_count =c()
  for(j in 1:length(genes)){
    bool_temp = tot_genes == genes[j]
    tot_count = c(tot_count, sum(bool_temp))
  }
  max_genes = data.frame(row.names = genes)
  max_genes = cbind(max_genes,tot_count)
  order_idx = order(max_genes$tot_count, decreasing = T)
  rnmss = rownames(max_genes)[order_idx]
  max_genes = max_genes[order_idx,]
  nnn = data.frame(row.names = rnmss)
  nnn=cbind(nnn,max_genes)
  hedited_genes[[names(all_files)[i]]] = nnn
  write.csv(nnn, paste0("Max_edited_genes",names(all_files)[i]))
  ################################ AG and TC max editings
  for(k in 1:2){
    ifelse(test = k == 1, yes = etyp<- "AG", no = etyp<- "TC")
    booltp = mat$RNAE_t == etyp
    mat_Etype = mat[booltp,]
    tot_genes = mat_Etype$gene 
    genes = levels(factor(tot_genes))
    tot_count =c()
    for(j in 1:length(genes)){
      bool_temp = tot_genes == genes[j]
      tot_count = c(tot_count, sum(bool_temp))
    }
    if(length(genes)==0){
      nnn = NA
    }
    else {
      max_genes = data.frame(row.names = genes)
      max_genes = cbind(max_genes,tot_count)
      order_idx = order(max_genes$tot_count, decreasing = T)
      rnmss = rownames(max_genes)[order_idx]
      max_genes = max_genes[order_idx,]
      nnn = data.frame(row.names = rnmss)
      nnn=cbind(nnn,max_genes)
    }
    ifelse(test=k==1,yes=hedited_genes_AG[[names(all_files)[i]]] <- nnn, no<-hedited_genes_TC[[names(all_files)[i]]] <- nnn)
    
    write.csv(nnn, paste0(etyp,"_edit_Max_edited_genes",names(all_files)[i]))
  }
}
save(all_count_etypes,
     all_exon_intron_dist_each_etype,
     chrwise_etype_dist,
     hedited_genes, all_files, file,hedited_genes_AG, file  = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/all_lists.RData")
