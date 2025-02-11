#siginficant editings versus expressions.
# look at the zero editings
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_significant_editings_expression")
l = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot/logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="logplotsVsCor.pdf"]
bool11 = grepl(x = nopath,pattern = "AG and TC Variants Only",fixed = T)
nopath = nopath[bool11]
l = l[bool11]
c= 1
RPKM = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
rownames(RPKM) = RPKM$V1
genenames = RPKM$V2
genenames = genenames[c(-1,-2)]
my_cell_names = as.vector(as.character(RPKM[1,]))
my_cell_names = my_cell_names[c(-1,-2)]
ccle_cell_names = as.vector(as.character(RPKM[2,]))
ccle_cell_names = ccle_cell_names[c(-1,-2)]
RPKM = RPKM[c(-1,-2)]
RPKM = RPKM[c(-1,-2),]
names(genenames)  = rownames(RPKM)
for(i in 1:ncol(RPKM)){
  RPKM[,i] = as.numeric(RPKM[,i])
}
RPKM = log2(RPKM+1)
#c=16,
for (i in l){
  #i= l[1]
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_significant_editings_expression/",substr(justname, 1, nchar(justname)- 13), ".pdf")
  za =c("All Variants", "AG Variants", "AG and TC Variants Only")
  vec = c(grepl(x = nopath[c],pattern = "All Variants"), grepl(x = nopath[c],pattern ="AG Variants"),grepl(x = nopath[c],pattern = "AG and TC Variants Only"))
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = za[vec], fixed = T ))-2)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations/", drugname,"_correlation.csv")
  
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  drug_mat = drug_mat[,!is.na(drug_mat[2,])]
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
  drug_mat = drug_mat[-1,]
  rnms = rownames(drug_mat)
  drug_mat=as.data.frame(apply(drug_mat, 2, FUN =  as.numeric),row.names = rnms)
  cor_mat = as.data.frame(apply(cor_mat, 2, FUN =  as.numeric),row.names = rnms[-1])
  mat = read.csv(file=i, header = T, stringsAsFactors = F, row.names = 1, check.names = F)
  bb = is.na(mat$`Pearson R`)
  cc  =is.na(mat$`Spearman R`)
  mat = mat[!(bb|cc),]
  mat2 = mat[5:8]
  mat3 = colSums(mat2<0.05)
  
  if(nrow(mat2)==0){
    c= c+1
    next()
  }
  #dev.off()
  if(sum(mat3)==0){
    c = c+1
    next()
    
  }
  else{
    pdf(pdf_name,paper = "a4r")
    mat2  = mat2[rowSums(mat2<0.05)>0,]
    sel_rnms = rownames(mat2)
    
    temp_drug_mat = drug_mat[sel_rnms,]
    idx_sel_gnms = gregexpr(text = sel_rnms, pattern = "|" ,fixed = T )
    abs_idx_genes = list()
    c2 =1
    for(kj in 1:length(idx_sel_gnms)){
      idx_temp = idx_sel_gnms[[kj]]
      idx_temp = idx_temp[c(4,5)]
      idx_temp[1] = idx_temp[1]+1
      idx_temp[2] = idx_temp[2]-1
      sel_rnms[kj] = substr(sel_rnms[kj], idx_temp[1], idx_temp[2])
      at = which(toupper(genenames) == toupper(sel_rnms[kj]))
      ay = which(grepl(x = toupper(genenames), pattern = toupper(sel_rnms[kj]),fixed = T))
      gene_not_found = c()
      if((length(at)==0)& (length(ay)==0)){
        gene_not_found = c(gene_not_found, rownames(temp_drug_mat[kj]))
        abs_idx_genes[[c2]] = NA
      }
      else if (length(at)==0){
        abs_idx_genes[[c2]] = ay
      }
      else {
        abs_idx_genes[[c2]] =  at
        
      }
      c2 = c2+1
    }
    temp_rpkm = RPKM[unlist(abs_idx_genes)[!is.na(unlist(abs_idx_genes))],]
    cell_nms = intersect(colnames(temp_rpkm),colnames(drug_mat))
    temp_rpkm = temp_rpkm[cell_nms]
    if(length(cell_nms)< length(colnames(drug_mat))){
      notpresent = colnames(drug_mat)[!(cell_nms)]
      rownm = rownames(temp_rpkm)
      for(jk in 1:length(notpresent)){
        temp_rpkm = cbind(temp_rpkm,NA)
      }
      colnames(temp_rpkm)[(ncol(temp_rpkm)-length(notpresent)+1):(ncol(temp_rpkm))] = notpresent 
    }
    rownames(temp_rpkm)
    output_frame = rbind(temp_rpkm, temp_drug_mat)
    cor_pear_vec =c()
    cor_spear_vec =c()
    p_pear_vec =c()
    p_spear_vec =c()
    p_pear_BH_vec =c()
    p_pear_Be_vec =c()
    p_spear_BH_vec =c()
    p_spear_Be_vec =c()
    
    for(j in 1:nrow(temp_drug_mat)){
      edit_perctg = temp_drug_mat[j,]
      edit_perctg = edit_perctg[order(edit_perctg,decreasing = F)]
      expr_indices = abs_idx_genes[[j]]
      if(length(expr_indices)>1){
        expr_rpkm = RPKM[expr_indices,]
        expr_rpkm = colMeans(expr_rpkm)
        expr_rpkm = expr_rpkm[colnames(edit_perctg)]
        cor_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$estimate
        cor_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$estimate
        p_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$p.value
        p_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$p.value
        p_pear_BH = p.adjust(p = p_pear, method = "BH" )
        p_pear_Be = p.adjust(p=p_pear, method = "bonferroni")
        p_spear_BH = p.adjust(p = p_spear, method = "BH" )
        p_spear_Be = p.adjust(p=p_pear, method = "bonferroni")
        cor_pear_vec =c(cor_pear_vec,rep(cor_pear,length(expr_indices)))
        cor_spear_vec =c(cor_spear_vec,rep(cor_spear,length(expr_indices)))
        p_pear_vec =c(p_pear_vec,rep(p_pear,length(expr_indices)))
        p_spear_vec =c(p_spear_vec,rep(p_spear,length(expr_indices)))
        p_pear_BH_vec =c(p_pear_BH_vec,rep(p_pear_BH,length(expr_indices)))
        p_pear_Be_vec =c(p_pear_Be_vec,rep(p_pear_Be,length(expr_indices)))
        p_spear_BH_vec =c(p_spear_BH_vec,rep(p_spear_BH,length(expr_indices)))
        p_spear_Be_vec =c(p_spear_Be_vec,rep(p_spear_BH,length(expr_indices)))
        x = as.vector(as.numeric(edit_perctg))
        y = as.vector(as.numeric(expr_rpkm))
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = "Editing Percentage", ylab = paste0("log2(RPKM) ", sel_rnms[j]," Gene"),main =paste0(rownames(temp_drug_mat)[j], "\n", "Pearson R = ", cor_pear, " BH_Pearson = " ,p_pear_BH, "Bonf_Pearson = ",p_pear_Be,"\n","Spearman Rho = ", cor_spear, " BH_Spearman = " ,p_spear_BH, "Bonf_Spearman = ",p_spear_Be),cex.main = 0.65) 
        lines(pr.fit~x)
        next()
      }
      else if(is.na(expr_indices)){
        cor_pear_vec =c(cor_pear_vec,"gene name not found in expression array")
        cor_spear_vec =c(cor_spear_vec,"gene name not found in expression array")
        p_pear_vec =c(p_pear_vec,"gene name not found in expression array")
        p_spear_vec =c(p_spear_vec,"gene name not found in expression array")
        p_pear_BH_vec =c(p_pear_BH_vec,"gene name not found in expression array")
        p_pear_Be_vec =c(p_pear_Be_vec,"gene name not found in expression array")
        p_spear_BH_vec =c(p_spear_BH_vec,"gene name not found in expression array")
        p_spear_Be_vec =c(p_spear_Be_vec,"gene name not found in expression array")
        next()
        
      }
      else {
        expr_rpkm = RPKM[expr_indices,]
      }
      expr_rpkm = expr_rpkm[colnames(edit_perctg)]
      cor_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$estimate
      cor_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$estimate
      p_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$p.value
      p_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$p.value
      p_pear_BH = p.adjust(p = p_pear, method = "BH" )
      p_pear_Be = p.adjust(p=p_pear, method = "bonferroni")
      p_spear_BH = p.adjust(p = p_spear, method = "BH" )
      p_spear_Be = p.adjust(p=p_pear, method = "bonferroni")
      cor_pear_vec =c(cor_pear_vec,cor_pear)
      cor_spear_vec =c(cor_spear_vec,cor_spear)
      p_pear_vec =c(p_pear_vec,p_pear)
      p_spear_vec =c(p_spear_vec,p_spear)
      p_pear_BH_vec =c(p_pear_BH_vec,p_pear_BH)
      p_pear_Be_vec =c(p_pear_Be_vec,p_pear_Be)
      p_spear_BH_vec =c(p_spear_BH_vec,p_spear_BH )
      p_spear_Be_vec =c(p_spear_Be_vec,p_spear_Be)
      x = as.vector(as.numeric(edit_perctg))
      y = as.vector(as.numeric(expr_rpkm))
      fit = lm(y~x)
      pr.fit = predict(fit)
      plot(x,y, xlab = "Editing Percentage", ylab = paste0("log2(RPKM) ", sel_rnms[j]," Gene"),main =paste0(rownames(temp_drug_mat)[j], "\n", "Pearson R = ", cor_pear, " BH_Pearson = " ,p_pear_BH, "Bonf_Pearson = ",p_pear_Be,"\n","Spearman Rho = ", cor_spear, " BH_Spearman = " ,p_spear_BH, "Bonf_Spearman = ",p_spear_Be),cex.main = 0.65) 
      lines(pr.fit~x)
    }
    cor_pear_vec =c(cor_pear_vec,cor_pear_vec[!duplicated(cor_pear_vec)])
    cor_spear_vec =c(cor_spear_vec,cor_spear_vec[!duplicated(cor_spear_vec)])
    p_pear_vec =c(p_pear_vec,p_pear[!duplicated(p_pear_vec)])
    p_spear_vec =c(p_spear_vec,p_spear_vec[!duplicated(p_spear_vec)])
    p_pear_BH_vec =c(p_pear_BH_vec,p_pear_BH_vec[!duplicated(p_pear_BH_vec)])
    p_pear_Be_vec =c(p_pear_Be_vec,p_pear_Be_vec[!duplicated(p_pear_Be_vec)])
    p_spear_BH_vec =c(p_spear_BH_vec,p_spear_BH_vec[!duplicated(p_spear_BH_vec)])
    p_spear_Be_vec =c(p_spear_Be_vec,p_spear_Be_vec[!duplicated(p_spear_Be_vec)])
    rnms_output = rownames(output_frame) 
    output_frame = cbind(output_frame,cor_pear_vec ,
                         cor_spear_vec ,
                         p_pear_vec ,
                         p_spear_vec ,
                         p_pear_BH_vec ,
                         p_pear_Be_vec ,
                         p_spear_BH_vec ,
                         p_spear_Be_vec )
    fin_gnms = genenames[rownames(output_frame)]
    rownames(output_frame) =paste0(rownames(output_frame),"/",fin_gnms)
    write.csv(output_frame, paste0(drugname,"editing_vs_expression_log2RPKM.csv"))
    dev.off()
  }
 c= c+1
}

##### Introns
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_significant_editings_expression")
l = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlations_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlations_logplot/logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlations_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="logplotsVsCor.pdf"]
bool11 = grepl(x = nopath,pattern = "AG and TC Variants Only",fixed = T)
nopath = nopath[bool11]
l = l[bool11]
c= 1
RPKM = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
rownames(RPKM) = RPKM$V1
genenames = RPKM$V2
genenames = genenames[c(-1,-2)]
my_cell_names = as.vector(as.character(RPKM[1,]))
my_cell_names = my_cell_names[c(-1,-2)]
ccle_cell_names = as.vector(as.character(RPKM[2,]))
ccle_cell_names = ccle_cell_names[c(-1,-2)]
RPKM = RPKM[c(-1,-2)]
RPKM = RPKM[c(-1,-2),]
names(genenames)  = rownames(RPKM)
for(i in 1:ncol(RPKM)){
  RPKM[,i] = as.numeric(RPKM[,i])
}
RPKM = log2(RPKM+1)
#c=16,
for (i in l){
  #i= l[1]
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_significant_editings_expression/",substr(justname, 1, nchar(justname)- 13), ".pdf")
  za =c("All Variants", "AG Variants", "AG and TC Variants Only")
  vec = c(grepl(x = nopath[c],pattern = "All Variants"), grepl(x = nopath[c],pattern ="AG Variants"),grepl(x = nopath[c],pattern = "AG and TC Variants Only"))
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = za[vec], fixed = T ))-2)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations/", drugname,"_correlation.csv")
  
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  drug_mat = drug_mat[,!is.na(drug_mat[2,])]
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
  drug_mat = drug_mat[-1,]
  rnms = rownames(drug_mat)
  drug_mat=as.data.frame(apply(drug_mat, 2, FUN =  as.numeric),row.names = rnms)
  cor_mat = as.data.frame(apply(cor_mat, 2, FUN =  as.numeric),row.names = rnms[-1])
  mat = read.csv(file=i, header = T, stringsAsFactors = F, row.names = 1, check.names = F)
  bb = is.na(mat$`Pearson R`)
  cc  =is.na(mat$`Spearman R`)
  mat = mat[!(bb|cc),]
  mat2 = mat[5:8]
  mat3 = colSums(mat2<0.05)
  
  if(nrow(mat2)==0){
    c= c+1
    next()
  }
  #dev.off()
  if(sum(mat3)==0){
    c = c+1
    next()
    
  }
  else{
    pdf(pdf_name,paper = "a4r")
    mat2  = mat2[rowSums(mat2<0.05)>0,]
    sel_rnms = rownames(mat2)
    
    temp_drug_mat = drug_mat[sel_rnms,]
    idx_sel_gnms = gregexpr(text = sel_rnms, pattern = "|" ,fixed = T )
    abs_idx_genes = list()
    c2 =1
    for(kj in 1:length(idx_sel_gnms)){
      idx_temp = idx_sel_gnms[[kj]]
      idx_temp = idx_temp[c(4,5)]
      idx_temp[1] = idx_temp[1]+1
      idx_temp[2] = idx_temp[2]-1
      sel_rnms[kj] = substr(sel_rnms[kj], idx_temp[1], idx_temp[2])
      at = which(toupper(genenames) == toupper(sel_rnms[kj]))
      ay = which(grepl(x = toupper(genenames), pattern = toupper(sel_rnms[kj]),fixed = T))
      gene_not_found = c()
      if((length(at)==0)& (length(ay)==0)){
        gene_not_found = c(gene_not_found, rownames(temp_drug_mat[kj]))
        abs_idx_genes[[c2]] = NA
      }
      else if (length(at)==0){
        abs_idx_genes[[c2]] = ay
      }
      else {
        abs_idx_genes[[c2]] =  at
        
      }
      c2 = c2+1
    }
    temp_rpkm = RPKM[unlist(abs_idx_genes)[!is.na(unlist(abs_idx_genes))],]
    cell_nms = intersect(colnames(temp_rpkm),colnames(drug_mat))
    temp_rpkm = temp_rpkm[cell_nms]
    if(length(cell_nms)< length(colnames(drug_mat))){
      notpresent = colnames(drug_mat)[!(cell_nms)]
      rownm = rownames(temp_rpkm)
      for(jk in 1:length(notpresent)){
        temp_rpkm = cbind(temp_rpkm,NA)
      }
      colnames(temp_rpkm)[(ncol(temp_rpkm)-length(notpresent)+1):(ncol(temp_rpkm))] = notpresent 
    }
    rownames(temp_rpkm)
    output_frame = rbind(temp_rpkm, temp_drug_mat)
    cor_pear_vec =c()
    cor_spear_vec =c()
    p_pear_vec =c()
    p_spear_vec =c()
    p_pear_BH_vec =c()
    p_pear_Be_vec =c()
    p_spear_BH_vec =c()
    p_spear_Be_vec =c()
    
    for(j in 1:nrow(temp_drug_mat)){
      edit_perctg = temp_drug_mat[j,]
      edit_perctg = edit_perctg[order(edit_perctg,decreasing = F)]
      expr_indices = abs_idx_genes[[j]]
      if(length(expr_indices)>1){
        expr_rpkm = RPKM[expr_indices,]
        expr_rpkm = colMeans(expr_rpkm)
        expr_rpkm = expr_rpkm[colnames(edit_perctg)]
        cor_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$estimate
        cor_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$estimate
        p_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$p.value
        p_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$p.value
        p_pear_BH = p.adjust(p = p_pear, method = "BH" )
        p_pear_Be = p.adjust(p=p_pear, method = "bonferroni")
        p_spear_BH = p.adjust(p = p_spear, method = "BH" )
        p_spear_Be = p.adjust(p=p_pear, method = "bonferroni")
        cor_pear_vec =c(cor_pear_vec,rep(cor_pear,length(expr_indices)))
        cor_spear_vec =c(cor_spear_vec,rep(cor_spear,length(expr_indices)))
        p_pear_vec =c(p_pear_vec,rep(p_pear,length(expr_indices)))
        p_spear_vec =c(p_spear_vec,rep(p_spear,length(expr_indices)))
        p_pear_BH_vec =c(p_pear_BH_vec,rep(p_pear_BH,length(expr_indices)))
        p_pear_Be_vec =c(p_pear_Be_vec,rep(p_pear_Be,length(expr_indices)))
        p_spear_BH_vec =c(p_spear_BH_vec,rep(p_spear_BH,length(expr_indices)))
        p_spear_Be_vec =c(p_spear_Be_vec,rep(p_spear_BH,length(expr_indices)))
        x = as.vector(as.numeric(edit_perctg))
        y = as.vector(as.numeric(expr_rpkm))
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = "Editing Percentage", ylab = paste0("log2(RPKM) ", sel_rnms[j]," Gene"),main =paste0(rownames(temp_drug_mat)[j], "\n", "Pearson R = ", cor_pear, " BH_Pearson = " ,p_pear_BH, "Bonf_Pearson = ",p_pear_Be,"\n","Spearman Rho = ", cor_spear, " BH_Spearman = " ,p_spear_BH, "Bonf_Spearman = ",p_spear_Be),cex.main = 0.65) 
        lines(pr.fit~x)
        next()
      }
      else if(is.na(expr_indices)){
        cor_pear_vec =c(cor_pear_vec,"gene name not found in expression array")
        cor_spear_vec =c(cor_spear_vec,"gene name not found in expression array")
        p_pear_vec =c(p_pear_vec,"gene name not found in expression array")
        p_spear_vec =c(p_spear_vec,"gene name not found in expression array")
        p_pear_BH_vec =c(p_pear_BH_vec,"gene name not found in expression array")
        p_pear_Be_vec =c(p_pear_Be_vec,"gene name not found in expression array")
        p_spear_BH_vec =c(p_spear_BH_vec,"gene name not found in expression array")
        p_spear_Be_vec =c(p_spear_Be_vec,"gene name not found in expression array")
        next()
        
      }
      else {
        expr_rpkm = RPKM[expr_indices,]
      }
      expr_rpkm = expr_rpkm[colnames(edit_perctg)]
      cor_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$estimate
      cor_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$estimate
      p_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$p.value
      p_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$p.value
      p_pear_BH = p.adjust(p = p_pear, method = "BH" )
      p_pear_Be = p.adjust(p=p_pear, method = "bonferroni")
      p_spear_BH = p.adjust(p = p_spear, method = "BH" )
      p_spear_Be = p.adjust(p=p_pear, method = "bonferroni")
      cor_pear_vec =c(cor_pear_vec,cor_pear)
      cor_spear_vec =c(cor_spear_vec,cor_spear)
      p_pear_vec =c(p_pear_vec,p_pear)
      p_spear_vec =c(p_spear_vec,p_spear)
      p_pear_BH_vec =c(p_pear_BH_vec,p_pear_BH)
      p_pear_Be_vec =c(p_pear_Be_vec,p_pear_Be)
      p_spear_BH_vec =c(p_spear_BH_vec,p_spear_BH )
      p_spear_Be_vec =c(p_spear_Be_vec,p_spear_Be)
      x = as.vector(as.numeric(edit_perctg))
      y = as.vector(as.numeric(expr_rpkm))
      fit = lm(y~x)
      pr.fit = predict(fit)
      plot(x,y, xlab = "Editing Percentage", ylab = paste0("log2(RPKM) ", sel_rnms[j]," Gene"),main =paste0(rownames(temp_drug_mat)[j], "\n", "Pearson R = ", cor_pear, " BH_Pearson = " ,p_pear_BH, "Bonf_Pearson = ",p_pear_Be,"\n","Spearman Rho = ", cor_spear, " BH_Spearman = " ,p_spear_BH, "Bonf_Spearman = ",p_spear_Be),cex.main = 0.65) 
      lines(pr.fit~x)
    }
    cor_pear_vec =c(cor_pear_vec,cor_pear_vec[!duplicated(cor_pear_vec)])
    cor_spear_vec =c(cor_spear_vec,cor_spear_vec[!duplicated(cor_spear_vec)])
    p_pear_vec =c(p_pear_vec,p_pear[!duplicated(p_pear_vec)])
    p_spear_vec =c(p_spear_vec,p_spear_vec[!duplicated(p_spear_vec)])
    p_pear_BH_vec =c(p_pear_BH_vec,p_pear_BH_vec[!duplicated(p_pear_BH_vec)])
    p_pear_Be_vec =c(p_pear_Be_vec,p_pear_Be_vec[!duplicated(p_pear_Be_vec)])
    p_spear_BH_vec =c(p_spear_BH_vec,p_spear_BH_vec[!duplicated(p_spear_BH_vec)])
    p_spear_Be_vec =c(p_spear_Be_vec,p_spear_Be_vec[!duplicated(p_spear_Be_vec)])
    rnms_output = rownames(output_frame) 
    #tryCatch({
    #output_frame = cbind(output_frame,cor_pear_vec ,
     #                    cor_spear_vec ,
    #                     p_pear_vec ,
    #                     p_spear_vec ,
    #                     p_pear_BH_vec ,
    #                     p_pear_Be_vec ,
    #                     p_spear_BH_vec ,
    #                     p_spear_Be_vec )
    #fin_gnms = genenames[rownames(output_frame)]
    #rownames(output_frame) =paste0(rownames(output_frame),"/",fin_gnms)
    #write.csv(output_frame, paste0(drugname,"editing_vs_expression_log2RPKM.csv"))
    #}, error = function(e){})
    dev.off()
  }
  c= c+1
}


### coding
dir.create("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_significant_editings_expression")
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_significant_editings_expression")
l = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlations_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlations_logplot/logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlations_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="logplotsVsCor.pdf"]
bool11 = grepl(x = nopath,pattern = "AG and TC Variants Only",fixed = T)
nopath = nopath[bool11]
l = l[bool11]
c= 1
RPKM = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
rownames(RPKM) = RPKM$V1
genenames = RPKM$V2
genenames = genenames[c(-1,-2)]
my_cell_names = as.vector(as.character(RPKM[1,]))
my_cell_names = my_cell_names[c(-1,-2)]
ccle_cell_names = as.vector(as.character(RPKM[2,]))
ccle_cell_names = ccle_cell_names[c(-1,-2)]
RPKM = RPKM[c(-1,-2)]
RPKM = RPKM[c(-1,-2),]
names(genenames)  = rownames(RPKM)
for(i in 1:ncol(RPKM)){
  RPKM[,i] = as.numeric(RPKM[,i])
}
RPKM = log2(RPKM+1)
#c=16,
for (i in l){
  #i= l[1]
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_significant_editings_expression/",substr(justname, 1, nchar(justname)- 13), ".pdf")
  za =c("All Variants", "AG Variants", "AG and TC Variants Only")
  vec = c(grepl(x = nopath[c],pattern = "All Variants"), grepl(x = nopath[c],pattern ="AG Variants"),grepl(x = nopath[c],pattern = "AG and TC Variants Only"))
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = za[vec], fixed = T ))-2)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations/", drugname,"_correlation.csv")
  
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  drug_mat = drug_mat[,!is.na(drug_mat[2,])]
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
  drug_mat = drug_mat[-1,]
  rnms = rownames(drug_mat)
  drug_mat=as.data.frame(apply(drug_mat, 2, FUN =  as.numeric),row.names = rnms)
  cor_mat = as.data.frame(apply(cor_mat, 2, FUN =  as.numeric),row.names = rnms[-1])
  mat = read.csv(file=i, header = T, stringsAsFactors = F, row.names = 1, check.names = F)
  bb = is.na(mat$`Pearson R`)
  cc  =is.na(mat$`Spearman R`)
  mat = mat[!(bb|cc),]
  mat2 = mat[5:8]
  mat3 = colSums(mat2<0.05)
  
  if(nrow(mat2)==0){
    c= c+1
    next()
  }
  #dev.off()
  if(sum(mat3)==0){
    c = c+1
    next()
    
  }
  else{
    pdf(pdf_name,paper = "a4r")
    mat2  = mat2[rowSums(mat2<0.05)>0,]
    sel_rnms = rownames(mat2)
    
    temp_drug_mat = drug_mat[sel_rnms,]
    idx_sel_gnms = gregexpr(text = sel_rnms, pattern = "|" ,fixed = T )
    abs_idx_genes = list()
    c2 =1
    for(kj in 1:length(idx_sel_gnms)){
      idx_temp = idx_sel_gnms[[kj]]
      idx_temp = idx_temp[c(4,5)]
      idx_temp[1] = idx_temp[1]+1
      idx_temp[2] = idx_temp[2]-1
      sel_rnms[kj] = substr(sel_rnms[kj], idx_temp[1], idx_temp[2])
      at = which(toupper(genenames) == toupper(sel_rnms[kj]))
      ay = which(grepl(x = toupper(genenames), pattern = toupper(sel_rnms[kj]),fixed = T))
      gene_not_found = c()
      if((length(at)==0)& (length(ay)==0)){
        gene_not_found = c(gene_not_found, rownames(temp_drug_mat[kj]))
        abs_idx_genes[[c2]] = NA
      }
      else if (length(at)==0){
        abs_idx_genes[[c2]] = ay
      }
      else {
        abs_idx_genes[[c2]] =  at
        
      }
      c2 = c2+1
    }
    temp_rpkm = RPKM[unlist(abs_idx_genes)[!is.na(unlist(abs_idx_genes))],]
    cell_nms = intersect(colnames(temp_rpkm),colnames(drug_mat))
    temp_rpkm = temp_rpkm[cell_nms]
    if(length(cell_nms)< length(colnames(drug_mat))){
      notpresent = colnames(drug_mat)[!(cell_nms)]
      rownm = rownames(temp_rpkm)
      for(jk in 1:length(notpresent)){
        temp_rpkm = cbind(temp_rpkm,NA)
      }
      colnames(temp_rpkm)[(ncol(temp_rpkm)-length(notpresent)+1):(ncol(temp_rpkm))] = notpresent 
    }
    rownames(temp_rpkm)
    output_frame = rbind(temp_rpkm, temp_drug_mat)
    cor_pear_vec =c()
    cor_spear_vec =c()
    p_pear_vec =c()
    p_spear_vec =c()
    p_pear_BH_vec =c()
    p_pear_Be_vec =c()
    p_spear_BH_vec =c()
    p_spear_Be_vec =c()
    
    for(j in 1:nrow(temp_drug_mat)){
      edit_perctg = temp_drug_mat[j,]
      edit_perctg = edit_perctg[order(edit_perctg,decreasing = F)]
      expr_indices = abs_idx_genes[[j]]
      if(length(expr_indices)>1){
        expr_rpkm = RPKM[expr_indices,]
        expr_rpkm = colMeans(expr_rpkm)
        expr_rpkm = expr_rpkm[colnames(edit_perctg)]
        cor_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$estimate
        cor_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$estimate
        p_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$p.value
        p_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$p.value
        p_pear_BH = p.adjust(p = p_pear, method = "BH" )
        p_pear_Be = p.adjust(p=p_pear, method = "bonferroni")
        p_spear_BH = p.adjust(p = p_spear, method = "BH" )
        p_spear_Be = p.adjust(p=p_pear, method = "bonferroni")
        cor_pear_vec =c(cor_pear_vec,rep(cor_pear,length(expr_indices)))
        cor_spear_vec =c(cor_spear_vec,rep(cor_spear,length(expr_indices)))
        p_pear_vec =c(p_pear_vec,rep(p_pear,length(expr_indices)))
        p_spear_vec =c(p_spear_vec,rep(p_spear,length(expr_indices)))
        p_pear_BH_vec =c(p_pear_BH_vec,rep(p_pear_BH,length(expr_indices)))
        p_pear_Be_vec =c(p_pear_Be_vec,rep(p_pear_Be,length(expr_indices)))
        p_spear_BH_vec =c(p_spear_BH_vec,rep(p_spear_BH,length(expr_indices)))
        p_spear_Be_vec =c(p_spear_Be_vec,rep(p_spear_BH,length(expr_indices)))
        x = as.vector(as.numeric(edit_perctg))
        y = as.vector(as.numeric(expr_rpkm))
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = "Editing Percentage", ylab = paste0("log2(RPKM) ", sel_rnms[j]," Gene"),main =paste0(rownames(temp_drug_mat)[j], "\n", "Pearson R = ", cor_pear, " BH_Pearson = " ,p_pear_BH, "Bonf_Pearson = ",p_pear_Be,"\n","Spearman Rho = ", cor_spear, " BH_Spearman = " ,p_spear_BH, "Bonf_Spearman = ",p_spear_Be),cex.main = 0.65) 
        lines(pr.fit~x)
        next()
      }
      else if(is.na(expr_indices)){
        cor_pear_vec =c(cor_pear_vec,"gene name not found in expression array")
        cor_spear_vec =c(cor_spear_vec,"gene name not found in expression array")
        p_pear_vec =c(p_pear_vec,"gene name not found in expression array")
        p_spear_vec =c(p_spear_vec,"gene name not found in expression array")
        p_pear_BH_vec =c(p_pear_BH_vec,"gene name not found in expression array")
        p_pear_Be_vec =c(p_pear_Be_vec,"gene name not found in expression array")
        p_spear_BH_vec =c(p_spear_BH_vec,"gene name not found in expression array")
        p_spear_Be_vec =c(p_spear_Be_vec,"gene name not found in expression array")
        next()
        
      }
      else {
        expr_rpkm = RPKM[expr_indices,]
      }
      expr_rpkm = expr_rpkm[colnames(edit_perctg)]
      cor_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$estimate
      cor_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$estimate
      p_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$p.value
      p_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$p.value
      p_pear_BH = p.adjust(p = p_pear, method = "BH" )
      p_pear_Be = p.adjust(p=p_pear, method = "bonferroni")
      p_spear_BH = p.adjust(p = p_spear, method = "BH" )
      p_spear_Be = p.adjust(p=p_pear, method = "bonferroni")
      cor_pear_vec =c(cor_pear_vec,cor_pear)
      cor_spear_vec =c(cor_spear_vec,cor_spear)
      p_pear_vec =c(p_pear_vec,p_pear)
      p_spear_vec =c(p_spear_vec,p_spear)
      p_pear_BH_vec =c(p_pear_BH_vec,p_pear_BH)
      p_pear_Be_vec =c(p_pear_Be_vec,p_pear_Be)
      p_spear_BH_vec =c(p_spear_BH_vec,p_spear_BH )
      p_spear_Be_vec =c(p_spear_Be_vec,p_spear_Be)
      x = as.vector(as.numeric(edit_perctg))
      y = as.vector(as.numeric(expr_rpkm))
      fit = lm(y~x)
      pr.fit = predict(fit)
      plot(x,y, xlab = "Editing Percentage", ylab = paste0("log2(RPKM) ", sel_rnms[j]," Gene"),main =paste0(rownames(temp_drug_mat)[j], "\n", "Pearson R = ", cor_pear, " BH_Pearson = " ,p_pear_BH, "Bonf_Pearson = ",p_pear_Be,"\n","Spearman Rho = ", cor_spear, " BH_Spearman = " ,p_spear_BH, "Bonf_Spearman = ",p_spear_Be),cex.main = 0.65) 
      lines(pr.fit~x)
    }
    cor_pear_vec =c(cor_pear_vec,cor_pear_vec[!duplicated(cor_pear_vec)])
    cor_spear_vec =c(cor_spear_vec,cor_spear_vec[!duplicated(cor_spear_vec)])
    p_pear_vec =c(p_pear_vec,p_pear[!duplicated(p_pear_vec)])
    p_spear_vec =c(p_spear_vec,p_spear_vec[!duplicated(p_spear_vec)])
    p_pear_BH_vec =c(p_pear_BH_vec,p_pear_BH_vec[!duplicated(p_pear_BH_vec)])
    p_pear_Be_vec =c(p_pear_Be_vec,p_pear_Be_vec[!duplicated(p_pear_Be_vec)])
    p_spear_BH_vec =c(p_spear_BH_vec,p_spear_BH_vec[!duplicated(p_spear_BH_vec)])
    p_spear_Be_vec =c(p_spear_Be_vec,p_spear_Be_vec[!duplicated(p_spear_Be_vec)])
    rnms_output = rownames(output_frame) 
    tryCatch({
    output_frame = cbind(output_frame,cor_pear_vec ,
                        cor_spear_vec ,
                         p_pear_vec ,
                         p_spear_vec ,
                         p_pear_BH_vec ,
                         p_pear_Be_vec ,
                         p_spear_BH_vec ,
                         p_spear_Be_vec )
    fin_gnms = genenames[rownames(output_frame)]
    rownames(output_frame) =paste0(rownames(output_frame),"/",fin_gnms)
    write.csv(output_frame, paste0(drugname,"editing_vs_expression_log2RPKM.csv"))
    }, error = function(e){})
    dev.off()
  }
  c= c+1
}

dir.create("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_significant_editings_expression")
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_significant_editings_expression")
l = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlations_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlations_logplot/logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlations_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="logplotsVsCor.pdf"]
bool11 = grepl(x = nopath,pattern = "AG and TC Variants Only",fixed = T)
nopath = nopath[bool11]
l = l[bool11]
c= 1
RPKM = read.csv("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
rownames(RPKM) = RPKM$V1
genenames = RPKM$V2
genenames = genenames[c(-1,-2)]
my_cell_names = as.vector(as.character(RPKM[1,]))
my_cell_names = my_cell_names[c(-1,-2)]
ccle_cell_names = as.vector(as.character(RPKM[2,]))
ccle_cell_names = ccle_cell_names[c(-1,-2)]
RPKM = RPKM[c(-1,-2)]
RPKM = RPKM[c(-1,-2),]
names(genenames)  = rownames(RPKM)
for(i in 1:ncol(RPKM)){
  RPKM[,i] = as.numeric(RPKM[,i])
}
RPKM = log2(RPKM+1)
#c=16,
for (i in l){
  #i= l[1]
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_significant_editings_expression/",substr(justname, 1, nchar(justname)- 13), ".pdf")
  za =c("All Variants", "AG Variants", "AG and TC Variants Only")
  vec = c(grepl(x = nopath[c],pattern = "All Variants"), grepl(x = nopath[c],pattern ="AG Variants"),grepl(x = nopath[c],pattern = "AG and TC Variants Only"))
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = za[vec], fixed = T ))-2)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations/", drugname,"_correlation.csv")
  
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  drug_mat = drug_mat[,!is.na(drug_mat[2,])]
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
  drug_mat = drug_mat[-1,]
  rnms = rownames(drug_mat)
  drug_mat=as.data.frame(apply(drug_mat, 2, FUN =  as.numeric),row.names = rnms)
  cor_mat = as.data.frame(apply(cor_mat, 2, FUN =  as.numeric),row.names = rnms[-1])
  mat = read.csv(file=i, header = T, stringsAsFactors = F, row.names = 1, check.names = F)
  bb = is.na(mat$`Pearson R`)
  cc  =is.na(mat$`Spearman R`)
  mat = mat[!(bb|cc),]
  mat2 = mat[5:8]
  mat3 = colSums(mat2<0.05)
  
  if(nrow(mat2)==0){
    c= c+1
    next()
  }
  #dev.off()
  if(sum(mat3)==0){
    c = c+1
    next()
    
  }
  else{
    pdf(pdf_name,paper = "a4r")
    mat2  = mat2[rowSums(mat2<0.05)>0,]
    sel_rnms = rownames(mat2)
    
    temp_drug_mat = drug_mat[sel_rnms,]
    idx_sel_gnms = gregexpr(text = sel_rnms, pattern = "|" ,fixed = T )
    abs_idx_genes = list()
    c2 =1
    for(kj in 1:length(idx_sel_gnms)){
      idx_temp = idx_sel_gnms[[kj]]
      idx_temp = idx_temp[c(4,5)]
      idx_temp[1] = idx_temp[1]+1
      idx_temp[2] = idx_temp[2]-1
      sel_rnms[kj] = substr(sel_rnms[kj], idx_temp[1], idx_temp[2])
      at = which(toupper(genenames) == toupper(sel_rnms[kj]))
      ay = which(grepl(x = toupper(genenames), pattern = toupper(sel_rnms[kj]),fixed = T))
      gene_not_found = c()
      if((length(at)==0)& (length(ay)==0)){
        gene_not_found = c(gene_not_found, rownames(temp_drug_mat[kj]))
        abs_idx_genes[[c2]] = NA
      }
      else if (length(at)==0){
        abs_idx_genes[[c2]] = ay
      }
      else {
        abs_idx_genes[[c2]] =  at
        
      }
      c2 = c2+1
    }
    temp_rpkm = RPKM[unlist(abs_idx_genes)[!is.na(unlist(abs_idx_genes))],]
    cell_nms = intersect(colnames(temp_rpkm),colnames(drug_mat))
    temp_rpkm = temp_rpkm[cell_nms]
    if(length(cell_nms)< length(colnames(drug_mat))){
      notpresent = colnames(drug_mat)[!(cell_nms)]
      rownm = rownames(temp_rpkm)
      for(jk in 1:length(notpresent)){
        temp_rpkm = cbind(temp_rpkm,NA)
      }
      colnames(temp_rpkm)[(ncol(temp_rpkm)-length(notpresent)+1):(ncol(temp_rpkm))] = notpresent 
    }
    rownames(temp_rpkm)
    output_frame = rbind(temp_rpkm, temp_drug_mat)
    cor_pear_vec =c()
    cor_spear_vec =c()
    p_pear_vec =c()
    p_spear_vec =c()
    p_pear_BH_vec =c()
    p_pear_Be_vec =c()
    p_spear_BH_vec =c()
    p_spear_Be_vec =c()
    
    for(j in 1:nrow(temp_drug_mat)){
      edit_perctg = temp_drug_mat[j,]
      edit_perctg = edit_perctg[order(edit_perctg,decreasing = F)]
      expr_indices = abs_idx_genes[[j]]
      if(length(expr_indices)>1){
        expr_rpkm = RPKM[expr_indices,]
        expr_rpkm = colMeans(expr_rpkm)
        expr_rpkm = expr_rpkm[colnames(edit_perctg)]
        cor_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$estimate
        cor_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$estimate
        p_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$p.value
        p_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$p.value
        p_pear_BH = p.adjust(p = p_pear, method = "BH" )
        p_pear_Be = p.adjust(p=p_pear, method = "bonferroni")
        p_spear_BH = p.adjust(p = p_spear, method = "BH" )
        p_spear_Be = p.adjust(p=p_pear, method = "bonferroni")
        cor_pear_vec =c(cor_pear_vec,rep(cor_pear,length(expr_indices)))
        cor_spear_vec =c(cor_spear_vec,rep(cor_spear,length(expr_indices)))
        p_pear_vec =c(p_pear_vec,rep(p_pear,length(expr_indices)))
        p_spear_vec =c(p_spear_vec,rep(p_spear,length(expr_indices)))
        p_pear_BH_vec =c(p_pear_BH_vec,rep(p_pear_BH,length(expr_indices)))
        p_pear_Be_vec =c(p_pear_Be_vec,rep(p_pear_Be,length(expr_indices)))
        p_spear_BH_vec =c(p_spear_BH_vec,rep(p_spear_BH,length(expr_indices)))
        p_spear_Be_vec =c(p_spear_Be_vec,rep(p_spear_BH,length(expr_indices)))
        x = as.vector(as.numeric(edit_perctg))
        y = as.vector(as.numeric(expr_rpkm))
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = "Editing Percentage", ylab = paste0("log2(RPKM) ", sel_rnms[j]," Gene"),main =paste0(rownames(temp_drug_mat)[j], "\n", "Pearson R = ", cor_pear, " BH_Pearson = " ,p_pear_BH, "Bonf_Pearson = ",p_pear_Be,"\n","Spearman Rho = ", cor_spear, " BH_Spearman = " ,p_spear_BH, "Bonf_Spearman = ",p_spear_Be),cex.main = 0.65) 
        lines(pr.fit~x)
        next()
      }
      else if(is.na(expr_indices)){
        cor_pear_vec =c(cor_pear_vec,"gene name not found in expression array")
        cor_spear_vec =c(cor_spear_vec,"gene name not found in expression array")
        p_pear_vec =c(p_pear_vec,"gene name not found in expression array")
        p_spear_vec =c(p_spear_vec,"gene name not found in expression array")
        p_pear_BH_vec =c(p_pear_BH_vec,"gene name not found in expression array")
        p_pear_Be_vec =c(p_pear_Be_vec,"gene name not found in expression array")
        p_spear_BH_vec =c(p_spear_BH_vec,"gene name not found in expression array")
        p_spear_Be_vec =c(p_spear_Be_vec,"gene name not found in expression array")
        next()
        
      }
      else {
        expr_rpkm = RPKM[expr_indices,]
      }
      expr_rpkm = expr_rpkm[colnames(edit_perctg)]
      cor_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$estimate
      cor_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$estimate
      p_pear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "pearson",use = "na.or.complete")$p.value
      p_spear = cor.test(as.vector(as.numeric(edit_perctg)),as.vector(as.numeric(expr_rpkm)),method = "spearman",use = "na.or.complete")$p.value
      p_pear_BH = p.adjust(p = p_pear, method = "BH" )
      p_pear_Be = p.adjust(p=p_pear, method = "bonferroni")
      p_spear_BH = p.adjust(p = p_spear, method = "BH" )
      p_spear_Be = p.adjust(p=p_pear, method = "bonferroni")
      cor_pear_vec =c(cor_pear_vec,cor_pear)
      cor_spear_vec =c(cor_spear_vec,cor_spear)
      p_pear_vec =c(p_pear_vec,p_pear)
      p_spear_vec =c(p_spear_vec,p_spear)
      p_pear_BH_vec =c(p_pear_BH_vec,p_pear_BH)
      p_pear_Be_vec =c(p_pear_Be_vec,p_pear_Be)
      p_spear_BH_vec =c(p_spear_BH_vec,p_spear_BH )
      p_spear_Be_vec =c(p_spear_Be_vec,p_spear_Be)
      x = as.vector(as.numeric(edit_perctg))
      y = as.vector(as.numeric(expr_rpkm))
      fit = lm(y~x)
      pr.fit = predict(fit)
      plot(x,y, xlab = "Editing Percentage", ylab = paste0("log2(RPKM) ", sel_rnms[j]," Gene"),main =paste0(rownames(temp_drug_mat)[j], "\n", "Pearson R = ", cor_pear, " BH_Pearson = " ,p_pear_BH, "Bonf_Pearson = ",p_pear_Be,"\n","Spearman Rho = ", cor_spear, " BH_Spearman = " ,p_spear_BH, "Bonf_Spearman = ",p_spear_Be),cex.main = 0.65) 
      lines(pr.fit~x)
    }
    cor_pear_vec =c(cor_pear_vec,cor_pear_vec[!duplicated(cor_pear_vec)])
    cor_spear_vec =c(cor_spear_vec,cor_spear_vec[!duplicated(cor_spear_vec)])
    p_pear_vec =c(p_pear_vec,p_pear[!duplicated(p_pear_vec)])
    p_spear_vec =c(p_spear_vec,p_spear_vec[!duplicated(p_spear_vec)])
    p_pear_BH_vec =c(p_pear_BH_vec,p_pear_BH_vec[!duplicated(p_pear_BH_vec)])
    p_pear_Be_vec =c(p_pear_Be_vec,p_pear_Be_vec[!duplicated(p_pear_Be_vec)])
    p_spear_BH_vec =c(p_spear_BH_vec,p_spear_BH_vec[!duplicated(p_spear_BH_vec)])
    p_spear_Be_vec =c(p_spear_Be_vec,p_spear_Be_vec[!duplicated(p_spear_Be_vec)])
    rnms_output = rownames(output_frame) 
    tryCatch({
      output_frame = cbind(output_frame,cor_pear_vec ,
                           cor_spear_vec ,
                           p_pear_vec ,
                           p_spear_vec ,
                           p_pear_BH_vec ,
                           p_pear_Be_vec ,
                           p_spear_BH_vec ,
                           p_spear_Be_vec )
      fin_gnms = genenames[rownames(output_frame)]
      rownames(output_frame) =paste0(rownames(output_frame),"/",fin_gnms)
      write.csv(output_frame, paste0(drugname,"editing_vs_expression_log2RPKM.csv"))
    }, error = function(e){})
    dev.off()
  }
  c= c+1
}
