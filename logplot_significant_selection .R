setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot/logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="logplotsVsCor.pdf"]
c= 1
#c=16
for (i in l){
  #i= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot/Irinotecan AG and TC Variants Only _zeros_omitted_logplots.csv"
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlation_significants_selected/",substr(justname, 1, nchar(justname)- 13), "significant_selected.pdf")
  za =c("All Variants", "AG Variants", "AG and TC Variants Only")
  vec = c(grepl(x = nopath[c],pattern = "All Variants"), grepl(x = nopath[c],pattern ="AG Variants"),grepl(x = nopath[c],pattern = "AG and TC Variants Only"))
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = za[vec], fixed = T ))-2)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations/", drugname,"_correlation.csv")
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
  rnms = rownames(drug_mat)
  drug_mat=as.data.frame(apply(drug_mat, 2, FUN =  as.numeric),row.names = rnms)
  cor_mat = as.data.frame(apply(cor_mat, 2, FUN =  as.numeric),row.names = rnms[-1])
  mat = read.csv(file=i, header = T, stringsAsFactors = F, row.names = 1, check.names = F)
  mat2 = mat[5:8]
  mat3 = colSums(mat2<0.05)
  bb = is.na(mat2$`Pearson R`)
  cc  =is.na(mat2$`Spearman R`)
  mat2 = mat2[!(bb|cc),]
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
    mat2 = mat2[mat3>0]
    sig_names = list()
    c2 = 1
    for(j in colnames(mat2)){
      temp = mat2[c(j)]
      bool = temp<0.05
      sig = rownames(temp)[bool]
      temp = temp[temp<0.05]
      
      sig_names[[c2]] = sig
      names(sig_names)[c2] = j 
      c3 = 1
      for(k in sig){
        x = as.vector(as.numeric(drug_mat[c("Area Under Curve"), ]))
        y = as.vector(as.numeric(drug_mat[c(k),]))
        R = cor_mat[c(k),ifelse(test = grepl(x = j, pattern = "Pearson", fixed = T), yes = 1,no=3)]
        y=y[!is.na(x)]
        x =x[!is.na(x)]
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = paste0(drugname," Area Under Curve"), ylab = "RNA editing",main =paste0(j, "\n",k, " ", " ", data_type, " \nR = ",R, " corrected p = ", temp[c3]),cex.main = 0.65) 
        lines(pr.fit~x)
        c3 = c3+1
      } 
      c2= c2+1
    }
    if(length(sig_names)>1){
      names = intersect(sig_names[[1]],sig_names[[2]])
      if(length(names)==0){
        
      }
      else{
        plot.new()
        title(paste(paste(names,collapse = "\n"), "\nare significant for both Bonferroni and BH corrections"), cex.main = 0.75, outer=F)
      }
    }
    dev.off()
  }
  c= c+1
  
}

############# without zeros selection
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot/Without_zeros_logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="Without_zeros_logplotsVsCor.pdf"]
c= 1
#c=16
for (i in l){
  #i= l[1]
  i= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot/Irinotecan AG and TC Variants Only _zeros_omitted_logplots.csv"
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlation_significants_selected_without_zero/",substr(justname, 1, nchar(justname)- 13), "significant_selected.pdf")
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = " ", fixed = T ))-1)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted/", drugname,"_correlation.csv")
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
  rnms = rownames(drug_mat)
  drug_mat=as.data.frame(apply(drug_mat, 2, FUN =  as.numeric),row.names = rnms)
  cor_mat = as.data.frame(apply(cor_mat, 2, FUN =  as.numeric),row.names = rnms[-1])
  mat = read.csv(file=i, header = T, stringsAsFactors = F, row.names = 1, check.names = F)
  mat2 = mat[5:8]
  mat3 = colSums(mat2<0.05)
  #dev.off()(
  if(sum(mat3)==0){
    c = c+1
    next()
    
    
  }
  else{
    pdf(pdf_name,paper = "a4r")
    mat2 = mat2[mat3>0]
    sig_names = list()
    c2 = 1
    for(j in colnames(mat2)){
      temp = mat2[c(j)]
      bool = temp<0.05
      sig = rownames(temp)[bool]
      temp = temp[temp<0.05]
      sig_names[[c2]] = sig
      names(sig_names)[c2] = j 
      c3 =1
      for(k in sig){
        print(temp[c3])
        print(k)
        print(j)
        x = as.vector(as.numeric(drug_mat[c("Area Under Curve"), ]))
        y = as.vector(as.numeric(drug_mat[c(k),]))
        xy=rbind(x,y)
        booool = colSums(is.na(xy))
        booool = booool==0
        xy = xy[,booool]
        x= as.vector(as.numeric(xy[1,]))
        y= as.vector(as.numeric(xy[2,]))
        R = cor_mat[c(k),ifelse(test = grepl(x = j, pattern = "Pearson", fixed = T), yes = 1,no=3)]
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = paste0(drugname," Area Under Curve"), ylab = "RNA editing",main =paste0(j, "\n",k, " ", " ", data_type, " \nR = ",R, " corrected p = ", temp[c3]),cex.main = 0.65) 
        lines(pr.fit~x)
        c3 = c3+1
      } 
      c2= c2+1
    }
    if(length(sig_names)>1){
      names = intersect(sig_names[[1]],sig_names[[2]])
      if(length(names)==0){
        
      }
      else{
        plot.new()
        title(paste(paste(names,collapse = "\n"), "\nare significant for both Bonferroni and BH corrections"), cex.main = 0.75, outer=F)
      }
    }
    dev.off()
  }
  c= c+1
  
}
dev.off(
)


######## 3'UTR 
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot/logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="logplotsVsCor.pdf"]
c= 1
#c=16
for (i in l){
  #i= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot/Irinotecan AG and TC Variants Only _zeros_omitted_logplots.csv"
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlation_significants_selected/",substr(justname, 1, nchar(justname)- 13), "significant_selected.pdf")
  za =c("All Variants", "AG Variants", "AG and TC Variants Only")
  vec = c(grepl(x = nopath[c],pattern = "All Variants"), grepl(x = nopath[c],pattern ="AG Variants"),grepl(x = nopath[c],pattern = "AG and TC Variants Only"))
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = za[vec], fixed = T ))-2)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations/", drugname,"_correlation.csv")
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
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
    mat2 = mat2[mat3>0]
    sig_names = list()
    c2 = 1
    for(j in colnames(mat2)){
      temp = mat2[c(j)]
      bool = temp<0.05
      sig = rownames(temp)[bool]
      temp = temp[temp<0.05]
      
      sig_names[[c2]] = sig
      names(sig_names)[c2] = j 
      c3 = 1
      for(k in sig){
        x = as.vector(as.numeric(drug_mat[c("Area Under Curve"), ]))
        y = as.vector(as.numeric(drug_mat[c(k),]))
        R = cor_mat[c(k),ifelse(test = grepl(x = j, pattern = "Pearson", fixed = T), yes = 1,no=3)]
        y=y[!is.na(x)]
        x =x[!is.na(x)]
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = paste0(drugname," Area Under Curve"), ylab = "RNA editing",main =paste0(j, "\n",k, " ", " ", data_type, " \nR = ",R, " corrected p = ", temp[c3]),cex.main = 0.65) 
        lines(pr.fit~x)
        c3 = c3+1
      } 
      c2= c2+1
    }
    if(length(sig_names)>1){
      names = intersect(sig_names[[1]],sig_names[[2]])
      if(length(names)==0){
        
      }
      else{
        plot.new()
        title(paste(paste(names,collapse = "\n"), "\nare significant for both Bonferroni and BH corrections"), cex.main = 0.75, outer=F)
      }
    }
    dev.off()
  }
  c= c+1
  
}

dev.off()

############## Introns
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlations_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlations_logplot/logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlations_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="logplotsVsCor.pdf"]
c= 1
#c=16
for (i in l){
  #i= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot/Irinotecan AG and TC Variants Only _zeros_omitted_logplots.csv"
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlation_significants_selected/",substr(justname, 1, nchar(justname)- 13), "significant_selected.pdf")
  za =c("All Variants", "AG Variants", "AG and TC Variants Only")
  vec = c(grepl(x = nopath[c],pattern = "All Variants"), grepl(x = nopath[c],pattern ="AG Variants"),grepl(x = nopath[c],pattern = "AG and TC Variants Only"))
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = za[vec], fixed = T ))-2)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations/", drugname,"_correlation.csv")
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
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
    mat2 = mat2[mat3>0]
    sig_names = list()
    c2 = 1
    for(j in colnames(mat2)){
      temp = mat2[c(j)]
      bool = temp<0.05
      sig = rownames(temp)[bool]
      temp = temp[temp<0.05]
      
      sig_names[[c2]] = sig
      names(sig_names)[c2] = j 
      c3 = 1
      for(k in sig){
        x = as.vector(as.numeric(drug_mat[c("Area Under Curve"), ]))
        y = as.vector(as.numeric(drug_mat[c(k),]))
        R = cor_mat[c(k),ifelse(test = grepl(x = j, pattern = "Pearson", fixed = T), yes = 1,no=3)]
        y=y[!is.na(x)]
        x =x[!is.na(x)]
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = paste0(drugname," Area Under Curve"), ylab = "RNA editing",main =paste0(j, "\n",k, " ", " ", data_type, " \nR = ",R, " corrected p = ", temp[c3]),cex.main = 0.65) 
        lines(pr.fit~x)
        c3 = c3+1
      } 
      c2= c2+1
    }
    if(length(sig_names)>1){
      names = intersect(sig_names[[1]],sig_names[[2]])
      if(length(names)==0){
        
      }
      else{
        plot.new()
        title(paste(paste(names,collapse = "\n"), "\nare significant for both Bonferroni and BH corrections"), cex.main = 0.75, outer=F)
      }
    }
    dev.off()
  }
  c= c+1
  
}

dev.off()

######### coding region
dir.create("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlation_significants_selected")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlations_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlations_logplot/logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlations_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="logplotsVsCor.pdf"]
c= 1
#c=16
for (i in l){
  #i= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot/Irinotecan AG and TC Variants Only _zeros_omitted_logplots.csv"
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlation_significants_selected/",substr(justname, 1, nchar(justname)- 13), "significant_selected.pdf")
  za =c("All Variants", "AG Variants", "AG and TC Variants Only")
  vec = c(grepl(x = nopath[c],pattern = "All Variants"), grepl(x = nopath[c],pattern ="AG Variants"),grepl(x = nopath[c],pattern = "AG and TC Variants Only"))
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = za[vec], fixed = T ))-2)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations/", drugname,"_correlation.csv")
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
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
    mat2 = mat2[mat3>0]
    sig_names = list()
    c2 = 1
    for(j in colnames(mat2)){
      temp = mat2[c(j)]
      bool = temp<0.05
      sig = rownames(temp)[bool]
      temp = temp[temp<0.05]
      
      sig_names[[c2]] = sig
      names(sig_names)[c2] = j 
      c3 = 1
      for(k in sig){
        x = as.vector(as.numeric(drug_mat[c("Area Under Curve"), ]))
        y = as.vector(as.numeric(drug_mat[c(k),]))
        R = cor_mat[c(k),ifelse(test = grepl(x = j, pattern = "Pearson", fixed = T), yes = 1,no=3)]
        y=y[!is.na(x)]
        x =x[!is.na(x)]
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = paste0(drugname," Area Under Curve"), ylab = "RNA editing",main =paste0(j, "\n",k, " ", " ", data_type, " \nR = ",R, " corrected p = ", temp[c3]),cex.main = 0.65) 
        lines(pr.fit~x)
        c3 = c3+1
      } 
      c2= c2+1
    }
    if(length(sig_names)>1){
      names = intersect(sig_names[[1]],sig_names[[2]])
      if(length(names)==0){
        
      }
      else{
        plot.new()
        title(paste(paste(names,collapse = "\n"), "\nare significant for both Bonferroni and BH corrections"), cex.main = 0.75, outer=F)
      }
    }
    dev.off()
  }
  c= c+1
  
}

dev.off()
#5UTR
dir.create("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlation_significants_selected")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlations_logplot",all.files = F, full.names = T)
l = l[l!= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlations_logplot/logplotsVsCor.pdf"]
nopath = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlations_logplot",all.files = F, full.names = F)
nopath = nopath[nopath!="logplotsVsCor.pdf"]
c= 1
#c=16
for (i in l){
  #i= "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot/Irinotecan AG and TC Variants Only _zeros_omitted_logplots.csv"
  justname = nopath[c]
  pdf_name = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlation_significants_selected/",substr(justname, 1, nchar(justname)- 13), "significant_selected.pdf")
  za =c("All Variants", "AG Variants", "AG and TC Variants Only")
  vec = c(grepl(x = nopath[c],pattern = "All Variants"), grepl(x = nopath[c],pattern ="AG Variants"),grepl(x = nopath[c],pattern = "AG and TC Variants Only"))
  drugname = substr(nopath[c], 1, (regexpr(text = nopath[c],pattern = za[vec], fixed = T ))-2)
  spaces = gregexpr(text = justname,pattern =  " ", fixed =T)[[1]]
  last_space = spaces[length(spaces)]
  data_type = substr(justname, nchar(drugname)+2, last_space-1)
  corr_path = paste0("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations/", drugname,"_correlation.csv")
  drug_mat  = read.csv(corr_path, row.names = 1, header = T, check.names = F, stringsAsFactors = F)
  cor_mat = drug_mat[3:nrow(drug_mat),(ncol(drug_mat)-3):ncol(drug_mat)]
  drug_mat = drug_mat[2:nrow(drug_mat),1:(ncol(drug_mat)-4)]
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
    mat2 = mat2[mat3>0]
    sig_names = list()
    c2 = 1
    for(j in colnames(mat2)){
      temp = mat2[c(j)]
      bool = temp<0.05
      sig = rownames(temp)[bool]
      temp = temp[temp<0.05]
      
      sig_names[[c2]] = sig
      names(sig_names)[c2] = j 
      c3 = 1
      for(k in sig){
        x = as.vector(as.numeric(drug_mat[c("Area Under Curve"), ]))
        y = as.vector(as.numeric(drug_mat[c(k),]))
        R = cor_mat[c(k),ifelse(test = grepl(x = j, pattern = "Pearson", fixed = T), yes = 1,no=3)]
        y=y[!is.na(x)]
        x =x[!is.na(x)]
        fit = lm(y~x)
        pr.fit = predict(fit)
        plot(x,y, xlab = paste0(drugname," Area Under Curve"), ylab = "RNA editing",main =paste0(j, "\n",k, " ", " ", data_type, " \nR = ",R, " corrected p = ", temp[c3]),cex.main = 0.65) 
        lines(pr.fit~x)
        c3 = c3+1
      } 
      c2= c2+1
    }
    if(length(sig_names)>1){
      names = intersect(sig_names[[1]],sig_names[[2]])
      if(length(names)==0){
        
      }
      else{
        plot.new()
        title(paste(paste(names,collapse = "\n"), "\nare significant for both Bonferroni and BH corrections"), cex.main = 0.75, outer=F)
      }
    }
    dev.off()
  }
  c= c+1
  
}

dev.off()
