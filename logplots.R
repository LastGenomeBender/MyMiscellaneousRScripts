#dir.create("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compound_correlations_logplot")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = T)
ll = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = F)
c=1
pdf(file = "logplotsVsCor.pdf")
for(i in l ){
  #i= l[1]
  drug_name = substr(ll[c], 1, (nchar(ll[c])-nchar("_correlation.csv")))
  data = read.csv(file = i, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  data2 = data[c((ncol(data)-3):ncol(data))]
  data2 = data2[c(-1,-2),]
  data2 = apply(data2, 2, as.numeric)
  data2 = as.data.frame(data2,row.names = rownames(data)[c(-1,-2)])
  data_AG = data2[grepl("\\|AG\\|", rownames(data2)),]
  data_TC_AG = data2[((grepl("\\|AG\\|", rownames(data2))) | (grepl("\\|TC\\|", rownames(data2))) ),]
  datas = list()
  datas[[1]] = data2
  datas[[2]] = data_AG
  datas[[3]] = data_TC_AG
  names_dataset = c("All Variants", "AG Variants", "AG and TC Variants Only")
  e=1
  for(j in 1:length(datas)){
    nm_dtst = names_dataset[e]
    #j = 1
    mat = datas[[j]]
    bb = is.na(mat$`Pearson R`)
    cc  =is.na(mat$`Spearman R`)
    mat = mat[!(bb|cc),]
    if(nrow(mat)==0){
      print(e)
      e = e+1
      next()
    }
    mat$BH_Pearson = p.adjust(mat$`Pearson Pval`, method = "BH")
    mat$Bonferroni_Pearson = p.adjust(mat$`Pearson Pval`, method = "bonferroni")
    mat$BH_Spearman = p.adjust(mat$`Spearman Pval`, method = "BH")
    mat$Bonferroni_Spearman = p.adjust(mat$`Spearman Pval`, method = "bonferroni")
    minus_log10=cbind(-log10(mat[,c("Pearson Pval","Spearman Pval")]), -log10(mat[c((ncol(mat)-3):ncol(mat))]))
    colnames(minus_log10) = paste("-log10(",colnames(minus_log10),")")
    mat = cbind(mat, minus_log10)
    fname = paste(drug_name, nm_dtst ,"logplots.csv")
    write.csv(mat, fname)
    for(k in 1:2){
      #k=1
      if ( k ==1 ){
        x_axis = mat$`Pearson R`
        y_axis_coors = c(9,11,12)
      } else if(k==2){
        x_axis = mat$`Spearman R`
        y_axis_coors = c(10,13,14)
      }
      for(o in y_axis_coors){
        #o = 10
        
        plot(as.numeric(x_axis), as.numeric(mat[,o]), xlab= ifelse(test = k==1, yes = "Pearson R", no = "Spearman Rho"), ylab = colnames(mat)[o])
        abline(h = -log10(0.05))
        title(paste0(drug_name, " ", nm_dtst))
      }
    }
    
    e = e+1
  }
  c = c+1
}
dev.off()

############ without zeros
dir.create("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted",all.files = F, full.names = T)
ll = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted",all.files = F, full.names = F)
c=1
pdf(file = "Without_zeros_logplotsVsCor.pdf")
for(i in l ){
  #i= l[1]
  drug_name = substr(ll[c], 1, (nchar(ll[c])-nchar("_correlation.csv")))
  data = read.csv(file = i, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  data2 = data[c((ncol(data)-3):ncol(data))]
  data2 = data2[c(-1,-2),]
  data2 = apply(data2, 2, as.numeric)
  data2 = as.data.frame(data2,row.names = rownames(data)[c(-1,-2)])
  data_AG = data2[grepl("\\|AG\\|", rownames(data2)),]
  data_TC_AG = data2[((grepl("\\|AG\\|", rownames(data2))) | (grepl("\\|TC\\|", rownames(data2))) ),]
  datas = list()
  datas[[1]] = data2
  datas[[2]] = data_AG
  datas[[3]] = data_TC_AG
  names_dataset = c("All Variants", "AG Variants", "AG and TC Variants Only")
  e=1
  for(j in 1:length(datas)){
    nm_dtst = names_dataset[e]
    #j = 1
    mat = datas[[j]]
    bb = is.na(mat$`Pearson R`)
    cc  =is.na(mat$`Spearman R`)
    dd = mat[!(bb|cc),]
    if(nrow(mat)==0){
      print(e)
      e = e+1
      next()
    }
    mat$BH_Pearson = p.adjust(mat$`Pearson Pval`, method = "BH")
    mat$Bonferroni_Pearson = p.adjust(mat$`Pearson Pval`, method = "bonferroni")
    mat$BH_Spearman = p.adjust(mat$`Spearman Pval`, method = "BH")
    mat$Bonferroni_Spearman = p.adjust(mat$`Spearman Pval`, method = "bonferroni")
    minus_log10=cbind(-log10(mat[,c("Pearson Pval","Spearman Pval")]), -log10(mat[c((ncol(mat)-3):ncol(mat))]))
    colnames(minus_log10) = paste("-log10(",colnames(minus_log10),")")
    mat = cbind(mat, minus_log10)
    fname = paste(drug_name, nm_dtst,"_zeros_omitted_logplots.csv")
    write.csv(mat, fname)
    for(k in 1:2){
      #k=1
      if ( k ==1 ){
        x_axis = mat$`Pearson R`
        y_axis_coors = c(9,11,12)
      } else if(k==2){
        x_axis = mat$`Spearman R`
        y_axis_coors = c(10,13,14)
      }
      for(o in y_axis_coors){
        #o = 10
        
        plot(as.numeric(x_axis), as.numeric(mat[,o]), xlab= ifelse(test = k==1, yes = "Pearson R", no = "Spearman Rho"), ylab = colnames(mat)[o])
        abline(h = -log10(0.05))
        title(paste0(drug_name, " ", nm_dtst))
      }
    }
    print(e)
    e = e+1
  }
  c = c+1
}
dev.off()
#############3UTR only

dir.create("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = T)
ll = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = F)
c=1
pdf(file = "logplotsVsCor.pdf")
for(i in l ){
  #i= l[1]
  drug_name = substr(ll[c], 1, (nchar(ll[c])-nchar("_correlation.csv")))
  data = read.csv(file = i, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  data2 = data[c((ncol(data)-3):ncol(data))]
  data2 = data2[c(-1,-2),]
  data2 = apply(data2, 2, as.numeric)
  data2 = as.data.frame(data2,row.names = rownames(data)[c(-1,-2)])
  data2 = data2[grepl(x = rownames(data2),pattern = "3_prime_UTR_variant",fixed = T),]
  data_AG = data2[grepl("\\|AG\\|", rownames(data2)),]
  data_TC_AG = data2[((grepl("\\|AG\\|", rownames(data2))) | (grepl("\\|TC\\|", rownames(data2))) ),]
  datas = list()
  datas[[1]] = data2
  datas[[2]] = data_AG
  datas[[3]] = data_TC_AG
  names_dataset = c("All Variants", "AG Variants", "AG and TC Variants Only")
  e=1
  for(j in 1:length(datas)){
    nm_dtst = names_dataset[e]
    #j = 1
    mat = datas[[j]]
    bb = is.na(mat$`Pearson R`)
    cc  =is.na(mat$`Spearman R`)
    mat = mat[!(bb|cc),]
    if(nrow(mat)==0){
      print(e)
      e = e+1
      next()
    }
    mat$BH_Pearson = p.adjust(mat$`Pearson Pval`, method = "BH")
    mat$Bonferroni_Pearson = p.adjust(mat$`Pearson Pval`, method = "bonferroni")
    mat$BH_Spearman = p.adjust(mat$`Spearman Pval`, method = "BH")
    mat$Bonferroni_Spearman = p.adjust(mat$`Spearman Pval`, method = "bonferroni")
    minus_log10=cbind(-log10(mat[,c("Pearson Pval","Spearman Pval")]), -log10(mat[c((ncol(mat)-3):ncol(mat))]))
    colnames(minus_log10) = paste("-log10(",colnames(minus_log10),")")
    mat = cbind(mat, minus_log10)
    fname = paste(drug_name, nm_dtst ,"logplots.csv")
    write.csv(mat, fname)
    for(k in 1:2){
      #k=1
      if ( k ==1 ){
        x_axis = mat$`Pearson R`
        y_axis_coors = c(9,11,12)
      } else if(k==2){
        x_axis = mat$`Spearman R`
        y_axis_coors = c(10,13,14)
      }
      for(o in y_axis_coors){
        #o = 10
        
        plot(as.numeric(x_axis), as.numeric(mat[,o]), xlab= ifelse(test = k==1, yes = "Pearson R", no = "Spearman Rho"), ylab = colnames(mat)[o])
        abline(h = -log10(0.05))
        title(paste0(drug_name, " ", nm_dtst))
      }
    }
    
    e = e+1
  }
  c = c+1
}
dev.off()

############ 3utr wo zeros

dir.create("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted",all.files = F, full.names = T)
ll = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations_Zeros_Omitted",all.files = F, full.names = F)
c=1
pdf(file = "Without_zeros_logplotsVsCor.pdf")
for(i in l ){
  #i= l[1]
  drug_name = substr(ll[c], 1, (nchar(ll[c])-nchar("_correlation.csv")))
  data = read.csv(file = i, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  data2 = data[c((ncol(data)-3):ncol(data))]
  data2 = data2[c(-1,-2),]
  data2 = apply(data2, 2, as.numeric)
  data2 = as.data.frame(data2,row.names = rownames(data)[c(-1,-2)])
  data_AG = data2[grepl("\\|AG\\|", rownames(data2)),]
  data_TC_AG = data2[((grepl("\\|AG\\|", rownames(data2))) | (grepl("\\|TC\\|", rownames(data2))) ),]
  datas = list()
  datas[[1]] = data2
  datas[[2]] = data_AG
  datas[[3]] = data_TC_AG
  names_dataset = c("All Variants", "AG Variants", "AG and TC Variants Only")
  e=1
  for(j in 1:length(datas)){
    nm_dtst = names_dataset[e]
    #j = 1
    mat = datas[[j]]
    bb = is.na(mat$`Pearson R`)
    cc  =is.na(mat$`Spearman R`)
    dd = mat[!(bb|cc),]
    if(nrow(mat)==0){
      print(e)
      e = e+1
      next()
    }
    mat$BH_Pearson = p.adjust(mat$`Pearson Pval`, method = "BH")
    mat$Bonferroni_Pearson = p.adjust(mat$`Pearson Pval`, method = "bonferroni")
    mat$BH_Spearman = p.adjust(mat$`Spearman Pval`, method = "BH")
    mat$Bonferroni_Spearman = p.adjust(mat$`Spearman Pval`, method = "bonferroni")
    minus_log10=cbind(-log10(mat[,c("Pearson Pval","Spearman Pval")]), -log10(mat[c((ncol(mat)-3):ncol(mat))]))
    colnames(minus_log10) = paste("-log10(",colnames(minus_log10),")")
    mat = cbind(mat, minus_log10)
    fname = paste(drug_name, nm_dtst,"_zeros_omitted_logplots.csv")
    write.csv(mat, fname)
    for(k in 1:2){
      #k=1
      if ( k ==1 ){
        x_axis = mat$`Pearson R`
        y_axis_coors = c(9,11,12)
      } else if(k==2){
        x_axis = mat$`Spearman R`
        y_axis_coors = c(10,13,14)
      }
      for(o in y_axis_coors){
        #o = 10
        
        plot(as.numeric(x_axis), as.numeric(mat[,o]), xlab= ifelse(test = k==1, yes = "Pearson R", no = "Spearman Rho"), ylab = colnames(mat)[o])
        abline(h = -log10(0.05))
        title(paste0(drug_name, " ", nm_dtst))
      }
    }
    print(e)
    e = e+1
  }
  c = c+1
}
dev.off()



dir.create("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlations_logplot")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = T)
ll = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = F)
c=1
pdf(file = "logplotsVsCor.pdf")
for(i in l ){
  #i= l[1]
  drug_name = substr(ll[c], 1, (nchar(ll[c])-nchar("_correlation.csv")))
  data = read.csv(file = i, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  data2 = data[c((ncol(data)-3):ncol(data))]
  data2 = data2[c(-1,-2),]
  data2 = apply(data2, 2, as.numeric)
  data2 = as.data.frame(data2,row.names = rownames(data)[c(-1,-2)])
  data2 = data2[grepl(x = rownames(data2),pattern = "intron_variant",fixed = T),]
  data_AG = data2[grepl("\\|AG\\|", rownames(data2)),]
  data_TC_AG = data2[((grepl("\\|AG\\|", rownames(data2))) | (grepl("\\|TC\\|", rownames(data2))) ),]
  datas = list()
  datas[[1]] = data2
  datas[[2]] = data_AG
  datas[[3]] = data_TC_AG
  names_dataset = c("All Variants", "AG Variants", "AG and TC Variants Only")
  e=1
  for(j in 1:length(datas)){
    nm_dtst = names_dataset[e]
    #j = 1
    mat = datas[[j]]
    bb = is.na(mat$`Pearson R`)
    cc  =is.na(mat$`Spearman R`)
    mat = mat[!(bb|cc),]
    if(nrow(mat)==0){
      print(e)
      e = e+1
      next()
    }
    mat$BH_Pearson = p.adjust(mat$`Pearson Pval`, method = "BH")
    mat$Bonferroni_Pearson = p.adjust(mat$`Pearson Pval`, method = "bonferroni")
    mat$BH_Spearman = p.adjust(mat$`Spearman Pval`, method = "BH")
    mat$Bonferroni_Spearman = p.adjust(mat$`Spearman Pval`, method = "bonferroni")
    minus_log10=cbind(-log10(mat[,c("Pearson Pval","Spearman Pval")]), -log10(mat[c((ncol(mat)-3):ncol(mat))]))
    colnames(minus_log10) = paste("-log10(",colnames(minus_log10),")")
    mat = cbind(mat, minus_log10)
    fname = paste(drug_name, nm_dtst ,"logplots.csv")
    write.csv(mat, fname)
    for(k in 1:2){
      #k=1
      if ( k ==1 ){
        x_axis = mat$`Pearson R`
        y_axis_coors = c(9,11,12)
      } else if(k==2){
        x_axis = mat$`Spearman R`
        y_axis_coors = c(10,13,14)
      }
      for(o in y_axis_coors){
        #o = 10
        
        plot(as.numeric(x_axis), as.numeric(mat[,o]), xlab= ifelse(test = k==1, yes = "Pearson R", no = "Spearman Rho"), ylab = colnames(mat)[o])
        abline(h = -log10(0.05))
        title(paste0(drug_name, " ", nm_dtst))
      }
    }
    
    e = e+1
  }
  c = c+1
}
dev.off()

#### coding
dir.create("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlations_logplot")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = T)
ll = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = F)
c=1
pdf(file = "logplotsVsCor.pdf")
for(i in l ){
  #i= l[1]
  drug_name = substr(ll[c], 1, (nchar(ll[c])-nchar("_correlation.csv")))
  data = read.csv(file = i, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  data2 = data[c((ncol(data)-3):ncol(data))]
  data2 = data2[c(-1,-2),]
  data2 = apply(data2, 2, as.numeric)
  data2 = as.data.frame(data2,row.names = rownames(data)[c(-1,-2)])
  data2 = data2[grepl(x = rownames(data2),pattern = "missense_variant",fixed = T)|grepl(x = rownames(data2),pattern = "gain",fixed = T),]
  data_AG = data2[grepl("\\|AG\\|", rownames(data2)),]
  data_TC_AG = data2[((grepl("\\|AG\\|", rownames(data2))) | (grepl("\\|TC\\|", rownames(data2))) ),]
  datas = list()
  datas[[1]] = data2
  datas[[2]] = data_AG
  datas[[3]] = data_TC_AG
  names_dataset = c("All Variants", "AG Variants", "AG and TC Variants Only")
  e=1
  for(j in 1:length(datas)){
    nm_dtst = names_dataset[e]
    #j = 1
    mat = datas[[j]]
    bb = is.na(mat$`Pearson R`)
    cc  =is.na(mat$`Spearman R`)
    mat = mat[!(bb|cc),]
    if(nrow(mat)==0){
      print(e)
      e = e+1
      next()
    }
    mat$BH_Pearson = p.adjust(mat$`Pearson Pval`, method = "BH")
    mat$Bonferroni_Pearson = p.adjust(mat$`Pearson Pval`, method = "bonferroni")
    mat$BH_Spearman = p.adjust(mat$`Spearman Pval`, method = "BH")
    mat$Bonferroni_Spearman = p.adjust(mat$`Spearman Pval`, method = "bonferroni")
    minus_log10=cbind(-log10(mat[,c("Pearson Pval","Spearman Pval")]), -log10(mat[c((ncol(mat)-3):ncol(mat))]))
    colnames(minus_log10) = paste("-log10(",colnames(minus_log10),")")
    mat = cbind(mat, minus_log10)
    fname = paste(drug_name, nm_dtst ,"logplots.csv")
    write.csv(mat, fname)
    for(k in 1:2){
      #k=1
      if ( k ==1 ){
        x_axis = mat$`Pearson R`
        y_axis_coors = c(9,11,12)
      } else if(k==2){
        x_axis = mat$`Spearman R`
        y_axis_coors = c(10,13,14)
      }
      for(o in y_axis_coors){
        #o = 10
        
        plot(as.numeric(x_axis), as.numeric(mat[,o]), xlab= ifelse(test = k==1, yes = "Pearson R", no = "Spearman Rho"), ylab = colnames(mat)[o])
        abline(h = -log10(0.05))
        title(paste0(drug_name, " ", nm_dtst))
      }
    }
    
    e = e+1
  }
  c = c+1
}
dev.off()

#5UTR
dir.create("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlations_logplot")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_Coompound_correlations_logplot")
l = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = T)
ll = list.files(path = "/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Compouds_correlations",all.files = F, full.names = F)
c=1
pdf(file = "logplotsVsCor.pdf")
for(i in l ){
  #i= l[1]
  drug_name = substr(ll[c], 1, (nchar(ll[c])-nchar("_correlation.csv")))
  data = read.csv(file = i, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  data2 = data[c((ncol(data)-3):ncol(data))]
  data2 = data2[c(-1,-2),]
  data2 = apply(data2, 2, as.numeric)
  data2 = as.data.frame(data2,row.names = rownames(data)[c(-1,-2)])
  data2 = data2[grepl(x = rownames(data2),pattern = "5_prime_UTR_variant",fixed = T),]
  data_AG = data2[grepl("\\|AG\\|", rownames(data2)),]
  data_TC_AG = data2[((grepl("\\|AG\\|", rownames(data2))) | (grepl("\\|TC\\|", rownames(data2))) ),]
  datas = list()
  datas[[1]] = data2
  datas[[2]] = data_AG
  datas[[3]] = data_TC_AG
  names_dataset = c("All Variants", "AG Variants", "AG and TC Variants Only")
  e=1
  for(j in 1:length(datas)){
    nm_dtst = names_dataset[e]
    #j = 1
    mat = datas[[j]]
    bb = is.na(mat$`Pearson R`)
    cc  =is.na(mat$`Spearman R`)
    mat = mat[!(bb|cc),]
    if(nrow(mat)==0){
      print(e)
      e = e+1
      next()
    }
    mat$BH_Pearson = p.adjust(mat$`Pearson Pval`, method = "BH")
    mat$Bonferroni_Pearson = p.adjust(mat$`Pearson Pval`, method = "bonferroni")
    mat$BH_Spearman = p.adjust(mat$`Spearman Pval`, method = "BH")
    mat$Bonferroni_Spearman = p.adjust(mat$`Spearman Pval`, method = "bonferroni")
    minus_log10=cbind(-log10(mat[,c("Pearson Pval","Spearman Pval")]), -log10(mat[c((ncol(mat)-3):ncol(mat))]))
    colnames(minus_log10) = paste("-log10(",colnames(minus_log10),")")
    mat = cbind(mat, minus_log10)
    fname = paste(drug_name, nm_dtst ,"logplots.csv")
    write.csv(mat, fname)
    for(k in 1:2){
      #k=1
      if ( k ==1 ){
        x_axis = mat$`Pearson R`
        y_axis_coors = c(9,11,12)
      } else if(k==2){
        x_axis = mat$`Spearman R`
        y_axis_coors = c(10,13,14)
      }
      for(o in y_axis_coors){
        #o = 10
        
        plot(as.numeric(x_axis), as.numeric(mat[,o]), xlab= ifelse(test = k==1, yes = "Pearson R", no = "Spearman Rho"), ylab = colnames(mat)[o])
        abline(h = -log10(0.05))
        title(paste0(drug_name, " ", nm_dtst))
      }
    }
    
    e = e+1
  }
  c = c+1
}
dev.off()
