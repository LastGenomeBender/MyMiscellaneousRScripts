setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/Compound_correlations_logplot")
l = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/Compouds_correlations",all.files = F, full.names = T)
ll = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/Compouds_correlations",all.files = F, full.names = F)
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
setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/Compouds_correlations_Zeros_Omitted_logplot")
l = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/Compouds_correlations_Zeros_Omitted",all.files = F, full.names = T)
ll = list.files(path = "/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/Compouds_correlations_Zeros_Omitted",all.files = F, full.names = F)
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
