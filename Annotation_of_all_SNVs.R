###### functional annotation of GIREMI output
for(j in list.files("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/A-to-I GIREMI Colon CCLE", full.names = F)){
  setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/A-to-I GIREMI Colon CCLE")
  library(tibble)
  library(vcfR)
  name = substr(j,1,(nchar(j)-27))
  gir=read.csv(j,row.names = 1, header = T, stringsAsFactors = F, check.names = F)
  vcfflnm = paste0("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/dbsn_and_funcannot_vcf/",name,"_dbsnp_n_fun_annotated.vcf")
  vcf = read.vcfR(vcfflnm)
  a = getINFO(vcf)
  b = getPOS(vcf)
  c = getCHROM(vcf)
  ala = cbind(a,b,c)
  i1 = seq(1,nrow(gir),1)
  for (i in i1){
    chr = gir$chr[i]
    pos = gir$coordinate[i]
    bool = (c==chr&b==pos)
    c_hat = ala[bool,]
    info = c_hat[1]
    anno=unlist(strsplit(strsplit(info, split = "ANN")[[1]][2], split=","))
    anno2 = strsplit(anno, split="\\|")
    anno4 = sapply(X = anno2, function(x) paste(x[2],x[4]))
    anno4 = anno4[!duplicated(anno4)]
    anno5 = paste(anno4, collapse = ",")
    gir$func_anno[i] = anno5
  }
  fffname = paste0(name,"_final_funcanno_usethis.csv" )
  setwd("/media/ahadli/Data/Farid/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/a-to-i_all_func_annotated")
  write.csv(gir, file = fffname)
}
