##### mutational profile ofthe crc cell lines of ccle

# tp53, APC, BRAF, KRAS, AXIN1 ,-2, beta catenin, smad4

## read the point mutation data,
point_mut = read.table("/home/ahadli/Desktop/CCLE_DepMap_18q3_maf_20180718.txt",header = T,stringsAsFactors = F, check.names = F,sep = "\t")
edit = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv", row.names = 1, header = T, stringsAsFactors = F, check.names = F)
editnms = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/cell_line_names.csv", header = T, stringsAsFactors = F, row.names = 1, check.names = F)
editnms = editnms[colnames(edit),]
class(gregexpr(text="aaadaad", "d", fixed = T)[[1]])
vec_expnms = point_mut$Tumor_Sample_Barcode[!duplicated(point_mut$Tumor_Sample_Barcode)]
underline_indices = regexpr(text = vec_expnms, pattern = "_", fixed = T)
raw_expnms = substr( vec_expnms, 1,(underline_indices -1))
raw_expnms=gsub("[^[:alnum:] ]", "", raw_expnms)
raw_expnms = toupper(raw_expnms)
raw_editnms = gsub("[^[:alnum:] ]", "", editnms)
raw_editnms = toupper(raw_editnms)
commons=intersect(raw_expnms, raw_editnms)
commons = match(commons, raw_expnms)
commons =  vec_expnms[commons]
point_mut_mycells = data.frame()
my_ccle = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv",header = T,stringsAsFactors = F, check.names = F)
myccle2 = as.character(my_ccle[2,c(4:ncol(my_ccle))])
for(i in 1:length(myccle2)){
  bool =  point_mut$Tumor_Sample_Barcode ==  myccle2[i]
  point_mut_mycells =rbind(point_mut_mycells ,point_mut[bool,])
}

my_genes = c("TP53", "APC","BRAF", "KRAS","AXIN1", "AXIN2","CTNNB1", "SMAD4")

point_mut_mycells_indices = match(point_mut_mycells$Hugo_Symbol,my_genes)
point_mut_mycells_mygenes = point_mut_mycells[!is.na(point_mut_mycells_indices),]

cells =  levels(factor(point_mut_mycells_mygenes$Tumor_Sample_Barcode))
muts_by_cell = list()
for(i in 1:length(cells)){
  bool = point_mut_mycells_mygenes$Tumor_Sample_Barcode == cells[i]
  muts_by_cell[[i]] = point_mut_mycells_mygenes[bool,]
}
#### keep only the ones whcih are either isDeleterious or Ttcgahptspot
for(i in 1:length(muts_by_cell)){
  temp = muts_by_cell[[i]]
  bool = temp$isDeleterious | temp$isTCGAhotspot
  muts_by_cell[[i]] = temp[bool,]
}

output_mutation = data.frame(row.names = my_genes)
for(i im 1:length(muts_by_cell)){
  temp = muts_by_cell[[i]]
  for(j in 1:length(my_genes)){
    temp2 = temp[my_genes[j],]
    apply(temp2, 1 , FUN = paste0(x$))
  }
}

test = as.data.frame(muts_by_cell)

point_mut_mycells$Genome_Change
point_mut_mycells$ExAC_AF
point_mut_mycells$SangerWES_AC
point_mut_mycells$Variant_Classification
point_mut_mycells$Codon_Change
point_mut_mycells$Protein_Change
