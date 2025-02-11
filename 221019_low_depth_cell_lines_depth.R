sig_depth=colnames(read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv",row.names = 1,header = T))
all_depth=colnames(read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/EditingRatio_all_cell_lines.csv",row.names = 1,header = T))
nonsig=setdiff(all_depth,sig_depth)
names(all_depth) = seq(1,58)

named_nonsig = all_depth[all_depth%in%nonsig]
names(named_nonsig)

cnms=read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_analysis/A-to-I GIREMI Colon CCLE/cell_line_names.csv",header = T,row.names = 1,stringsAsFactors = F)
nonsig_cnms = cnms[names(named_nonsig),]
nonsig
nonsig_cnms
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/221019_nonsig_depth_from_other_databases")
write(nonsig_cnms,"names_of_celllines_with_low_depth.csv")
