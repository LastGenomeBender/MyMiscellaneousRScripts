setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/11.11.2019_AHR_edits_BAM_readcount/bam_readcount_outs")
names = list.files(full.names = F,include.dirs = F,pattern = ".txt")
temp=read.table(names[1],sep = "\t",stringsAsFactors = F,check.names = F)
rnnames = paste0(temp$V1,"|",temp$V2)
final_out = data.frame(row.names = rnnames)
for(i in 1:length(names)){
  bamread_count = read.table(names[i],sep = "\t",stringsAsFactors = F,check.names = F)

  chck = bamread_count$V3=="A"
  chck2 = bamread_count$V3[!is.na(bamread_count$V3)]=="T"
  chck3 = sum(chck2|chck)
  #if(chck3==nrow(bamread_count)){
    rownames(bamread_count) = paste0(bamread_count$V1, "|",bamread_count$V2)
    bamread_count = bamread_count[rnnames,]
    rownames(bamread_count) = rnnames
    g_counts= c()#vapply(X = strsplit(x = bamread_count$V8,fixed = T,split = ":"), FUN = function(X) X[2],FUN.VALUE = "1")
    for( j in 1:nrow(bamread_count)){
      if(bamread_count$V3[j]=="A"&(!is.na(bamread_count$V3[j]))){
        g_counts = c(g_counts,strsplit(x = bamread_count[j,c("V8")],fixed = T,split = ":")[[1]][2])
      }
      else if(bamread_count$V3[j]=="T"&(!is.na(bamread_count$V3[j]))){
        g_counts = c(g_counts,strsplit(x = bamread_count[j,c("V7")],fixed = T,split = ":")[[1]][2])
      }
      else{
        g_counts = c(g_counts,NA)
      }
    }
    out = paste0(g_counts,"/",bamread_count$V4) 
    final_out= cbind(final_out,out)
  #}
  #else{ 
  #}
}
out_names = gsub("_dup.txt","_final_funcanno_usethis.csv",names,fixed = T)
colnames(final_out) = out_names
write.csv(final_out,"/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/11.11.2019_AHR_edits_BAM_readcount/AHR_editing_analysis/AHR_g_total_reads.csv")

################# check if the 0 editings in AHR is indeed NA or 0 reads :)

editings=read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/significant_editing.csv",header = T,row.names = 1, check.names = F,stringsAsFactors = F)
editing_AHR = editings[grepl(pattern = "|AHR|",x = rownames(editings),fixed = T),] 
zero_edits=colnames(editings[,editing_AHR==0])
AHR_read = final_out[15,zero_edits]

out_names

zero_edits
