# function requires BSgenome Package, and USCS hg19(!!! 665Mb) reference genome acquired via BSgenome Package
# The code will try to install these automatically, but problems may arise, most probably in Ubuntu OS. Consult me in that case and
# we will figure out mising configurations in your PC
# 
# You need to input the path to miRNA seed sequences (/for/example/seeds.csv), path to your selected output of giremi,
# and the path to output folder.
# 
# Fucntions outputs 3 files,
# 1) "raw_mirBind_before_editing.csv" ==>A matrix of UNEDITED sequences with miRNA names in columns, gene names in rows, and binding locations
# of miRNA for 8-,7-,6-mers as "chr:start-end|chr:start-end|chr:start-end", respectively. if that miRNA
# does not bind to the specified region you will see NA instead. If one or 2 of the "mers" does not bind
# you will see "NA|chr:start-end|chr:start-end".
# 
# 2)"raw_mirBind_after_editing.csv" ==> The Matrix of EDITED sequences, in the same format as 1)
# 
# 3)"gained_and_lost_miRNAs.csv"==> A matrix with genes in rows. First column will include the miRNAs that CANNOT bind after editing(lost),
# and second row will include miRNAs that CAN bind after editing.

mirBindGiremi = function(path_to_mirseeds, path_to_editings, path_to_output_folder){
  
  tryCatch({library(BiocManager)},error = function(e) {install.packages("BiocManager")})  
  library(BiocManager)
  tryCatch({library(BSgenome)},error = function(e) {BiocManager::install("BSgenome")})  
  library(BSgenome)
  tryCatch({library(BSgenome.Hsapiens.UCSC.hg19)},error = function(e) {getBSgenome(genome = "BSgenome.Hsapiens.UCSC.hg19",masked = F)})  
  library(BSgenome.Hsapiens.UCSC.hg19)
  setwd(path_to_output_folder)
  miRNA_seeds = read.csv(path_to_mirseeds, stringsAsFactors = F, header = T, row.names = 1,check.names = F)
  sig_edits = rownames(read.csv(path_to_editings, stringsAsFactors = F, header = T, row.names = 1,check.names = F))
  sig_edits_splt = strsplit(sig_edits, split = "|", fixed=T)
  getsq = function(x) {
    chr = paste0("chr" , x[1])
    start = as.numeric(x[2]) - 7
    end = as.numeric(x[2]) +7
    strand = x[3]
    if(x[4]=="TC" & strand == "+"){
      strand == "+"
    }
    else if(x[4]=="TC" & strand == "-"){
      strand == "-"
    }
    sequences = getSeq(Hsapiens, names = chr, start = start, end = end)
  }
  gene_seqs = lapply(sig_edits_splt, FUN = function(x)getsq(x))
  gene_seqs_edited = lapply(gene_seqs,FUN=function(x) replaceAt(x = x,at = 8,value = "G"))
  #### modify the input file
  new_seeds = data.frame() 
  for( i in 1:nrow(miRNA_seeds)){
    temp = miRNA_seeds[i, ]
    new_seeds = rbind(new_seeds,as.character(temp[c(4,5,8,10,12)]), as.character(temp[c(13,14,17,19,21)]),stringsAsFactors = F)
  }
  colnames(new_seeds)=colnames(temp[c(4,5,8,10,12)])
  new_seeds = new_seeds[rowSums(new_seeds=="")==0,]
  pdict = c()
  for(i in 1:nrow(new_seeds)){
    pdict = c(pdict, as.character(new_seeds[i,c(3,4,5)]))
  }
  pdict_set = RNAStringSet(pdict)
  pdict_set = complement(DNAStringSet(pdict_set))
  matchFun = function(pdict,gene_list){
    res = lapply(gene_list,FUN = function(x) matchPDict(pdict,x))
  }
  bef_edit = matchFun(pdict = pdict_set, gene_seqs)
  after_edit = matchFun(pdict = pdict_set,gene_seqs_edited)
  fin_data = data.frame()
  for ( i in 1:length(bef_edit)){
    temp = bef_edit[[i]]
    results_big = c()
    counter = 1
    results = c()
    for(j in 1:length(temp)){
      temp2 = temp[[j]]
      if(length(temp2)==0){
        results = c(results, NA)
      }
      else if(start(temp2)>8|end(temp2)<8){
        results = c(results, NA)
      }
      else{
        out = paste0("chr",sig_edits_splt[[i]][1],":",(as.numeric(sig_edits_splt[[i]][2])-7+start(temp2)-1),"-",(as.numeric(sig_edits_splt[[i]][2])-7+end(temp2)-1))
        results = c(results, out)
      }
      counter = counter+1
      if(counter>3){
        if(sum(is.na(results))!=3){
          results_big = c(results_big, paste0(results,collapse = "|"))
          results = c()
          counter = 1
        }
        else{
          results_big = c(results_big, NA)
          results = c()
          counter = 1
        }
      }
    }
    fin_data = rbind(fin_data, results_big,stringsAsFactors=F)
  }
  mir_bef_edit=fin_data[!duplicated(colnames(fin_data))]
  colnames(mir_bef_edit) = new_seeds$Mature1_ID
  mir_bef_edit=mir_bef_edit[!duplicated(colnames(mir_bef_edit))]
  rownames(mir_bef_edit) = sig_edits
  log_mir_bef_edit = !is.na(mir_bef_edit)
  bindingMirs = c()
  for(i in 1:nrow(mir_bef_edit)){
    aqa=paste0(colnames(mir_bef_edit[log_mir_bef_edit[i,]]),collapse = "|")
    if(aqa==""){
      bindingMirs=c(bindingMirs,"NA")
      next()
    }
    bindingMirs=c(bindingMirs, aqa)
    
  }
  mir_bef_edit$bindingMirs = bindingMirs
  fin_data = data.frame()
  for ( i in 1:length(after_edit)){
    temp = after_edit[[i]]
    results_big = c()
    counter = 1
    results = c()
    for(j in 1:length(temp)){
      temp2 = temp[[j]]
      if(length(temp2)==0){
        results = c(results, NA)
      }
      else if(start(temp2)>8|end(temp2)<8){
        results = c(results, NA)
      }
      else{
        out = paste0("chr",sig_edits_splt[[i]][1],":",(as.numeric(sig_edits_splt[[i]][2])-7+start(temp2)-1),"-",(as.numeric(sig_edits_splt[[i]][2])-7+end(temp2)-1))
        results = c(results, out)
      }
      counter = counter+1
      if(counter>3){
        if(sum(is.na(results))!=3){
          results_big = c(results_big, paste0(results,collapse = "|"))
          results = c()
          counter = 1
        }
        else{
          results_big = c(results_big, NA)
          results = c()
          counter = 1
        }
      }
    }
    fin_data = rbind(fin_data, results_big,stringsAsFactors=F)
  }
  mir_after_edit=fin_data[!duplicated(colnames(fin_data))]
  colnames(mir_after_edit) = new_seeds$Mature1_ID
  mir_after_edit=mir_after_edit[!duplicated(colnames(mir_after_edit))]
  rownames(mir_after_edit) = sig_edits
  log_mir_after_edit = !is.na(mir_after_edit)
  bindingMirs = c()
  for(i in 1:nrow(mir_after_edit)){
    aqa=paste0(colnames(mir_after_edit[log_mir_after_edit[i,]]),collapse = "|")
    if(aqa==""){
      bindingMirs=c(bindingMirs,"NA")
      next()
    }
    bindingMirs=c(bindingMirs, aqa)
  }
  mir_after_edit$bindingMirs = bindingMirs
  strsplit(mir_after_edit$bindingMirs,split="|",fixed=T)
  funn = function(i){paste0(setdiff(strsplit(mir_after_edit$bindingMirs,split="|",fixed=T)[[i]],strsplit(mir_bef_edit$bindingMirs,split="|",fixed=T)[[i]]),collapse = "|")}
  gained_mirs = unlist(lapply(1:64, FUN = function(x) funn(x)))
  funn2 = function(i){paste0(setdiff(strsplit(mir_bef_edit$bindingMirs,split="|",fixed=T)[[i]],strsplit(mir_after_edit$bindingMirs,split="|",fixed=T)[[i]]),collapse = "|")}
  lost_mirs =  unlist(lapply(1:64, FUN = function(x) funn2(x)))
  finally = as.data.frame(cbind(lost_mirs,gained_mirs))
  finally[finally== ""] = NA
  rownames(finally) = sig_edits
  return(finally)
  write.csv(finally, "gained_and_lost_miRNAs.csv")
  write.csv(mir_bef_edit, "raw_mirBind_before_editing.csv")
  write.csv(mir_bef_edit, "raw_mirBind_after_editing.csv")
  
}


mirBindGiremi(path_to_mirseeds = "",path_to_editings = "",path_to_output_folder = "")
