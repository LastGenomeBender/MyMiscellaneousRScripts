# AgiMicroRMANormalizeAnnotate(), Farid Ahadli.
# You need to download AgiMicroRna package of bioconductor if you have not done so.
#Code ouputs gives log2(RMA) matrix with probe name and their annottaion to mirBase v21. 
#backgroundcorr logical, decides if to perform RMA normalization with background correction or not
#Agidirectory string, where your raw data is. Please uzip them, and delete the zips themselves. 
#Only raw data text files should exist in the Agidirectory
#outputdirectory string, where the result is saved. CanNOT be the same as Agidirectory.
#outputfname string, name of the final output. Ex:"test.csv". YOUR ANNOTATION FILE(AgiMicroRna_annotation_file.csv) HAS TO BE IN TEST OUTPUT DIRECTORY.
# example: AgiMicroRMANormalizeAnnotate(backgroundcorr = F,Agidirectory = "/home/ahadli/Desktop/film/tempGSE55139", outputfname = "test.csv", outputdirectory = "/home/ahadli/Desktop/film")
AgiMicroRMANormalizeAnnotate<-function(backgroundcorr, Agidirectory, outputdirectory, outputfname){
  oldwd = getwd()
  setwd(Agidirectory)
  a=list.files()
  b= seq(1,length(a))
  setwd(outputdirectory)
  sink("TargetFile.txt")
  sink()
  for(i in 0:length(a)){
    if(i==0){
      sink("TargetFile.txt", append = T)
      cat("FileName","\t","Treatment","\t","GErep","\t","Subject" , "\n") 
      sink()
    }
    else if(i == length(a)){
      sink("TargetFile.txt", append = T)
      cat(a[i],"\t",1,"\t",1,"\t",b[i])
      sink()
    } 
    else{
      sink("TargetFile.txt", append = T)
      cat(a[i],"\t",1,"\t",1,"\t",b[i],"\n")
      sink()
    }
  }
  library(AgiMicroRna)
  table = readTargets(infile = "TargetFile.txt", verbose = T)
  setwd(Agidirectory)
  dd = readMicroRnaAFE(table)
  dd.RMA = rmaMicroRna(dd, normalize = T, background = backgroundcorr)
  esetPre =  esetMicroRna(dd.RMA, table, makePLOT = F, verbose = T)
  setwd(outputdirectory)
  writeEset(eset = esetPre, ddPROC = dd.RMA, targets = table, verbose = F)
  matrix = read.table("ProcessedData.txt", sep = "\t",row.names=1, header = T, stringsAsFactors = F)
  annot = read.csv("AgiMicroRna_annotation_file.csv",row.names=1, header = T, stringsAsFactors = F)
  matrix$newAnnot = matrix$GENE
  for(i in 1:nrow(matrix)){
    new_a = annot[rownames(matrix)[i],1]
    new_a = as.character(new_a)
    matrix$newAnnot[i] = new_a
  }
  write.csv(matrix,outputfname)
  file.remove("ProcessedData.txt")
  file.remove("TargetFile.txt")
  setwd(oldwd)
} 
