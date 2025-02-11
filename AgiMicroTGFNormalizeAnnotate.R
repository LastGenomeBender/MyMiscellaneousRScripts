# AgiMicroRMANormalizeAnnotate(), Farid Ahadli.
# You need to download AgiMicroRna package of bioconductor if you have not done so.
# Code ouputs gives log2(norm(TGS)) matrix with probe name and their annottaion to mirBase v21. 
# ifhalf logical, decides how to deal with negative values. if true, converts everything <0.5 to 0.5, otherwise a small positive constant (offset) is used and quantity (abs(min(data))+offset) is added to each signal.
# Agidirectory string, where your raw data is. Please uzip them, and delete the zips themselves. 
# Only raw data text files should exist in the Agidirectory
# outputdirectory string, where the result is saved. CanNOT be the same as Agidirectory.
# outputfname string, name of the final output. Ex:"test.csv". YOUR ANNOTATION FILE(AgiMicroRna_annotation_file.csv) HAS TO BE IN TEST OUTPUT DIRECTORY.
# normalizationMeth string, can be one of "none", "quantile", and "scale". Determines the cross array validation method
# example: AgiMicroRMANormalizeAnnotate(backgroundcorr = F,Agidirectory = "/home/ahadli/Desktop/film/tempGSE55139", outputfname = "test.csv", outputdirectory = "/home/ahadli/Desktop/film")
Agidirectory = "/home/ahadli/Desktop/film/GSE100488_RAW"
ifhalf = T
outputdirectory = "/home/ahadli/Desktop/film/tempGSE100488"
outputfname = "quantile.csv"
normalizationMeth = "quantile"
AgiMicroTGSNormalizeAnnotate<-function(ifhalf, Agidirectory, outputdirectory, outputfname,normalizationMeth){
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
  # dd = read.maimages(files = table, path = Agidirectory, green.only = T,source="agilent",other.columns=list(IsGeneDetected="gIsGeneDetected",
  #                                                                                                           IsSaturated="gIsSaturated",
  #                                                                                                           IsFeatNonUnifOF="gIsFeatNonUnifOL",
  #                                                                                                           IsFeatPopnOL="gIsFeatPopnOL",
  #                                                                                                           BGKmd="gBGMedianSignal"),annotation =  c( "ControlType", "ProbeName","GeneName"),
  #                    verbose=TRUE) 
  dd = readMicroRnaAFE(table, verbose = T)
  dd.TGS = tgsMicroRna(dd, half = ifhalf, makePLOT = F)
  ddNORM = tgsNormalization(dd.TGS, NORMmethod = normalizationMeth, makePLOTpre = F, makePLOTpost = F, table, verbose = T)
  ddNORM = filterMicroRna(ddNORM, dd, control = T, verbose = T,targets = table,writeout = F)
  esetPre =  esetMicroRna(ddNORM, table, makePLOT = F, verbose = T)
  setwd(outputdirectory)
  writeEset(eset = esetPre, ddPROC = ddNORM, targets = table, verbose = F)
  matrix = read.table("ProcessedData.txt", sep = "\t",row.names=1, header = T, stringsAsFactors = F)
  annot = read.csv("AgiMicroRna_annotation_file.csv",row.names=1, header = T, stringsAsFactors = F)
  # matrix$newAnnot = matrix$GENE
  # for(i in 1:nrow(matrix)){
  #   new_a = annot[rownames(matrix)[i],1]
  #   new_a = as.character(new_a)
  #   matrix$newAnnot[i] = new_a
  # }
  write.csv(matrix,outputfname)
  file.remove("ProcessedData.txt")
  file.remove("TargetFile.txt")
  setwd(oldwd)
} 
AgiMicroTGSNormalizeAnnotate(ifhalf = T, Agidirectory = "/home/ahadli/Desktop/film/GSE86241_RAW", outputdirectory = "/home/ahadli/Desktop/film/tempGSE86241",outputfname = "gse86241_ifhalTquatile.csv",normalizationMeth = "quantile")
