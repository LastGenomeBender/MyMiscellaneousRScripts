multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  x = read.table(file=filenames[1], header = T)
  x=x[,c(1,2)]
  y = read.table(file = filenames[2], header = T)
  y = y[,2]
  x = cbind(x,y)
  for (i in 3:(length(filenames)-1)){
    x2 = read.table(file=filenames[i], header = T)
    x2=x2[,2]
    x = cbind(x,x2)
  }
  manifest= read.delim(file=filenames[length(filenames)], header = T)
  fnames = manifest[,2]
  rrn = x[,1]
  x = x[,-1]
  rownames(x)=rrn
  for(i in 1:(length(filenames)-1)) {
    fname = filenames[i]
    fname = substring(fname, (nchar(mypath)+2), nchar(fname))
    for(j in 1 : length(fnames)){
      if (fname == (as.character(fnames[j])[1])){
        colnames(x)[i] = as.character(manifest[j,7])[1]
      }
    }
  }
  return(x)
  }

