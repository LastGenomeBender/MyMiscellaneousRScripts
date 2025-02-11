RPMmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  x = read.table(file=filenames[1], header = T)
  x=x[,c(1,3)]
  ordered_x = x[,2]
  ordered_x = ordered_x[ordered_x!=0]
  n_x= min(ordered_x)
  neww =  x[,2]
  neww[neww == 0] = n_x
  x[,2] = neww
  y = read.table(file = filenames[2], header = T)
  y = y[,3]
  ordered_y = y[y!=0]
  n_y = min(ordered_y)
  neww =  y
  neww[neww == 0] = n_y
  y = neww
  x = cbind(x,y)
  for (i in 3:(length(filenames)-1)){
    x2 = read.table(file=filenames[i], header = T)
    x2=x2[,3]
    ordered_x2 = x2[x2!=0]
    n_x2=min(ordered_x2)
    neww =  x2
    neww[neww == 0] = n_x2
    x2 = neww
    x = cbind(x,x2)
  }
  manifest= read.delim(file=filenames[length(filenames)], header = T, stringsAsFactors = F)
  fnames = manifest[,2]
  rrn = x[,1]
  x = x[,-1]
  rownames(x)=rrn
  for(i in 1:(length(filenames)-1)) {
    fname = filenames[i]
    fname = substring(fname, (nchar(mypath)+2), nchar(fname))
    for(j in 1 : length(fnames)){
      if (fname == (fnames[j])){
        colnames(x)[i] = manifest[j,7]
      }
    }
  }
  return(x)
}

