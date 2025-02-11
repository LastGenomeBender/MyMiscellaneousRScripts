#This function takes a count matrix (data.frame) and returns a data.matrix with standardized values.
standardize <- function(counts){
  standardized = data.frame(row.names = rownames(counts))
  for(i in 1:(length(rownames(counts)))){
    row = 0
    for(k in 1:(length(counts))){
      row = c(row, counts[i,k])
    }
    row = row[-1]
    mean = mean(row)
    std.dev = sd(row)
    for(j in 1:length(colnames(counts))){
      standardized[i,j] = (counts[i,j] - mean)/std.dev
    }
  }
  colnames(standardized) = colnames(counts)
  return(standardized)
}

