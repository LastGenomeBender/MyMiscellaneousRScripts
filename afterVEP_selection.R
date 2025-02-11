######
setwd("D:/Farid/AOG lab/transepigenomics")
x = read.delim("SW1116vsSW620_vep.txt", sep  ="\t", stringsAsFactors = F)
l = (substr(x$Consequence,1,8) == "missense")
x= x[l,]
x= x[x$Existing_variation=="-",]
t =c()
t = c(T)
for(i in 2: nrow(x)){
  if(x$SYMBOL[i]==x$SYMBOL[i-1]&x$Location[i]==x$Location[i-1]){
    j = F
    t= c(t,j)
  }
  else{
    j=T
    t= c(t,j)
  }
}
x = x[t,]
write.csv(x,file="SW1116vsSW620_selected.csv")


