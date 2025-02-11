check = c(TRUE)
t1rnames= rownames(limma.res_t1)
t2rnames = rownames(limma.res_t2)
for(i in 1:100){
  for(j in 1:100){
    if(t1rnames[i]==t2rnames[j]){
      if((limma.res_t1[i,1]<0 && limma.res_t2[j,1]<0) || (limma.res_t1[i,1]>0 && limma.res_t2[j,1]>0)){
        check= c(check, TRUE)
      }
      else{
        check = c(check, FALSE)
      }
    }
  }
}
check = check [-1]
AtSameDirection = sum(check) == 71