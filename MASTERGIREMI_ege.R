setwd("/media/egedede/HD-B1/masterGIREMI/unsorted_bams") 
bams = list.files(path="/media/egedede/HD-B1/masterGIREMI/unsorted_bams", full.names=F)

for(j in 1:length(bams)){
i = bams[j]
cat("running for the file", i, "started at", as.character(Sys.time()),"\n")
name = substring(i,1,(nchar(i)-4))
cmd = paste0("samtools rmdup -S /media/egedede/HD-B1/masterGIREMI/unsorted_bams/" ,name,".bam '/media/egedede/HD-B1/masterGIREMI/sorted_rmdup_bams/",name,"_dup.bam'")
print(cmd)
a = system(cmd, wait = T, intern = T)
#bams now are preprocessed.time for pileup, vcf calling
cmd = paste0("samtools index /media/egedede/HD-B1/masterGIREMI/sorted_rmdup_bams/",name,"_dup.bam > /media/egedede/HD-B1/masterGIREMI/sorted_rmdup_bams/",name,"_dup.bam.bai")
print(cmd)
system(cmd, intern=T, wait = T)
cmd = paste0("samtools idxstats '/media/egedede/HD-B1/masterGIREMI/sorted_rmdup_bams/",name,"_dup.bam' | cut -f1 |time parallel -j12 'samtools mpileup -uf /media/egedede/HD-B1/masterGIREMI/HG19_Broad_variant.fasta " ," -r {} /media/egedede/HD-B1/masterGIREMI/sorted_rmdup_bams/",name,"_dup.bam  -o /media/egedede/\"1TBVOL1\"/masterGIREMI/pileup/",name,"-{}.pileup'") 
print(cmd)
a = system(cmd, intern=T, wait = T)

#######
setwd("/media/egedede/1TBVOL1/masterGIREMI/pileup")
fnames = system("ls", intern=T,wait=T)
names = substr(fnames,1,(nchar(fnames)-7))
for(i in 1:length(fnames)){
	j = fnames[i]
	nm = names[i]	
	cmd = paste0("/home/egedede/anaconda2/bin/bcftools call -c -v ",j, " | /home/egedede/anaconda2/bin/bcftools filter -i \"DP>= 5\" >  " , "'/media/egedede/1TBVOL1/masterGIREMI/bcftools_out/",nm,".vcf'")
	cat(cmd)
	system(cmd, intern=T, wait = TRUE)
}

#######

setwd("/media/egedede/1TBVOL1/masterGIREMI/bcftools_out")
fnames = list.files(path = "/media/egedede/1TBVOL1/masterGIREMI/bcftools_out", full.names=T)
print(fnames)
for(i in 1:length(fnames)){
   cmd <- paste0("gzip ", fnames[i])
   print(cmd)
   a <- system(cmd, intern=T, wait = TRUE)
   print(a)
}

#######
setwd("/media/egedede/1TBVOL1/masterGIREMI/bcftools_out")
fnames = system("ls", intern=T,wait=T)
mm = paste(fnames, collapse = ' ')
cmd = paste0("vcf-concat ", mm," > final.vcf")
print(cmd)
system(cmd,intern=T,wait=T)
cmd = paste0("/home/egedede/anaconda2/bin/bcftools sort final.vcf | /home/egedede/anaconda2/bin/bcftools norm -d both > /media/egedede/HD-B1/masterGIREMI/",name,"_first.vcf")
print(cmd)
system(cmd, intern = T, wait = T)

#######


setwd("/home/egedede/snpEff_latest_core/snpEff")
###snpEff fuctional annotation
cmd = paste0("java -Xmx4g -jar '/home/egedede/snpEff_latest_core/snpEff/snpEff.jar' hg19 '/media/egedede/HD-B1/masterGIREMI/",name,"_first.vcf' >"  , " '/media/egedede/HD-B1/masterGIREMI/",name,"_functionally_anno.vcf'")
print(cmd)
a = system(cmd, wait = T, intern = T)
### snpSift dbSNP annotation
cmd = paste0("java -jar '/home/egedede/snpEff_latest_core/snpEff/SnpSift.jar'  annotate  -id '/media/egedede/HD-B1/masterGIREMI/GRCH37/dbsnp_135.b37.vcf'  '/media/egedede/HD-B1/masterGIREMI/",name,"_functionally_anno.vcf'  >", " '/media/egedede/HD-B1/masterGIREMI/dbnsp_and_func_vcf/",name,"_dbsnp_n_fun_annotated.vcf'")
print(cmd)
a = system(cmd, wait = T, intern = T)
#### only dbsnp cvf generation
cmd = paste0("java -jar '/home/egedede/snpEff_latest_core/snpEff/SnpSift.jar' filter  \"(ID =~ 'rs' )\" " , " '/media/egedede/HD-B1/masterGIREMI/dbnsp_and_func_vcf/",name,"_dbsnp_n_fun_annotated.vcf' > ", " '/media/egedede/HD-B1/masterGIREMI/onlysnp_vcf/",name,"_onlysnp.vcf'")
cat(cmd)
a = system(cmd, wait = T, intern = T)

####### mark_smp procedurec and GIREMI
setwd("/home/egedede/Desktop/giremi-master")
cmd = paste0("python  '/home/egedede/Desktop/giremi-master/mark_snp.py' -s ", " '/media/egedede/HD-B1/masterGIREMI/onlysnp_vcf/",name,"_onlysnp.vcf'"  ," -i '/media/egedede/HD-B1/masterGIREMI/dbnsp_and_func_vcf/",name,"_dbsnp_n_fun_annotated.vcf'   -g ", "  '/media/egedede/HD-B1/masterGIREMI/ref_gen.txt' > '/media/egedede/HD-B1/masterGIREMI/giremi_input_list/",name,"_giremi_input.txt'")
print(cmd)
a = system(cmd, wait = T, intern = T)
cmd=paste0("'/home/egedede/Desktop/giremi-master/giremi'  -f '/media/egedede/HD-B1/masterGIREMI/HG19_Broad_variant.fasta'  -l '/media/egedede/HD-B1/masterGIREMI/giremi_input_list/",name,"_giremi_input.txt'  -o  '/media/egedede/HD-B1/masterGIREMI/giremi_outputs/",name, "_output.lst'" ," '/media/egedede/HD-B1/masterGIREMI/sorted_rmdup_bams/",name,"_dup.bam'")
print(cmd)
a = system(cmd, wait = T, intern = T)
###### functional annotation of GIREMI output
setwd("/media/egedede/HD-B1/masterGIREMI/giremi_outputs")
library(tibble)
library(vcfR)
gfilenm= paste0(name, "_output.lst.res")
gir=read.delim(gfilenm, sep = "\t")
vcfflnm = paste0("/media/egedede/HD-B1/masterGIREMI/dbnsp_and_func_vcf/",name,"_dbsnp_n_fun_annotated.vcf")
vcf = read.vcfR(vcfflnm)
a = getINFO(vcf)
b = getPOS(vcf)
c = getCHROM(vcf)
ala = cbind(a,b,c)
func_anno = rep(x="anno", times = nrow(gir))
gir=add_column(gir, func_anno, .after = "chr")
bool = gir$ifRNAE>0
i1 = seq(1,nrow(gir),1)
i1 = i1[bool]
for (i in i1){
  chr = gir$chr[i]
  pos = gir$coordinate[i]
  bool = (c==chr&b==pos)
  c_hat = ala[bool,]
  info = c_hat[1]
  anno=unlist(strsplit(strsplit(info, split = "ANN")[[1]][2], split=","))
  anno2 = strsplit(anno, split="\\|")
  anno4 = sapply(X = anno2, function(x) paste(x[2],x[4]))
  anno4 = anno4[!duplicated(anno4)]
  anno5 = paste(anno4, collapse = ",")
  gir$func_anno[i] = anno5
}
fffname = paste0(name,"_final_funcanno_usethis.csv" )
write.csv(gir, file = fffname)
######
print("rm /media/egedede/HD-B1/masterGIREMI/sorted_rmdup_bams/*")
system("rm /media/egedede/HD-B1/masterGIREMI/sorted_rmdup_bams/*", wait=T, intern=T)
#### delete unnecessary vcf files
cmd = paste0(" rm '/media/egedede/HD-B1/masterGIREMI/",name,"_functionally_anno.vcf'")
print(cmd)
a = system(cmd, wait = T, intern = T)
cmd = paste0(" rm '/media/egedede/HD-B1/masterGIREMI/",name,"_first.vcf'")
print(cmd)
a = system(cmd, wait = T, intern = T)
#######removals
system("sudo rm /media/egedede/HD-B1/masterGIREMI/onlysnp_vcf/*", wait= T, intern = T)
system("sudo rm  /media/egedede/\"1TBVOL1\"/masterGIREMI/bcftools_out/* ",wait= T, intern = T)
system("sudo rm /media/egedede/\"1TBVOL1\"/masterGIREMI/pileup/* ",wait= T, intern = T)
cat("running for the file", i, "finished at", as.character(Sys.time()),"\n")

}


