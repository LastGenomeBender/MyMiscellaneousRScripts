
install.packages("vcfR")
library(vcfR)
setwd("D:/Farid/AOG lab/transepigenomics/New folder")
a=read.vcfR("HCT116vsCOLO-678.vcf")
b=read.vcfR("HCT116vsHT55.vcf")
c=read.vcfR("HCT116vsSW1116.vcf")
d=read.vcfR("LS411NvsCOLO-678.vcf")
e=read.vcfR("LS411NvsHT55.vcf")
f=read.vcfR("LS411NvsSW1116.vcf")
g=read.vcfR("LS513vsCOLO-678.vcf")
h=read.vcfR("LS513vsHT55.vcf")
i=read.vcfR("LS513vsSW1116.vcf")

la=((a@fix[,"REF"]=="A") & (a@fix[,"ALT"]=="G"))
lb = ((b@fix[,"REF"]=="A") & (b@fix[,"ALT"]=="G"))
lc= ((c@fix[,"REF"]=="A") & (c@fix[,"ALT"]=="G"))
ld= ((d@fix[,"REF"]=="A") & (d@fix[,"ALT"]=="G"))
le= ((e@fix[,"REF"]=="A") & (e@fix[,"ALT"]=="G"))
lf= ((f@fix[,"REF"]=="A") & (f@fix[,"ALT"]=="G"))
lg= ((g@fix[,"REF"]=="A") & (g@fix[,"ALT"]=="G"))
lh= ((h@fix[,"REF"]=="A") & (h@fix[,"ALT"]=="G"))
li= ((i@fix[,"REF"]=="A") & (i@fix[,"ALT"]=="G"))

a@fix=a@fix[la,]
a@gt = a@gt[la,]
b@fix=b@fix[lb,]
b@gt = b@gt[lb,]
c@fix=c@fix[lc,]
c@gt = c@gt[lc,]
d@fix=d@fix[ld,]
d@gt = d@gt[ld,]
e@fix=e@fix[le,]
e@gt = e@gt[le,]
f@fix=f@fix[lf,]
f@gt = f@gt[lf,]
g@fix=g@fix[lg,]
g@gt = g@gt[lg,]
h@fix=h@fix[lh,]
h@gt = h@gt[lh,]
i@fix=i@fix[li,]
i@gt = i@gt[li,]


write.vcf(a, file="HCT116vsCOLO-678_filtered.vcf")
write.vcf(b, file="HCT116vsHT55_filtered.vcf")
write.vcf(c, file="HCT116vsSW1116_filtered.vcf")
write.vcf(d, file="LS411NvsCOLO-678_filtered.vcf")
write.vcf(e, file="LS411NvsHT55_filtered.vcf")
write.vcf(f, file="LS411NvsSW1116_filtered.vcf")
write.vcf(g, file="LS513vsCOLO-678_filtered.vcf")
write.vcf(h, file="LS513vsHT55_filtered.vcf")
write.vcf(i, file="LS513vsSW1116_filtered.vcf")
nrow(a@fix)
