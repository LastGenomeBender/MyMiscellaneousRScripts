install.packages("vcfR")
library(vcfR)
setwd("D:/Farid/AOG lab/transepigenomics/5FU-resistance")
a=read.vcfR("C2BBE1vsHCT15.vcf")
b=read.vcfR("C2BBE1vsLS123.vcf")
c=read.vcfR("C2BBE1vsSW620.vcf")
d=read.vcfR("SNUC1vsHCT15.vcf")
e=read.vcfR("SNUC1vsLS123.vcf")
f=read.vcfR("SNUC1vsSW620.vcf")
g=read.vcfR("SW1116vsHCT-15.vcf")
h=read.vcfR("SW1116vsLS123.vcf")
i=read.vcfR("SW1116vsSW620.vcf")

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

write.vcf(a, file="C2BBE1vsHCT15_filtered.vcf")
write.vcf(b, file="C2BBE1vsLS123_filtered.vcf")
write.vcf(c, file="C2BBE1vsSW620_filtered.vcf")
write.vcf(d, file="SNUC1vsHCT15_filtered.vcf")
write.vcf(e, file="SNUC1vsLS123_filtered.vcf")
write.vcf(f, file="SNUC1vsSW620_filtered.vcf")
write.vcf(g, file="SW1116vsHCT-15_filtered.vcf")
write.vcf(h, file="SW1116vsLS123_filtered.vcf")
write.vcf(i, file="SW1116vsSW620_filtered.vcf")
nrow(a@fix)
 