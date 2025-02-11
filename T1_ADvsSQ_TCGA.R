t1_sq_cnames= rep("SQ", times=length(t1_squamous_rep2))
t1_ad_cnames =  rep("AD", times=length(t1_AD_rep2))
t1.ADvsSq_counts = data.frame(t1_AD_rep2,t1_squamous_rep2)
colnames(t1.ADvsSq_counts)= c(t1_ad_cnames,t1_sq_cnames)
group <- factor(colnames(t1.ADvsSq_counts))
table(group)
design <- model.matrix(~ 0 + group)
design
## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design

fit <- lmFit(t1.ADvsSq_counts,design)
names(fit)
dim(fit$coefficients)
cont.matrix <- makeContrasts(ADvsSq=AD - SQ,levels=design)
cont.matrix
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)
boxplot(t1.ADvsSq_counts[,1:50], xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(data.matrix(t1.ADvsSq_counts[,1:10])),col="blue")

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

limma.res_t1 <- topTable(fit.cont, sort.by="p",n=100)
write.csv(limma.res_t1, file="T1_ADvsSCC_100mostSig_TCGA_rep2.csv")
