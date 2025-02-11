t2_sq_cnames= rep("SQ", times=length(t2_squamous_rep2))
t2_ad_cnames =  rep("AD", times=length(t2_AD_rep2))
t2.ADvsSq_counts = data.frame(t2_AD_rep2,t2_squamous_rep2)
colnames(t2.ADvsSq_counts)= c(t2_ad_cnames,t2_sq_cnames)
group <- factor(colnames(t2.ADvsSq_counts))
table(group)
design <- model.matrix(~ 0 + group)
design
## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design

fit <- lmFit(t2.ADvsSq_counts,design)
names(fit)
dim(fit$coefficients)
cont.matrix <- makeContrasts(ADvsSq=AD - SQ,levels=design)
cont.matrix
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)
boxplot(t2.ADvsSq_counts[,1:50], xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(data.matrix(t2.ADvsSq_counts[,1:10])),col="blue")

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

limma.res_t2 <- topTable(fit.cont, sort.by="p",n=100)
write.csv(limma.res_t2, file="T2_ADvsSCC_100mostSig_TCGA_rep2.csv")

head(limma.res_t2)


