
t1_AD_primary_cnames= rep("ADpri", times=length(t1_AD_p_mi))
t1_AD_normal_cnames =  rep("ADnorm", times=length(t1_AD_n_mi))
t1.ADprivsADnorm = data.frame(t1_AD_p_mi,t1_AD_n_mi)
colnames(t1.ADprivsADnorm)= c(t1_AD_primary_cnames,t1_AD_normal_cnames)
filtered_t1.ADprivsADnorm = t1.ADprivsADnorm[rowSums(t1.ADprivsADnorm<=10)<10,]
y <- DGEList(filtered_t1.ADprivsADnorm)
group <- factor(colnames(filtered_t1.ADprivsADnorm))
table(group)
design <- model.matrix(~ 0 + group)
design
y <- calcNormFactors(y)
v <- voom(y,design,plot = TRUE)
v
fit <- lmFit(v)
names(fit)

boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

cont.matrix <- makeContrasts(PrimaryvsNormal= groupADpri- groupADnorm,levels=design)
cont.matrix
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)


summa.fit <- decideTests(fit.cont)
summary(summa.fit)

limma.res <- topTable(fit.cont,coef="PrimaryvsNormal",sort.by="p",n="100")
head(limma.res)
write.csv(limma.res,"T1_NormalvsTumor_AD_edgeR.csv")
