library("limma")
library("edgeR")
library(gplots)
library("org.Mm.eg.db")
library("RColorBrewer")
library("Glimma")
#set the working directory with setwd()
seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
sampleinfo <- read.delim("SampleInfo_Corrected.txt")
head(seqdata)
dim(seqdata)
sampleinfo
#only counts
countdata <- seqdata[,-(1:2)]
#assagin ids as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)
colnames(countdata)

# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
##Check if sample names in count is int same order with sampleindo
table(colnames(countdata)==sampleinfo$SampleName)

##When there are biological replicates in each group, in this case we have a sample size of 2 in each group, we favour filtering on a minimum counts per million threshold present in at least 2 samples

# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
head(keep)

dim(counts.keep)

#As a general rule, a good threshold can be chosen by identifying the CPM that corresponds to a count of 10, which in this case is about 0.5.
plot(myCPM[,1],countdata[,1])

y <- DGEList(counts.keep)
# have a look at y
y

# See what slots are stored in y
names(y)

# Library size information is stored in the samples slot
y$samples

y$samples$lib.size

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#multidemsional plot for Principle component
plotMDS(y)

# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$CellType)
## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)

# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
# Add a title
title("Cell type")

# Similarly for status
levels(sampleinfo$Status)

col.status <- c("blue","red","dark green")[sampleinfo$Status]
col.status

plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

# Dimension 3 appears to separate pregnant samples from the rest. Dim4?
#plotMDS(y,dim=c(3,4),col=col.status,pch=char.celltype,cex=2)
#legend("topright",legend=levels(sampleinfo$Status),col=cols,pch=16)
#legend("bottomright",legend=levels(sampleinfo$CellType),pch=c(1,4))

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")

# Save the heatmap
#png(file="High_var_genes.heatmap.png")
#heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
#dev.off()

#for interactive MDSplot
labels <- paste(sampleinfo$SampleName, sampleinfo$CellType, sampleinfo$Status)
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group <- factor(group)
glMDSPlot(y, labels=labels, groups=group, folder="mds")

# Apply normalisation(TMM norm) to DGEList object(for compositional )
y <- calcNormFactors(y)
#This will update the normalisation factors in the DGEList object (their default values are 1). Take a look at the normalisation factors for these samples.

y$samples

# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
design
## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
v

par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

# Fit the linear model (DiffGeneExp)
fit <- lmFit(v)
names(fit)

cont.matrix <- makeContrasts(B.PregVsLac=basal.pregnant - basal.lactate,L.PregVsLac=luminal.pregnant-luminal.lactate,levels=design)
cont.matrix
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)


summa.fit <- decideTests(fit.cont)
summary(summa.fit)

##for venn diagram
vennDiagram(summa.fit,include=c("up", "down"),counts.col=c("red", "blue"),circle.col = c("red", "blue", "green3"))

#summarize the DE
topTable(fit.cont,coef="B.PregVsLac",sort.by="p")

##ANNOTATION
columns(org.Mm.eg.db)
ann <- select(org.Mm.eg.db,keys=rownames(fit.cont),columns=c("ENTREZID","SYMBOL","GENENAME"))
head(ann)
table(ann$ENTREZID==rownames(fit.cont))

fit.cont$genes <- ann
topTable(fit.cont,coef="B.PregVsLac",sort.by="p")

##to get all differentially expresses genes
limma.res <- topTable(fit.cont,coef="B.PregVsLac",sort.by="p",n="Inf")

#to write into csv
write.csv(limma.res,file="B.PregVsLacResults.csv",row.names=FALSE)
fit.cont

##to discard genes below a threshold
# Let's decide that we are only interested in genes that have a absolute logFC of 1.
# This corresponds to a fold change of 2, or 0.5 (i.e. double or half).
# We can perform a treat analysis which ranks our genes according to p-value AND logFC.
# This is easy to do after our analysis, we just give the treat function the fit.cont object and specify our cut-off.
fit.treat <- treat(fit.cont,lfc=1)
res.treat <- decideTests(fit.treat)
summary(res.treat)
dim(res.treat)
