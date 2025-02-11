#### test the CYPIA1 levels between edit and non edit types.

editings = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/significant_editing.csv", header = T, row.names = 1, check.names = F, stringsAsFactors = F)
rpkm_data = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/CCLE_colon_54cell_expression_RPKM.csv",header = T, row.names = 1, check.names = F, stringsAsFactors = F) 

#### there are 33 cells with the expression value while 36 cells with editings. lets intersect them
genes = strsplit(x = rownames(editings),split = "|",fixed = T)
genes2 = unlist(lapply(genes, function(X) X[5])) 
new_mat = data.frame()
rownames(rpkm_data) = rpkm_data$V1
genenames = rpkm_data$V2
genenames = genenames[c(-1,-2)]
my_cell_names = as.vector(as.character(rpkm_data[1,]))
my_cell_names = my_cell_names[c(-1,-2)]
ccle_cell_names = as.vector(as.character(rpkm_data[2,]))
ccle_cell_names = ccle_cell_names[c(-1,-2)]
rpkm_data = rpkm_data[c(-1,-2)]
rpkm_data = rpkm_data[c(-1,-2),]
names(genenames)  = rownames(rpkm_data)
diff = setdiff(colnames(editings),colnames(rpkm_data))
AHR_editings_rppa_selected = editings[grepl(pattern ="7|17384354",x =  rownames(editings),fixed = T),]
diff = setdiff(colnames(editings),colnames(rpkm_data))
AHR_editings_rppa_selected = AHR_editings_rppa_selected[colnames(rpkm_data)]
AHR_noneditings = colnames(AHR_editings_rppa_selected[,AHR_editings_rppa_selected==0])
AHR_editings = setdiff(colnames(AHR_editings_rppa_selected),y = AHR_noneditings)
for(i in 1:ncol(rpkm_data)){
  rpkm_data[,i] = as.numeric(rpkm_data[,i])
}
rpkm_data = log2(rpkm_data+1)
CYP1A1=rpkm_data[genenames=="CYP1A1",]

boxplot(as.vector(as.numeric(CYP1A1[AHR_editings])),as.vector(as.numeric(CYP1A1[AHR_noneditings])) )
hist(as.numeric(CYP1A1),breaks = 30)
stripchart(x=list(as.numeric(CYP1A1[AHR_editings]),as.numeric(CYP1A1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(CYP1A1[AHR_editings]),as.numeric(CYP1A1[AHR_noneditings]))
wilcox.test(as.numeric(CYP1A1[AHR_editings]),as.numeric(CYP1A1[AHR_noneditings]))


# ALDH1 stemness marker
ALDH1=rpkm_data[genenames=="ALDH1A1",]
boxplot(as.vector(as.numeric(ALDH1[AHR_editings])),as.vector(as.numeric(ALDH1[AHR_noneditings])) )
hist(as.numeric(ALDH1),breaks = 30)
stripchart(x=list(as.numeric(ALDH1[AHR_editings]),as.numeric(ALDH1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(ALDH1[AHR_editings]),as.numeric(ALDH1[AHR_noneditings]))
wilcox.test(as.numeric(ALDH1[AHR_editings]),as.numeric(ALDH1[AHR_noneditings]))

########### ABCG2 stemness marker
ABCG2=rpkm_data[genenames=="ABCG2",]
boxplot(as.vector(as.numeric(ABCG2[AHR_editings])),as.vector(as.numeric(ABCG2[AHR_noneditings])) )
hist(as.numeric(ABCG2),breaks = 30)
stripchart(x=list(as.numeric(ABCG2[AHR_editings]),as.numeric(ABCG2[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(ABCG2[AHR_editings]),as.numeric(ABCG2[AHR_noneditings]))
wilcox.test(as.numeric(ABCG2[AHR_editings]),as.numeric(ABCG2[AHR_noneditings]))

#########  KLF4  stemness marker
KLF4=rpkm_data[genenames=="KLF4",]
boxplot(as.vector(as.numeric(KLF4[AHR_editings])),as.vector(as.numeric(KLF4[AHR_noneditings])) )
hist(as.numeric(KLF4),breaks = 30)
stripchart(x=list(as.numeric(KLF4[AHR_editings]),as.numeric(KLF4[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(KLF4[AHR_editings]),as.numeric(KLF4[AHR_noneditings]))
wilcox.test(as.numeric(KLF4[AHR_editings]),as.numeric(KLF4[AHR_noneditings]))

######## CDH1 - E-cadherin, Ephitelial marker

CDH1=rpkm_data[genenames=="CDH1",]
boxplot(as.vector(as.numeric(CDH1[AHR_editings])),as.vector(as.numeric(CDH1[AHR_noneditings])) )
hist(as.numeric(CDH1),breaks = 30)
stripchart(x=list(as.numeric(CDH1[AHR_editings]),as.numeric(CDH1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(CDH1[AHR_editings]),as.numeric(CDH1[AHR_noneditings]))
wilcox.test(as.numeric(CDH1[AHR_editings]),as.numeric(CDH1[AHR_noneditings]))

########  check the AHR mRNA levels.
AHR=rpkm_data[genenames=="AHR",]
boxplot(as.vector(as.numeric(AHR[AHR_editings])),as.vector(as.numeric(AHR[AHR_noneditings])) )
hist(as.numeric(AHR),breaks = 30)
stripchart(x=list(as.numeric(AHR[AHR_editings]),as.numeric(AHR[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(AHR[AHR_editings]),as.numeric(AHR[AHR_noneditings]))
wilcox.test(as.numeric(AHR[AHR_editings]),as.numeric(AHR[AHR_noneditings]))

#### ADAR1 gene
ADAR=rpkm_data[genenames=="ADAR",]
boxplot(as.vector(as.numeric(ADAR[AHR_editings])),as.vector(as.numeric(ADAR[AHR_noneditings])) )
hist(as.numeric(ADAR),breaks = 30)
stripchart(x=list(as.numeric(ADAR[AHR_editings]),as.numeric(ADAR[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(ADAR[AHR_editings]),as.numeric(ADAR[AHR_noneditings]))
wilcox.test(as.numeric(ADAR[AHR_editings]),as.numeric(ADAR[AHR_noneditings]))
 #### ADAR2 gene
ADARB1=rpkm_data[genenames=="ADARB1",]
boxplot(as.vector(as.numeric(ADARB1[AHR_editings])),as.vector(as.numeric(ADARB1[AHR_noneditings])) )
hist(as.numeric(ADARB1),breaks = 30)
stripchart(x=list(as.numeric(ADARB1[AHR_editings]),as.numeric(ADARB1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(ADARB1[AHR_editings]),as.numeric(ADARB1[AHR_noneditings]))
wilcox.test(as.numeric(ADARB1[AHR_editings]),as.numeric(ADARB1[AHR_noneditings]))

##### ADAR3 gene
ADARB2=rpkm_data[genenames=="ADARB2",]
boxplot(as.vector(as.numeric(ADARB2[AHR_editings])),as.vector(as.numeric(ADARB2[AHR_noneditings])) )
hist(as.numeric(ADARB2),breaks = 30)
stripchart(x=list(as.numeric(ADARB2[AHR_editings]),as.numeric(ADARB2[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(ADARB2[AHR_editings]),as.numeric(ADARB2[AHR_noneditings]))
wilcox.test(as.numeric(ADARB2[AHR_editings]),as.numeric(ADARB2[AHR_noneditings]))

### ITGB1 gene *****Significant

ITGB1=rpkm_data[genenames=="ITGB1",]
boxplot(as.vector(as.numeric(ITGB1[AHR_editings])),as.vector(as.numeric(ITGB1[AHR_noneditings])) )
hist(as.numeric(ITGB1),breaks = 30)
stripchart(x=list(as.numeric(ITGB1[AHR_editings]),as.numeric(ITGB1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(ITGB1[AHR_editings]),as.numeric(ITGB1[AHR_noneditings]))
wilcox.test(as.numeric(ITGB1[AHR_editings]),as.numeric(ITGB1[AHR_noneditings]))


#OCT4 gene
POU5F1=rpkm_data[genenames=="POU5F1",]
boxplot(as.vector(as.numeric(POU5F1[AHR_editings])),as.vector(as.numeric(POU5F1[AHR_noneditings])) )
hist(as.numeric(POU5F1),breaks = 30)
stripchart(x=list(as.numeric(POU5F1[AHR_editings]),as.numeric(POU5F1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(POU5F1[AHR_editings]),as.numeric(POU5F1[AHR_noneditings]))
wilcox.test(as.numeric(POU5F1[AHR_editings]),as.numeric(POU5F1[AHR_noneditings]))

## nanog gene

NANOG=rpkm_data[genenames=="NANOG",]
boxplot(as.vector(as.numeric(NANOG[AHR_editings])),as.vector(as.numeric(NANOG[AHR_noneditings])) )
hist(as.numeric(NANOG),breaks = 30)
stripchart(x=list(as.numeric(NANOG[AHR_editings]),as.numeric(NANOG[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(NANOG[AHR_editings]),as.numeric(NANOG[AHR_noneditings]))
wilcox.test(as.numeric(NANOG[AHR_editings]),as.numeric(NANOG[AHR_noneditings]))

#####   CYP1A2
CYP1A2=rpkm_data[genenames=="CYP1A2",]

boxplot(as.vector(as.numeric(CYP1A2[AHR_editings])),as.vector(as.numeric(CYP1A2[AHR_noneditings])) )
hist(as.numeric(CYP1A2),breaks = 30)
stripchart(x=list(as.numeric(CYP1A2[AHR_editings]),as.numeric(CYP1A2[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(CYP1A2[AHR_editings]),as.numeric(CYP1A2[AHR_noneditings]))
wilcox.test(as.numeric(CYP1A2[AHR_editings]),as.numeric(CYP1A2[AHR_noneditings]))

#######
CYP1B1=rpkm_data[genenames=="CYP1B1",]

boxplot(as.vector(as.numeric(CYP1B1[AHR_editings])),as.vector(as.numeric(CYP1B1[AHR_noneditings])) )
hist(as.numeric(CYP1B1),breaks = 30)
stripchart(x=list(as.numeric(CYP1B1[AHR_editings]),as.numeric(CYP1B1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(CYP1B1[AHR_editings]),as.numeric(CYP1B1[AHR_noneditings]))
wilcox.test(as.numeric(CYP1B1[AHR_editings]),as.numeric(CYP1B1[AHR_noneditings]))

#####MMP2
MMP2=rpkm_data[genenames=="MMP2",]

boxplot(as.vector(as.numeric(MMP2[AHR_editings])),as.vector(as.numeric(MMP2[AHR_noneditings])) )
hist(as.numeric(MMP2),breaks = 30)
stripchart(x=list(as.numeric(MMP2[AHR_editings]),as.numeric(MMP2[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(MMP2[AHR_editings]),as.numeric(MMP2[AHR_noneditings]))
wilcox.test(as.numeric(MMP2[AHR_editings]),as.numeric(MMP2[AHR_noneditings]))

###### cathepsin B
CTSB=rpkm_data[genenames=="CTSB",]

boxplot(as.vector(as.numeric(CTSB[AHR_editings])),as.vector(as.numeric(CTSB[AHR_noneditings])) )
hist(as.numeric(CTSB),breaks = 30)
stripchart(x=list(as.numeric(CTSB[AHR_editings]),as.numeric(CTSB[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(CTSB[AHR_editings]),as.numeric(CTSB[AHR_noneditings]))
wilcox.test(as.numeric(CTSB[AHR_editings]),as.numeric(CTSB[AHR_noneditings]))

#############  Cathepsin D **** significant in wlicox, hid in non-edit
CTSD=rpkm_data[genenames=="CTSD",]

boxplot(as.vector(as.numeric(CTSD[AHR_editings])),as.vector(as.numeric(CTSD[AHR_noneditings])) )
hist(as.numeric(CTSD),breaks = 30)
stripchart(x=list(as.numeric(CTSD[AHR_editings]),as.numeric(CTSD[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(CTSD[AHR_editings]),as.numeric(CTSD[AHR_noneditings]))
wilcox.test(as.numeric(CTSD[AHR_editings]),as.numeric(CTSD[AHR_noneditings]))



#### TGFB3
TGFB3=rpkm_data[genenames=="TGFB3",]

boxplot(as.vector(as.numeric(TGFB3[AHR_editings])),as.vector(as.numeric(TGFB3[AHR_noneditings])) )
hist(as.numeric(TGFB3),breaks = 30)
stripchart(x=list(as.numeric(TGFB3[AHR_editings]),as.numeric(TGFB3[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(TGFB3[AHR_editings]),as.numeric(TGFB3[AHR_noneditings]))
wilcox.test(as.numeric(TGFB3[AHR_editings]),as.numeric(TGFB3[AHR_noneditings]))

############
TGFB1=rpkm_data[genenames=="TGFB1",]

boxplot(as.vector(as.numeric(TGFB1[AHR_editings])),as.vector(as.numeric(TGFB1[AHR_noneditings])) )
hist(as.numeric(TGFB1),breaks = 30)
stripchart(x=list(as.numeric(TGFB1[AHR_editings]),as.numeric(TGFB1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(TGFB1[AHR_editings]),as.numeric(TGFB1[AHR_noneditings]))
wilcox.test(as.numeric(TGFB1[AHR_editings]),as.numeric(TGFB1[AHR_noneditings]))
#########snail
SNAI1=rpkm_data[genenames=="SNAI1",]

boxplot(as.vector(as.numeric(SNAI1[AHR_editings])),as.vector(as.numeric(SNAI1[AHR_noneditings])) )
hist(as.numeric(SNAI1),breaks = 30)
stripchart(x=list(as.numeric(SNAI1[AHR_editings]),as.numeric(SNAI1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(SNAI1[AHR_editings]),as.numeric(SNAI1[AHR_noneditings]))
wilcox.test(as.numeric(SNAI1[AHR_editings]),as.numeric(SNAI1[AHR_noneditings]))

###slug
SNAI2=rpkm_data[genenames=="SNAI2",]

boxplot(as.vector(as.numeric(SNAI2[AHR_editings])),as.vector(as.numeric(SNAI2[AHR_noneditings])) )
hist(as.numeric(SNAI2),breaks = 30)
stripchart(x=list(as.numeric(SNAI2[AHR_editings]),as.numeric(SNAI2[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(SNAI2[AHR_editings]),as.numeric(SNAI2[AHR_noneditings]))
wilcox.test(as.numeric(SNAI2[AHR_editings]),as.numeric(SNAI2[AHR_noneditings]))

###
HIF1A=rpkm_data[genenames=="HIF1A",]

boxplot(as.vector(as.numeric(HIF1A[AHR_editings])),as.vector(as.numeric(HIF1A[AHR_noneditings])) )
hist(as.numeric(HIF1A),breaks = 30)
stripchart(x=list(as.numeric(HIF1A[AHR_editings]),as.numeric(HIF1A[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(HIF1A[AHR_editings]),as.numeric(HIF1A[AHR_noneditings]))
wilcox.test(as.numeric(HIF1A[AHR_editings]),as.numeric(HIF1A[AHR_noneditings]))
###
NFATC1=rpkm_data[genenames=="NFATC1",]

boxplot(as.vector(as.numeric(NFATC1[AHR_editings])),as.vector(as.numeric(NFATC1[AHR_noneditings])) )
hist(as.numeric(NFATC1),breaks = 30)
stripchart(x=list(as.numeric(NFATC1[AHR_editings]),as.numeric(NFATC1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(NFATC1[AHR_editings]),as.numeric(NFATC1[AHR_noneditings]))
wilcox.test(as.numeric(NFATC1[AHR_editings]),as.numeric(NFATC1[AHR_noneditings]))
###
ENPP2=rpkm_data[genenames=="ENPP2",]

boxplot(as.vector(as.numeric(ENPP2[AHR_editings])),as.vector(as.numeric(ENPP2[AHR_noneditings])) )
hist(as.numeric(ENPP2),breaks = 30)
stripchart(x=list(as.numeric(ENPP2[AHR_editings]),as.numeric(ENPP2[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(ENPP2[AHR_editings]),as.numeric(ENPP2[AHR_noneditings]))
wilcox.test(as.numeric(ENPP2[AHR_editings]),as.numeric(ENPP2[AHR_noneditings]))

#####
TNF=rpkm_data[genenames=="TNF",]

boxplot(as.vector(as.numeric(TNF[AHR_editings])),as.vector(as.numeric(TNF[AHR_noneditings])) )
hist(as.numeric(TNF),breaks = 30)
stripchart(x=list(as.numeric(TNF[AHR_editings]),as.numeric(TNF[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(TNF[AHR_editings]),as.numeric(TNF[AHR_noneditings]))
wilcox.test(as.numeric(TNF[AHR_editings]),as.numeric(TNF[AHR_noneditings]))
####
JUN=rpkm_data[genenames=="JUN",]

boxplot(as.vector(as.numeric(JUN[AHR_editings])),as.vector(as.numeric(JUN[AHR_noneditings])) )
hist(as.numeric(JUN),breaks = 30)
stripchart(x=list(as.numeric(JUN[AHR_editings]),as.numeric(JUN[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(JUN[AHR_editings]),as.numeric(JUN[AHR_noneditings]))
wilcox.test(as.numeric(JUN[AHR_editings]),as.numeric(JUN[AHR_noneditings]))
#####HSPA4
HSPA4=rpkm_data[genenames=="HSPA4",]

boxplot(as.vector(as.numeric(HSPA4[AHR_editings])),as.vector(as.numeric(HSPA4[AHR_noneditings])) )
hist(as.numeric(HSPA4),breaks = 30)
stripchart(x=list(as.numeric(HSPA4[AHR_editings]),as.numeric(HSPA4[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(HSPA4[AHR_editings]),as.numeric(HSPA4[AHR_noneditings]))
wilcox.test(as.numeric(HSPA4[AHR_editings]),as.numeric(HSPA4[AHR_noneditings]))
######
IRF9=rpkm_data[genenames=="IRF9",]

boxplot(as.vector(as.numeric(IRF9[AHR_editings])),as.vector(as.numeric(IRF9[AHR_noneditings])) )
hist(as.numeric(IRF9),breaks = 30)
stripchart(x=list(as.numeric(IRF9[AHR_editings]),as.numeric(IRF9[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(IRF9[AHR_editings]),as.numeric(IRF9[AHR_noneditings]))
wilcox.test(as.numeric(IRF9[AHR_editings]),as.numeric(IRF9[AHR_noneditings]))
######

IL1B=rpkm_data[genenames=="IL1B",]

boxplot(as.vector(as.numeric(IL1B[AHR_editings])),as.vector(as.numeric(IL1B[AHR_noneditings])) )
hist(as.numeric(IL1B),breaks = 30)
stripchart(x=list(as.numeric(IL1B[AHR_editings]),as.numeric(IL1B[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(IL1B[AHR_editings]),as.numeric(IL1B[AHR_noneditings]))
wilcox.test(as.numeric(IL1B[AHR_editings]),as.numeric(IL1B[AHR_noneditings]))
###
JAK1=rpkm_data[genenames=="JAK1",]

boxplot(as.vector(as.numeric(JAK1[AHR_editings])),as.vector(as.numeric(JAK1[AHR_noneditings])) )
hist(as.numeric(JAK1),breaks = 30)
stripchart(x=list(as.numeric(JAK1[AHR_editings]),as.numeric(JAK1[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(JAK1[AHR_editings]),as.numeric(JAK1[AHR_noneditings]))
wilcox.test(as.numeric(JAK1[AHR_editings]),as.numeric(JAK1[AHR_noneditings]))
###
JAK2=rpkm_data[genenames=="JAK2",]

boxplot(as.vector(as.numeric(JAK2[AHR_editings])),as.vector(as.numeric(JAK2[AHR_noneditings])) )
hist(as.numeric(JAK2),breaks = 30)
stripchart(x=list(as.numeric(JAK2[AHR_editings]),as.numeric(JAK2[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(JAK2[AHR_editings]),as.numeric(JAK2[AHR_noneditings]))
wilcox.test(as.numeric(JAK2[AHR_editings]),as.numeric(JAK2[AHR_noneditings]))

####
JAK3=rpkm_data[genenames=="JAK3",]

boxplot(as.vector(as.numeric(JAK3[AHR_editings])),as.vector(as.numeric(JAK3[AHR_noneditings])) )
hist(as.numeric(JAK3),breaks = 30)
stripchart(x=list(as.numeric(JAK3[AHR_editings]),as.numeric(JAK3[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(JAK3[AHR_editings]),as.numeric(JAK3[AHR_noneditings]))
wilcox.test(as.numeric(JAK3[AHR_editings]),as.numeric(JAK3[AHR_noneditings]))
#####PPAR delta sig in both 

PPARD=rpkm_data[genenames=="PPARD",] 

boxplot(as.vector(as.numeric(PPARD[AHR_editings])),as.vector(as.numeric(PPARD[AHR_noneditings])) )
hist(as.numeric(PPARD),breaks = 30)
stripchart(x=list(as.numeric(PPARD[AHR_editings]),as.numeric(PPARD[AHR_noneditings])), method = "jitter",add = T,vertical = T)
?stripchart


t.test(as.numeric(PPARD[AHR_editings]),as.numeric(PPARD[AHR_noneditings]))
wilcox.test(as.numeric(PPARD[AHR_editings]),as.numeric(PPARD[AHR_noneditings]))

