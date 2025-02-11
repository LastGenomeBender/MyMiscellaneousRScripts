### cluster according to significant editings.

# 1) get significant editings

setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/3UTR_significant_editings_expression")
files = list.files()
files = files[grepl(x = files,pattern = ".csv",fixed = T)]
sig_edits = c()
for(i in files){
  sig_edits = c(sig_edits,rownames(read.csv(i, header = T,row.names = 1,check.names = F,stringsAsFactors = F)))
}
sig_edits = sig_edits[grepl(x = sig_edits,pattern = "|",fixed = T)]
head(sig_edits)
sig_edits = substr(sig_edits,1,regexpr(text = sig_edits,pattern = "/",fixed = T)-1)
sig_edits = sig_edits[!(duplicated(sig_edits))]
edits = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
sig_edits = edits[sig_edits,]

edit_dist = dist(sig_edits,method = "e")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff")
heatmap(as.matrix(sig_edits))
write.csv(sig_edits,file = "significant_editing.csv")
sig_edits = read.csv("significant_editing.csv",header = T,row.names = 1)
ss_edits = sig_edits
plot(hclust(edit_dist))
a = hclust(edit_dist)
sig_edit = data.frame()
for (i in 1:nrow(sig_edits)){
  std = sd(as.numeric(sig_edits[i,]))
  avg = mean(as.numeric(sig_edits[i,]))
  new = (as.numeric(sig_edits[i,]) - avg)/std
  sig_edit = rbind(sig_edit,new)
}
cor_mat=cor(sig_edits)
cor_mat2 = cor(t(sig_edits))
colnames(sig_edit) = colnames(sig_edits)
heatmap(cor_mat2)
rownames(sig_edit) = rownames(sig_edits)
library(RColorBrewer)
coul = rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(100))
heatmap(cor_mat2,verbose = T,margins = c(15,15))
  heatmap(as.matrix(sig_edit),col = coul,verbose = T,margins = c(15,15),scale = "row")
plot(heat.colors(2))
install.packages("RColorBrewer") 

library(gplots)

heatmap.2(as.matrix(cor_mat),col = coul,margins = c(20,20))
##### INTRONS
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Intron_significant_editings_expression")
files = list.files()
files = files[grepl(x = files,pattern = ".csv",fixed = T)]
sig_edits = c()
for(i in files){
  sig_edits = c(sig_edits,rownames(read.csv(i, header = T,row.names = 1,check.names = F,stringsAsFactors = F)))
}
sig_edits = sig_edits[grepl(x = sig_edits,pattern = "|",fixed = T)]
head(sig_edits)
sig_edits = substr(sig_edits,1,regexpr(text = sig_edits,pattern = "/",fixed = T)-1)
sig_edits = sig_edits[!(duplicated(sig_edits))]
edits = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
sig_edits = edits[sig_edits,]
std_sig_edits = apply(sig_edits, 1, FUN = function(x) )
edit_dist = dist(sig_edits,method = "e")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff")
heatmap(as.matrix(edits))
write.csv(sig_edits,file = "Introns_significant_editing.csv")
ss_edits = sig_edits
plot(hclust(edit_dist))
a = hclust(edit_dist)

heatmap(as.matrix(sig_edits),col = heat.colors(10000),verbose = T,margins = c(15,15))
plot(heat.colors(2))

###coding

setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/Coding_significant_editings_expression")
files = list.files()
files = files[grepl(x = files,pattern = ".csv",fixed = T)]
sig_edits = c()
for(i in files){
  sig_edits = c(sig_edits,rownames(read.csv(i, header = T,row.names = 1,check.names = F,stringsAsFactors = F)))
}
sig_edits = sig_edits[grepl(x = sig_edits,pattern = "|",fixed = T)]
head(sig_edits)
sig_edits = substr(sig_edits,1,regexpr(text = sig_edits,pattern = "/",fixed = T)-1)
sig_edits = sig_edits[!(duplicated(sig_edits))]
edits = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
sig_edits = edits[sig_edits,]
std_sig_edits = apply(sig_edits, 1, FUN = function(x) )
edit_dist = dist(sig_edits,method = "e")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff")
heatmap(as.matrix(edits))
write.csv(sig_edits,file = "Codings_significant_editing.csv")
ss_edits = sig_edits
plot(hclust(edit_dist))
a = hclust(edit_dist)

heatmap(as.matrix(sig_edits),col = heat.colors(50),verbose = T)
plot(heat.colors(2))
###5UTR
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/5UTR_significant_editings_expression")
files = list.files()
files = files[grepl(x = files,pattern = ".csv",fixed = T)]
sig_edits = c()
for(i in files){
  sig_edits = c(sig_edits,rownames(read.csv(i, header = T,row.names = 1,check.names = F,stringsAsFactors = F)))
}
sig_edits = sig_edits[grepl(x = sig_edits,pattern = "|",fixed = T)]
head(sig_edits)
sig_edits = substr(sig_edits,1,regexpr(text = sig_edits,pattern = "/",fixed = T)-1)
sig_edits = sig_edits[!(duplicated(sig_edits))]
edits = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff/EditingRatio_all_cell_lines.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
sig_edits = edits[sig_edits,]
std_sig_edits = apply(sig_edits, 1, FUN = function(x) )
edit_dist = dist(sig_edits,method = "e")
setwd("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/Rna Editing GIREMI/SNV_AG_Depth_Cutoff")
heatmap(as.matrix(edits))
write.csv(sig_edits,file = "5UTRs_significant_editing.csv")
ss_edits = sig_edits
plot(hclust(edit_dist))
a = hclust(edit_dist)

heatmap(as.matrix(sig_edits),col = heat.colors(50),verbose = T)
plot(heat.colors(2))

test=kmeans(t(sig_edits),centers = 3)
bb1=test$cluster
bb2 = test$cluster
names(bb1) == names (bb2)
plot(test)
str(test)
bb3=cbind(bb1,bb2)
