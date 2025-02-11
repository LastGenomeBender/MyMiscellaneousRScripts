###### TCGA patient ttest with the data of Han et.al
## read han etal

edit_tcga = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/[HAN et.al]A-to-I_TCGA/CRC_AG_annovar_nomutation.csv",header = F, stringsAsFactors = F, check.names = F,fill = T)
edit_tcga_trimmed = edit_tcga[1:229]
colnames(edit_tcga_trimmed) = edit_tcga_trimmed[1,]
rownames(edit_tcga_trimmed) = edit_tcga_trimmed$Samples
edit_tcga_trimmed = edit_tcga_trimmed[-1,-1]
rnms = rownames(edit_tcga_trimmed)
edit_tcga_trimmed = as.data.frame(sapply(edit_tcga_trimmed, FUN = as.numeric))
rownames(edit_tcga_trimmed) = rnms
AHR_tcga = edit_tcga_trimmed[c("chr7|17384354|UTR3|AHR|+|Alu|both"),]

###### read the patient data for COAD

drug_data = read.csv("colon_drug_data.csv",row.names = 1,header = T,stringsAsFactors = F, check.names = F)
drug_data = drug_data[drug_data$measure_of_response!="",]
drug_data = drug_data[grepl(pattern = "FU",x = drug_data$drug_name,fixed = T)|grepl(pattern = "Fluorouracil",x = drug_data$drug_name,useBytes = T,ignore.case = T),]
drug_data = drug_data[!duplicated(drug_data$bcr_patient_barcode),]
trimmed_tcga_barcode = substr(colnames(edit_tcga_trimmed),start = nchar("COAD-tumor--"),stop = nchar(colnames(edit_tcga_trimmed)))
indices = match(trimmed_tcga_barcode, drug_data$bcr_patient_barcode)
#indices = which(  drug_data$bcr_patient_barcode %in% trimmed_tcga_barcode)
AHR_tcga_samples_exists_in_both = AHR_tcga[!is.na(indices)] 
drug_data_exists_in_both = drug_data[indices[!is.na(indices)],]
drug_response = levels(factor(drug_data_exists_in_both$measure_of_response))
drug_response = drug_response[c(2,3,4,1)]
sub_response = list()
for(i in 1:length(drug_response)){
  sub_response[[i]] = drug_data_exists_in_both[drug_data_exists_in_both$measure_of_response==drug_response[i],]
}
##3test between responders and non-responders

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode,sub_response[[4]]$bcr_patient_barcode)],paired = F))



boxplot(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[4]]$bcr_patient_barcode)],paired = F))

barplot(as.numeric(AHR_tcga_samples_exists_in_both))


hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)

############### try another editing

edit_tcga = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/[HAN et.al]A-to-I_TCGA/CRC_AG_annovar_nomutation.csv",header = F, stringsAsFactors = F, check.names = F,fill = T)
edit_tcga_trimmed = edit_tcga[1:229]
colnames(edit_tcga_trimmed) = edit_tcga_trimmed[1,]
rownames(edit_tcga_trimmed) = edit_tcga_trimmed$Samples
edit_tcga_trimmed = edit_tcga_trimmed[-1,-1]
rnms = rownames(edit_tcga_trimmed)
edit_tcga_trimmed = as.data.frame(sapply(edit_tcga_trimmed, FUN = as.numeric))
rownames(edit_tcga_trimmed) = rnms
AHR_tcga = edit_tcga_trimmed[c("chr19|55900533|UTR3|RPL28|+|Alu|chimp"),]
AHR_tcga = AHR_tcga[,!is.na(AHR_tcga)]
AHR_tcga = AHR_tcga[,grepl(x=colnames(AHR_tcga),pattern = "COAD",fixed = T)]

hist(as.numeric(AHR_tcga))
###### read the patient data for COAD

drug_data = read.csv("colon_drug_data.csv",row.names = 1,header = T,stringsAsFactors = F, check.names = F)
drug_data = drug_data[drug_data$measure_of_response!="",]
trimmed_tcga_barcode = substr(colnames(AHR_tcga),start = nchar("COAD-tumor--"),stop = nchar(colnames(edit_tcga_trimmed)))
indices = match(trimmed_tcga_barcode, drug_data$bcr_patient_barcode)
AHR_tcga_samples_exists_in_both = AHR_tcga[,!is.na(indices)] 
drug_data_exists_in_both = drug_data[indices[!is.na(indices)],]
#drug_data = drug_data[grepl(pattern = "FU",x = drug_data$drug_name,fixed = T)|grepl(pattern = "Fluorouracil",x = drug_data$drug_name,useBytes = T,ignore.case = T),]
drug_data = drug_data_exists_in_both[!duplicated(drug_data$bcr_patient_barcode),]

drug_response = levels(factor(drug_data_exists_in_both$measure_of_response))
drug_response = drug_response[c(2,3,1)]
sub_response = list()

for(i in 1:length(drug_response)){
  sub_response[[i]] = drug_data_exists_in_both[drug_data_exists_in_both$measure_of_response==drug_response[i],]
}
##3test between responders and non-responders

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))



boxplot(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))
stripchart(list(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)])),vertical = T, pch=20,method = "jitter",add = T,col ="red")
barplot(as.numeric(AHR_tcga_samples_exists_in_both))


hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)


####  try another
edit_tcga = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/[HAN et.al]A-to-I_TCGA/CRC_AG_annovar_nomutation.csv",header = F, stringsAsFactors = F, check.names = F,fill = T)
edit_tcga_trimmed = edit_tcga[1:229]
colnames(edit_tcga_trimmed) = edit_tcga_trimmed[1,]
rownames(edit_tcga_trimmed) = edit_tcga_trimmed$Samples
edit_tcga_trimmed = edit_tcga_trimmed[-1,-1]
rnms = rownames(edit_tcga_trimmed)
edit_tcga_trimmed = as.data.frame(sapply(edit_tcga_trimmed, FUN = as.numeric))
rownames(edit_tcga_trimmed) = rnms
AHR_tcga = edit_tcga_trimmed[c("chr14|20835969|UTR3|TEP1|-|Alu|no-conserve"),]
AHR_tcga = AHR_tcga[,!is.na(AHR_tcga)]
#AHR_tcga = AHR_tcga[,grepl(x=colnames(AHR_tcga),pattern = "COAD",fixed = T)]

hist(as.numeric(AHR_tcga))
###### read the patient data for COAD

drug_data = read.csv("colon_drug_data.csv",row.names = 1,header = T,stringsAsFactors = F, check.names = F)
drug_data = drug_data[drug_data$measure_of_response!="",]
trimmed_tcga_barcode = substr(colnames(AHR_tcga),start = nchar("COAD-tumor--"),stop = nchar(colnames(edit_tcga_trimmed)))
indices = match(trimmed_tcga_barcode, drug_data$bcr_patient_barcode)
AHR_tcga_samples_exists_in_both = AHR_tcga[,!is.na(indices)] 
drug_data_exists_in_both = drug_data[indices[!is.na(indices)],]
#drug_data = drug_data[grepl(pattern = "FU",x = drug_data$drug_name,fixed = T)|grepl(pattern = "Fluorouracil",x = drug_data$drug_name,useBytes = T,ignore.case = T),]
drug_data = drug_data_exists_in_both[!duplicated(drug_data$bcr_patient_barcode),]

drug_response = levels(factor(drug_data_exists_in_both$measure_of_response))
drug_response = drug_response[c(2,3,1)]
sub_response = list()

for(i in 1:length(drug_response)){
  sub_response[[i]] = drug_data_exists_in_both[drug_data_exists_in_both$measure_of_response==drug_response[i],]
}
##3test between responders and non-responders

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))



boxplot(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))
stripchart(list(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)])),vertical = T, pch=20,method = "jitter",add = T,col ="red")
barplot(as.numeric(AHR_tcga_samples_exists_in_both))


hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)

####




edit_tcga = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/[HAN et.al]A-to-I_TCGA/CRC_AG_annovar_nomutation.csv",header = F, stringsAsFactors = F, check.names = F,fill = T)
edit_tcga_trimmed = edit_tcga[1:229]
colnames(edit_tcga_trimmed) = edit_tcga_trimmed[1,]
rownames(edit_tcga_trimmed) = edit_tcga_trimmed$Samples
edit_tcga_trimmed = edit_tcga_trimmed[-1,-1]
rnms = rownames(edit_tcga_trimmed)
edit_tcga_trimmed = as.data.frame(sapply(edit_tcga_trimmed, FUN = as.numeric))
rownames(edit_tcga_trimmed) = rnms
AHR_tcga = edit_tcga_trimmed[c("chr1|204524934|UTR3|MDM4|+|Alu|no-conserve"),]
AHR_tcga = AHR_tcga[,!is.na(AHR_tcga)]
#AHR_tcga = AHR_tcga[,grepl(x=colnames(AHR_tcga),pattern = "COAD",fixed = T)]

hist(as.numeric(AHR_tcga))
###### read the patient data for COAD

drug_data = read.csv("colon_drug_data.csv",row.names = 1,header = T,stringsAsFactors = F, check.names = F)
drug_data = drug_data[drug_data$measure_of_response!="",]
trimmed_tcga_barcode = substr(colnames(AHR_tcga),start = nchar("COAD-tumor--"),stop = nchar(colnames(edit_tcga_trimmed)))
indices = match(trimmed_tcga_barcode, drug_data$bcr_patient_barcode)
AHR_tcga_samples_exists_in_both = AHR_tcga[,!is.na(indices)] 
drug_data_exists_in_both = drug_data[indices[!is.na(indices)],]
#drug_data = drug_data[grepl(pattern = "FU",x = drug_data$drug_name,fixed = T)|grepl(pattern = "Fluorouracil",x = drug_data$drug_name,useBytes = T,ignore.case = T),]
drug_data_exists_in_both = drug_data_exists_in_both[!duplicated(drug_data_exists_in_both$bcr_patient_barcode),]

drug_response = levels(factor(drug_data_exists_in_both$measure_of_response))
drug_response = drug_response[c(2,3,4,1)]
sub_response = list()

for(i in 1:length(drug_response)){
  sub_response[[i]] = drug_data_exists_in_both[drug_data_exists_in_both$measure_of_response==drug_response[i],]
}
##3test between responders and non-responders

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))



boxplot(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))
stripchart(list(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)])),vertical = T, pch=20,method = "jitter",add = T,col ="red")
barplot(as.numeric(AHR_tcga_samples_exists_in_both))


hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)

################
edit_tcga = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/[HAN et.al]A-to-I_TCGA/CRC_AG_annovar_nomutation.csv",header = F, stringsAsFactors = F, check.names = F,fill = T)
edit_tcga_trimmed = edit_tcga[1:229]
colnames(edit_tcga_trimmed) = edit_tcga_trimmed[1,]
rownames(edit_tcga_trimmed) = edit_tcga_trimmed$Samples
edit_tcga_trimmed = edit_tcga_trimmed[-1,-1]
rnms = rownames(edit_tcga_trimmed)
edit_tcga_trimmed = as.data.frame(sapply(edit_tcga_trimmed, FUN = as.numeric))
rownames(edit_tcga_trimmed) = rnms
AHR_tcga = edit_tcga_trimmed[c("chr7|17384354|UTR3|AHR|+|Alu|both"),]
AHR_tcga = AHR_tcga[,!is.na(AHR_tcga)]
#AHR_tcga = AHR_tcga[,grepl(x=colnames(AHR_tcga),pattern = "COAD",fixed = T)]

hist(as.numeric(AHR_tcga))
###### read the patient data for COAD

drug_data = read.csv("colon_drug_data.csv",row.names = 1,header = T,stringsAsFactors = F, check.names = F)
drug_data = drug_data[drug_data$measure_of_response!="",]
trimmed_tcga_barcode = substr(colnames(AHR_tcga),start = nchar("COAD-tumor--"),stop = nchar(colnames(edit_tcga_trimmed)))
indices = match(trimmed_tcga_barcode, drug_data$bcr_patient_barcode)
AHR_tcga_samples_exists_in_both = AHR_tcga[,!is.na(indices)] 
drug_data_exists_in_both = drug_data[indices[!is.na(indices)],]
#drug_data = drug_data[grepl(pattern = "FU",x = drug_data$drug_name,fixed = T)|grepl(pattern = "Fluorouracil",x = drug_data$drug_name,useBytes = T,ignore.case = T),]
drug_data_exists_in_both = drug_data_exists_in_both[!duplicated(drug_data_exists_in_both$bcr_patient_barcode),]

drug_response = levels(factor(drug_data_exists_in_both$measure_of_response))
drug_response = drug_response[c(2,3,4,1)]
sub_response = list()

for(i in 1:length(drug_response)){
  sub_response[[i]] = drug_data_exists_in_both[drug_data_exists_in_both$measure_of_response==drug_response[i],]
}
##3test between responders and non-responders

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode)]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[4]]$bcr_patient_barcode)],paired = F))



boxplot(as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[1]]$bcr_patient_barcode,sub_response[[2]]$bcr_patient_barcode)]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[3]]$bcr_patient_barcode,sub_response[[4]]$bcr_patient_barcode)],paired = F))

barplot(as.numeric(AHR_tcga_samples_exists_in_both))


hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)


boxplot(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))
stripchart(list(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)])),vertical = T, pch=20,method = "jitter",add = T,col ="red")
barplot(as.numeric(AHR_tcga_samples_exists_in_both))


hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)


#######

edit_tcga = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/[HAN et.al]A-to-I_TCGA/CRC_AG_annovar_nomutation.csv",header = F, stringsAsFactors = F, check.names = F,fill = T)
edit_tcga_trimmed = edit_tcga[1:229]
colnames(edit_tcga_trimmed) = edit_tcga_trimmed[1,]
rownames(edit_tcga_trimmed) = edit_tcga_trimmed$Samples
edit_tcga_trimmed = edit_tcga_trimmed[-1,-1]
rnms = rownames(edit_tcga_trimmed)
edit_tcga_trimmed = as.data.frame(sapply(edit_tcga_trimmed, FUN = as.numeric))
rownames(edit_tcga_trimmed) = rnms
AHR_tcga = edit_tcga_trimmed[c("chr1|204524934|UTR3|MDM4|+|Alu|no-conserve"),]
AHR_tcga = AHR_tcga[,!is.na(AHR_tcga)]
#AHR_tcga = AHR_tcga[,grepl(x=colnames(AHR_tcga),pattern = "COAD",fixed = T)]

hist(as.numeric(AHR_tcga))
###### read the patient data for COAD

drug_data = read.csv("colon_drug_data.csv",row.names = 1,header = T,stringsAsFactors = F, check.names = F)
drug_data = drug_data[drug_data$measure_of_response!="",]
trimmed_tcga_barcode = substr(colnames(AHR_tcga),start = nchar("COAD-tumor--"),stop = nchar(colnames(edit_tcga_trimmed)))
indices = match(trimmed_tcga_barcode, drug_data$bcr_patient_barcode)
AHR_tcga_samples_exists_in_both = AHR_tcga[,!is.na(indices)] 
drug_data_exists_in_both = drug_data[indices[!is.na(indices)],]
#drug_data = drug_data[grepl(pattern = "FU",x = drug_data$drug_name,fixed = T)|grepl(pattern = "Fluorouracil",x = drug_data$drug_name,useBytes = T,ignore.case = T),]
drug_data_exists_in_both = drug_data_exists_in_both[!duplicated(drug_data_exists_in_both$bcr_patient_barcode),]

drug_response = levels(factor(drug_data_exists_in_both$measure_of_response))
drug_response = drug_response[c(2,3,4,1)]
sub_response = list()

for(i in 1:length(drug_response)){
  sub_response[[i]] = drug_data_exists_in_both[drug_data_exists_in_both$measure_of_response==drug_response[i],]
}
##3test between responders and non-responders

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))



boxplot(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))
stripchart(list(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)])),vertical = T, pch=20,method = "jitter",add = T,col ="red")
barplot(as.numeric(AHR_tcga_samples_exists_in_both))


hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)

################
edit_tcga = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/[HAN et.al]A-to-I_TCGA/CRC_AG_annovar_nomutation.csv",header = F, stringsAsFactors = F, check.names = F,fill = T)
edit_tcga_trimmed = edit_tcga[1:229]
colnames(edit_tcga_trimmed) = edit_tcga_trimmed[1,]
rownames(edit_tcga_trimmed) = edit_tcga_trimmed$Samples
edit_tcga_trimmed = edit_tcga_trimmed[-1,-1]
rnms = rownames(edit_tcga_trimmed)
edit_tcga_trimmed = as.data.frame(sapply(edit_tcga_trimmed, FUN = as.numeric))
rownames(edit_tcga_trimmed) = rnms
AHR_tcga = edit_tcga_trimmed[c("chr20|43706947|UTR3|STK4|+|Alu|both"),]
AHR_tcga = AHR_tcga[,!is.na(AHR_tcga)]
#AHR_tcga = AHR_tcga[,grepl(x=colnames(AHR_tcga),pattern = "COAD",fixed = T)]

hist(as.numeric(AHR_tcga))
###### read the patient data for COAD

drug_data = read.csv("colon_drug_data.csv",row.names = 1,header = T,stringsAsFactors = F, check.names = F)
drug_data = drug_data[drug_data$measure_of_response!="",]
trimmed_tcga_barcode = substr(colnames(AHR_tcga),start = nchar("COAD-tumor--"),stop = nchar(colnames(edit_tcga_trimmed)))
indices = match(trimmed_tcga_barcode, drug_data$bcr_patient_barcode)
AHR_tcga_samples_exists_in_both = AHR_tcga[,!is.na(indices)] 
drug_data_exists_in_both = drug_data[indices[!is.na(indices)],]
#drug_data = drug_data[grepl(pattern = "FU",x = drug_data$drug_name,fixed = T)|grepl(pattern = "Fluorouracil",x = drug_data$drug_name,useBytes = T,ignore.case = T),]
drug_data_exists_in_both = drug_data_exists_in_both[!duplicated(drug_data_exists_in_both$bcr_patient_barcode),]

drug_response = levels(factor(drug_data_exists_in_both$measure_of_response))
drug_response = drug_response[c(2,3,4,1)]
sub_response = list()

for(i in 1:length(drug_response)){
  sub_response[[i]] = drug_data_exists_in_both[drug_data_exists_in_both$measure_of_response==drug_response[i],]
}
##3test between responders and non-responders

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode,sub_response[[4]]$bcr_patient_barcode)],paired = F))

par(ncol)
for (i in 1:length(sub_response)) {

}

boxplot(as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[1]]$bcr_patient_barcode,sub_response[[2]]$bcr_patient_barcode)]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[3]]$bcr_patient_barcode,sub_response[[4]]$bcr_patient_barcode)],paired = F))

barplot(as.numeric(AHR_tcga_samples_exists_in_both))

boxplot(sub_response)
hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)


boxplot(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))
stripchart(list(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)])),vertical = T, pch=20,method = "jitter",add = T,col ="red")
barplot(as.numeric(AHR_tcga_samples_exists_in_both))


hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)




edit_tcga = read.csv("/media/ahadli/467DC90F6E3447F5/Ahadli/AOG_lab/[HAN et.al]A-to-I_TCGA/CRC_AG_annovar_nomutation.csv",header = F, stringsAsFactors = F, check.names = F,fill = T)
edit_tcga_trimmed = edit_tcga[1:229]
colnames(edit_tcga_trimmed) = edit_tcga_trimmed[1,]
rownames(edit_tcga_trimmed) = edit_tcga_trimmed$Samples
edit_tcga_trimmed = edit_tcga_trimmed[-1,-1]
rnms = rownames(edit_tcga_trimmed)
edit_tcga_trimmed = as.data.frame(sapply(edit_tcga_trimmed, FUN = as.numeric))
rownames(edit_tcga_trimmed) = rnms
AHR_tcga = edit_tcga_trimmed[c("chr7|17384354|UTR3|AHR|+|Alu|both"),]
AHR_tcga = AHR_tcga[,!is.na(AHR_tcga)]
#AHR_tcga = AHR_tcga[,grepl(x=colnames(AHR_tcga),pattern = "COAD",fixed = T)]

hist(as.numeric(AHR_tcga))
###### read the patient data for COAD

drug_data = read.csv("colon_drug_data.csv",row.names = 1,header = T,stringsAsFactors = F, check.names = F)
trimmed_tcga_barcode = substr(colnames(AHR_tcga),start = nchar("COAD-tumor--"),stop = nchar(colnames(edit_tcga_trimmed)))
indices = match(trimmed_tcga_barcode, drug_data$bcr_patient_barcode)
AHR_tcga_samples_exists_in_both = AHR_tcga[,!is.na(indices)] 
drug_data_exists_in_both = drug_data[indices[!is.na(indices)],]
hist(as.numeric(AHR_tcga_samples_exists_in_both))
#drug_data = drug_data[grepl(pattern = "FU",x = drug_data$drug_name,fixed = T)|grepl(pattern = "Fluorouracil",x = drug_data$drug_name,useBytes = T,ignore.case = T),]
drug_data_exists_in_both = drug_data_exists_in_both[!duplicated(drug_data_exists_in_both$bcr_patient_barcode),]

drug_response = levels(factor(drug_data_exists_in_both$measure_of_response))
drug_response = drug_response[c(2,3,4,1)]
sub_response = list()

for(i in 1:length(drug_response)){
  sub_response[[i]] = drug_data_exists_in_both[drug_data_exists_in_both$measure_of_response==drug_response[i],]
}
##3test between responders and non-responders

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))

colnames(AHR_tcga_samples_exists_in_both) = drug_data_exists_in_both$bcr_patient_barcode


wilcox.test(x = as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),y = as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode,sub_response[[4]]$bcr_patient_barcode)],paired = F))

par(ncol)
for (i in 1:length(sub_response)) {
  
}

boxplot(as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[1]]$bcr_patient_barcode,sub_response[[2]]$bcr_patient_barcode)]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[3]]$bcr_patient_barcode,sub_response[[4]]$bcr_patient_barcode)],paired = F))

barplot(as.numeric(AHR_tcga_samples_exists_in_both))

boxplot(sub_response)
hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)


boxplot(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)],paired = F))
stripchart(list(as.numeric(AHR_tcga_samples_exists_in_both[sub_response[[1]]$bcr_patient_barcode]),as.numeric(AHR_tcga_samples_exists_in_both[c(sub_response[[2]]$bcr_patient_barcode,sub_response[[3]]$bcr_patient_barcode)])),vertical = T, pch=20,method = "jitter",add = T,col ="red")
barplot(as.numeric(AHR_tcga_samples_exists_in_both))


hist(as.numeric(AHR_tcga_samples_exists_in_both),breaks = 4)
