##a-raf Ser-299 p positively regulates activity

install.packages("devtools")
library(devtools)
devtools::install_github("mariodeng/FirebrowseR")
library(FirebrowseR)
vignette("FirebrowseR")
?Samples.Clinical
coad_clinical = Samples.Clinical(format = "csv",cohort = "COAD")
library(BiocManager)
  install("TCGAbiolinks")
install.packages("devtools")
library(devtools)
devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')
library("TCGAbiolinks")
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "drug")
length(clinical$bcr_patient_barcode[!duplicated(clinical$bcr_patient_barcode)])
write.csv(clinical,"colon_drug_data.csv")
