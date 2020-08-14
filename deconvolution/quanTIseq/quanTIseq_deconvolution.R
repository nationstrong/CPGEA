library(immunedeconv)

rm(list = ls())

#set work directory
setwd('/afs/crc.nd.edu/user/y/ygu8/project/CPGEA')

#import dataset
df.norm.tpm <- read.csv('LuXin_6176/FOXA1_normal_expr_TPM.csv',
	header = T, stringsAsFactors = F, check.names = F)
head(df.norm.tpm[,1:5])

df.muta.tpm <- read.csv('LuXin_6176/FOXA1_mutation_expr_TPM.csv',
	header = T, stringsAsFactor = F, check.names = F)
head(df.muta.tpm[,1:5])

#dissect matrix
df.norm.mtx <- as.matrix(df.norm.tpm[,-c(1,2)])
row.names(df.norm.mtx) <- df.norm.tpm$SYMBOL
head(df.norm.mtx)

df.muta.mtx <- as.matrix(df.muta.tpm[,-c(1,2)])
row.names(df.muta.mtx) <- df.muta.tpm$SYMBOL
head(df.muta.mtx)

#immune deconvolution

#quanTIseq analysis

quanTIseq_norm <- deconvolute(df.norm.mtx,"quantiseq")
quanTIseq_muta <- deconvolute(df.muta.mtx,"quantiseq")

write.csv(quanTIseq_norm,file = "deconvolution/quanTIseq/quanTIseq_norm.csv")
write.csv(quanTIseq_muta,file = "deconvolution/quanTIseq/quanTIseq_muta.csv")

