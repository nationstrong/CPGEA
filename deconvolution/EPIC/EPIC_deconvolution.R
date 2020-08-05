library(immunedeconv)

rm(list = ls())

# set work director
setwd('/afs/crc.nd.edu/user/y/ygu8/project/CPGEA')

# import dataset
df.norm.tpm = read.csv('LuXin_6176/FOXA1_normal_expr_TPM.csv',
	header = T, stringsAsFactors = F, check.names = F)
head(df.norm.tpm)

df.muta.tpm = read.csv('LuXin_6176/FOXA1_mutation_expr_TPM.csv',
	header = T, stringsAsFactors = F, check.names = F)
head(df.muta.tpm)


# dissect matrix
df.norm.mtx = as.matrix(df.norm.tpm[,-c(1,2)])
rownames(df.norm.mtx) = make.unique(df.norm.tpm$SYMBOL)
head(df.norm.mtx)

df.muta.mtx = as.matrix(df.muta.tpm[,-c(1,2)])
rownames(df.muta.mtx) = make.unique(df.muta.tpm$SYMBOL)
head(df.muta.mtx)

# immune deconvolution

# EPIC analysis

EPIC.norm = deconvolute(df.norm.mtx, 'epic')
EPIC.muta = deconvolute(df.muta.mtx, 'epic')

write.csv(EPIC.norm, file = 'deconvolution/EPIC/EPIC_norm.csv')
write.csv(EPIC.muta, file = 'deconvolution/EPIC/EPIC_muta.csv')