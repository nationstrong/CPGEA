library(immunedeconv)
rm(list = ls())

# set work director
setwd('/home/guoqiang/Documents/bioinfo/CPGEA/')

# load CIBERSORT 
set_cibersort_binary("deconvolution/CIBERSORT/source/CIBERSORT.R")
set_cibersort_mat("deconvolution/CIBERSORT/source/LM22.txt")

# import dataset
df.norm.tpm = read.csv('LuXin_6176/FOXA1_normal_expr_TPM.csv',
	header = T, stringsAsFactors = F, check.names = F)
head(df.norm.tpm)

df.muta.tpm = read.csv('LuXin_6176/FOXA1_mutation_expr_TPM.csv',
	header = T, stringsAsFactors = F, check.names = F)
head(df.muta.tpm)


# dissect matrix
df.norm.mtx = as.matrix(df.norm.tpm[,-c(1,2)])
rownames(df.norm.mtx) = df.norm.tpm$SYMBOL
head(df.norm.mtx)

df.muta.mtx = as.matrix(df.muta.tpm[,-c(1,2)])
rownames(df.muta.mtx) = df.muta.tpm$SYMBOL
head(df.muta.mtx)

# immune deconvolution

# CIBERSORT analysis
nams = rownames(df.norm.mtx)
rownames(df.norm.mtx) = make.names(nams, unique=TRUE)

nams = rownames(df.muta.mtx)
rownames(df.muta.mtx) = make.names(nams, unique=TRUE)

CIBERSORT.norm = deconvolute(df.norm.mtx, "cibersort") 
CIBERSORT.muta = deconvolute(df.muta.mtx, "cibersort") 

# write into csv
write.csv(CIBERSORT.norm, file = 'deconvolution/CIBERSORT/CIBERSORT_norm.csv', row.names  = F)
write.csv(CIBERSORT.muta, file = 'deconvolution/CIBERSORT/CIBERSORT_muta.csv', row.names  = F)


# CIBERSROT in abs mode
CIBERSORT.abs.norm = deconvolute(df.norm.mtx, "cibersort_abs") 
CIBERSORT.abs.muta = deconvolute(df.muta.mtx, "cibersort_abs") 

# write into csv
write.csv(CIBERSORT.abs.norm, file = 'deconvolution/CIBERSORT/CIBERSORT_abs_norm.csv', row.names  = F)
write.csv(CIBERSORT.abs.muta, file = 'deconvolution/CIBERSORT/CIBERSORT_abs_muta.csv', row.names  = F)