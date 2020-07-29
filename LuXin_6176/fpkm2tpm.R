rm(list=ls())

# set work dir
setwd('/home/guoqiang/Documents/bioinfo/CPGEA/LuXin_6176/')

# import datasets
df.norm = read.csv('FOXA1_normal_expr_log2FPKM.csv', header = T, stringsAsFactors = F, check.names = F)
df.muta = read.csv('FOXA1_mutation_expr_log2FPKM.csv', header = T, stringsAsFactors = F, check.names = F)

# dissect matrix
df.norm.mtx = as.matrix(df.norm[,-c(1:2)])
df.muta.mtx = as.matrix(df.muta[,-c(1:2)])

# define function for transformation
fpkmToTpm <- function(fpkm) {
  (2^fpkm) * 1e6 / sum(2^fpkm)
}

# apply transformation on matrix
df.norm.tfm = apply(df.norm.mtx, 2, fpkmToTpm)
df.muta.tfm = apply(df.muta.mtx, 2, fpkmToTpm)

# integrate transformation and save
df.norm.tpm = data.frame(df.norm[,c(1,2)], df.norm.tfm)
df.muta.tpm = data.frame(df.muta[,c(1,2)], df.muta.tfm)
write.csv(df.norm.tpm, file = 'FOXA1_normal_expr_TPM.csv', row.names = F)
write.csv(df.muta.tpm, file = 'FOXA1_mutation_expr_TPM.csv', row.names = F)
