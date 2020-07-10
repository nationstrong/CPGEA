library(immunedeconv)

rm(list = ls())

# set work director
setwd('/home/guoqiang/Documents/bioinfo/CPGEA/')

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

# TIMER analysis
indi.timer = rep('PRAD', ncol(df.norm.mtx))
nams = rownames(df.norm.mtx)
rownames(df.norm.mtx) = make.names(nams, unique=TRUE)

TIMER.norm = deconvolute(df.norm.mtx, 'timer', 
	indications = indi.timer)

indi.timer = rep('PRAD', ncol(df.muta.mtx))
nams = rownames(df.muta.mtx)
rownames(df.muta.mtx) = make.names(nams, unique=TRUE)

TIMER.muta = deconvolute(df.muta.mtx, 'timer', 
	indications = indi.timer)

write.csv(TIMER.norm, file = 'deconvolution/TIMER/TIMER_norm.csv')
write.csv(TIMER.muta, file = 'deconvolution/TIMER/TIMER_muta.csv')