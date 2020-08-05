library(immunedeconv)

rm(list = ls())

#set work directory
setwd('/afs/crc.nd.edu/user/y/ygu8/project/CPGEA')

#import dataset
df.norm.tpm <- read.csv('LuXin_6176/FOXA1_normal_expr_TPM.csv',
	header = T, stringsAsFactors = F, check.names = F)
head(df.norm.tpm[])

df.muta.tpm <- read.csv('LuXin_6176/FOXA1_mutation_expr_TPM.csv',
	header = T, stringsAsFactor = F, check.names = F)
head(df.muta.tpm)

#dissect matrix
df.norm.mtx <- as.matrix(df.norm.tpm[,-c(1,2)])
row.names(df.norm.mtx) <- df.norm.tpm$SYMBOL
head(df.norm.mtx)

df.muta.mtx <- as.matrix(df.muta.tpm[,-c(1,2)])
row.names(df.muta.mtx) <- df.muta.tpm$SYMBOL
head(df.muta.mtx)

#immune deconvolution

#MCP-counter analysis

MCPcounter_norm <- deconvolute(df.norm.mtx,"mcp_counter")
MCPcounter_muta <- deconvolute(df.muta.mtx,"mcp_counter")

write.csv(MCPcounter_norm,file = "deconvolution/MCP-counter/MCP-counter_norm.csv")
write.csv(MCPcounter_muta,file = "deconvolution/MCP-counter/MCP-counter_muta.csv")


