rm(list = ls())

library(ggplot2)
library(reshape2)

#set work directory
setwd('/afs/crc.nd.edu/user/y/ygu8/project/CPGEA/deconvolution/MCP-counter')

#import datasets
df.norm <- read.csv('MCP-counter_norm.csv',
	header = T, stringsAsFactor = F, check.names = F)

df.muta <- read.csv('MCP-counter_muta.csv',
	header = T, stringsAsFactor = F, check.names = F)

rownames(df.norm) = df.norm$cell_type
rownames(df.muta) = df.muta$cell_type
df.norm <- df.norm[,-c(1,2)]
df.muta <- df.muta[,-c(1,2)]

df.norm <- t(df.norm)
df.muta <- t(df.muta)

#melt and merge datasets
df.norm.melt <- melt(df.norm)
colnames(df.norm.melt) <- c('patients','cell.type','gene.expression')
df.norm.melt$FOXA1 <- rep('normal',nrow(df.norm.melt))

df.muta.melt <- melt(df.muta)
colnames(df.muta.melt) <- c('patients','cell.type','gene.expression')
df.muta.melt$FOXA1 <- rep('mutant',nrow(df.muta.melt))

df.melt <- rbind(df.norm.melt,df.muta.melt)
head(df.melt)

#conduct t.test
p.list <- c()
for (i in c(1:ncol(df.norm))){
  tt <- t.test(df.norm[,i],df.muta[,i])
  tp <- round(tt$p.value,3)
  if(is.na(tp) == F){
   star = 'n.s'
    if(tp < 0.05){
     star = '*'
      if(tp < 0.01){
       star = '**'
        if(tp < 0.001){
         star = '***'
        }
      }
    }
  p.list <- c(p.list,paste0('p=',tp,'\n',star))
 }else{
 p.list <- c(p.list,"N.A")
 }
}

gene_names <- colnames(df.norm)
df.text <- data.frame(gene = gene_names, pval = p.list, val = rep(max(df.melt[,3])*1.1,ncol(df.norm)))

#Visualization of MCP-counter
jpeg(filename = 'MCP-counter_boxplot.jpg', width = 8, height = 6, units = 'in', res = 300)
ggplot(data = df.melt, aes(x = cell.type, y = gene.expression, fill = FOXA1)) + geom_boxplot() + 
  theme_minimal() + 
  geom_text(data = df.text,aes(x = gene, y = val, label = pval), inherit.aes = FALSE, size = 2.5) + 
  labs(title = 'MCP-counter deconvolution of CPGEA cohort', 
       x = '', y = 'gene expression(AU)',
       fill = 'FOXA1')+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()