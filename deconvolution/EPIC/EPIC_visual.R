rm(list=ls())

library(ggplot2)
library(reshape2)

# set work director
setwd('/afs/crc.nd.edu/user/y/ygu8/project/CPGEA/deconvolution/EPIC')

# import datasets
df.norm = read.csv('EPIC_norm.csv', row.names = 1,
                   header = T, stringsAsFactors = F, check.names = F)
df.muta = read.csv('EPIC_muta.csv', row.names = 1,
                   header = T, stringsAsFactors = F, check.names = F)

rownames(df.norm) = df.norm$cell_type
rownames(df.muta) = df.muta$cell_type
df.norm = df.norm[,c(-1)]
df.muta = df.muta[,c(-1)]

df.norm = t(df.norm)
head(df.norm)
df.muta = t(df.muta)
head(df.muta)


# melt and merge datasets
df.norm.melt = melt(df.norm)
colnames(df.norm.melt) = c('patient','cell.type','percentage')
df.norm.melt$FOXA1 = rep('normal', nrow(df.norm.melt))

df.muta.melt = melt(df.muta)
colnames(df.muta.melt) = c('patient','cell.type','percentage')
df.muta.melt$FOXA1 = rep('mutation', nrow(df.muta.melt))

df.melt = rbind(df.norm.melt, df.muta.melt)
head(df.melt)

# visualization of EPIC
jpeg(filename = 'EPIC_boxplot.jpg', width = 8, height = 6, units = 'in', res = 300)
ggplot(data = df.melt, aes(x = cell.type, y = percentage, fill = FOXA1)) + geom_boxplot() + 
  theme_minimal() + 
  labs(title = 'EPIC deconvolution of CPGEA cohort', 
       x = '', y = 'percentage',
       fill = 'FOXA1')+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()


