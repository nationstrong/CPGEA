rm(list=ls())

library(ggplot2)
library(reshape2)

# set work director
setwd('/afs/crc.nd.edu/user/y/ygu8/project/CPGEA/deconvolution/CIBERSORT')

# import datasets
df.norm = read.csv('CIBERSORT_abs_norm.csv', row.names = 1,
                   header = T, stringsAsFactors = F, check.names = F)
df.muta = read.csv('CIBERSORT_abs_muta.csv', row.names = 1,
                   header = T, stringsAsFactors = F, check.names = F)
df.norm = t(df.norm)
df.muta = t(df.muta)

# melt and merge datasets
df.norm.melt = melt(df.norm)
colnames(df.norm.melt) = c('patient','cell.type','percentage')
df.norm.melt$FOXA1 = rep('normal', nrow(df.norm.melt))

df.muta.melt = melt(df.muta)
colnames(df.muta.melt) = c('patient','cell.type','percentage')
df.muta.melt$FOXA1 = rep('mutation', nrow(df.muta.melt))

df.melt = rbind(df.norm.melt, df.muta.melt)

# visualization of CIBERSORT
jpeg(filename = 'CIBERSORT_abs_boxplot.jpg', width = 10, height = 7, units = 'in', res = 300)
ggplot(data = df.melt, aes(x = cell.type, y = percentage, fill = FOXA1)) + geom_boxplot() + 
  theme_minimal() + 
  labs(title = 'CIBERSORT abs deconvolution of CPGEA cohort', 
       x = '', y = 'absolute level (AU)',
       fill = 'FOXA1') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


