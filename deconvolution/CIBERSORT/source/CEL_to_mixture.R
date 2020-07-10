library(affy)
library(annotate)
library(hgu133plus2hsentrezgcdf)
library(org.Hs.eg.db)
Data<-ReadAffy(cdfname = "hgu133plus2hsentrezgcdf")
#eset<-rma(Data)
eset<-mas5(Data)
ID<-featureNames(eset)
ID2<-sub("_at","",ID)
GS <- as.matrix(getSYMBOL(ID2, 'org.Hs.eg'))
ematrix<-exprs(eset)
rows <- GS
cols =  c("GeneSymbol",colnames(ematrix))
ematrix <- cbind(rows,ematrix)
ematrix <- ematrix[which(ematrix[,1] != "NA"),] #remove NAs
ematrix <- ematrix[order(ematrix[,1]),] #sort by gene name 
ematrix <- rbind(cols, ematrix)
write.table(ematrix,file="NormalizedExpressionArray.customCDF.txt",sep="\t", col.names=F, row.names=F,quote=FALSE)
