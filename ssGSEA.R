

library(GSVA)
library(limma)
library(GSEABase)

rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())


ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)




exp = read.delim("1.txt", row.names = 1, header = T)
gene_set<-read.csv("mmc3.csv")[, 1:2]
head(gene_set)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
gsva_matrix<- gsva(as.matrix(exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix1<- t(scale(t(gsva_matrix)))
gsva_matrix1[gsva_matrix1< -2] <- -2
gsva_matrix1[gsva_matrix1>2] <- 2
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1)
write.csv(nor_gsva_matrix1,"nor_gsva_matrix.csv")
