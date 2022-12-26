library(survival)
library(survminer)
pFilter=0.001                                                        
rt=read.table("xl.txt",header=T,sep="\t",check.names=F,row.names=1)    

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 if(coxP<pFilter){
     sigGenes=c(sigGenes,i)
		 outTab=rbind(outTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
  }
}
write.table(outTab,file="tcgaUniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="tcgaUniSigExp.txt",sep="\t",row.names=F,quote=F)



rt <- read.table("tcgaUniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))


pdf(file="forest.pdf", width = 6,height = 4.5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))


xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()


library(glmnet)
library(survival)


rt=read.table("tcgaUniSigExp.txt",header=T,sep="\t",row.names=1)            

#rt$futime=rt$futime/365


x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000,alpha = 1)
pdf("lasso1.pdf")
plot(fit,xvar="lambda",label=T)
dev.off()
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000,alpha=1,nfolds=10)
pdf("lasso2.pdf")
plot(cvfit)
dev.off()


coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
lassoSummary=summary(coef)
geneCoef=cbind(Gene=lassoGene,Coef=actCoef,
                  id=lassoGene,
		              coef= lassoSummary$coef.int[,"exp(coef)"],
		              coef.95L=lassoSummary$coef.int[,"lower .95"],
		              coef.95H=lassoSummary$coef.int[,"upper .95"],
		              pvalue=lassoSummary$coef[,"Pr(>|z|)"])
		              
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)




trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="tcgaRisk.txt",sep="\t",quote=F,row.names=F)


rt=read.table("cs.txt",header=T,sep="\t",row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
rt$futime=rt$futime/365
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="geoRisk.txt",sep="\t",quote=F,row.names=F)





bioSurvival=function(inputFile=null,outFile=null){
		rt=read.table(inputFile,header=T,sep="\t")                   
	
		diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
		pValue=1-pchisq(diff$chisq,df=1)
		pValue=signif(pValue,4)
		pValue=format(pValue, scientific = TRUE)
		fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
		surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=TRUE,
		           pval=paste0("p=",pValue),
		           pval.size=4,
		           risk.table=TRUE,
		           legend.labs=c("High risk", "Low risk"),
		           legend.title="Risk",
		           xlab="Time(years)",
		           break.time.by = 1,
		           risk.table.title="",
		           palette=c("red", "blue"),
		           risk.table.height=.25)
		pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
		print(surPlot)
		dev.off()
}
bioSurvival(inputFile="tcgaRisk.txt",outFile="tcgaRisk.pdf")
bioSurvival(inputFile="geoRisk.txt",outFile="geoRisk.pdf")


library(pheatmap)

bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
		rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)   
		rt=rt[order(rt$riskScore),]                                           
		riskClass=rt[,"risk"]
		lowLength=length(riskClass[riskClass=="low"])
		highLength=length(riskClass[riskClass=="high"])
		line=rt[,"riskScore"]
		line[line>10]=10
		pdf(file=riskScoreFile,width = 10,height = 3.5)
		plot(line, type="p", pch=20,
		     xlab="Patients (increasing risk socre)", ylab="Risk score",
		     col=c(rep("green",lowLength),rep("red",highLength)) )
		abline(h=median(rt$riskScore),v=lowLength,lty=2)
		legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
		dev.off()
		
		
		color=as.vector(rt$fustat)
		color[color==1]="red"
		color[color==0]="green"
		pdf(file=survStatFile,width = 10,height = 3.5)
		plot(rt$futime, pch=19,
		     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
		     col=color)
		legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
		abline(v=lowLength,lty=2)
		dev.off()
		
		
		rt1=rt[c(3:(ncol(rt)-2))]
		rt1=log2(rt1+1)
		rt1=t(rt1)
		annotation=data.frame(type=rt[,ncol(rt)])
		rownames(annotation)=rownames(rt)
		pdf(file=heatmapFile,width = 10,height = 3.5)
		pheatmap(rt1, 
		         annotation=annotation, 
		         cluster_cols = FALSE,
		         fontsize_row=11,
		         show_colnames = F,
		         fontsize_col=3,
		         color = colorRampPalette(c("green", "black", "red"))(50) )
		dev.off()
}
bioRiskPlot(inputFile="geoRisk.txt",riskScoreFile="geo.riskScore.pdf",survStatFile="geo.survStat.pdf",heatmapFile="geo.heatmap.pdf")
bioRiskPlot(inputFile="tcgaRisk.txt",riskScoreFile="tcga.riskScore.pdf",survStatFile="tcga.survStat.pdf",heatmapFile="tcga.heatmap.pdf")




library(survivalROC)

bioROC=function(riskFile=null,cliFile=null,outFile=null){
		risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
		cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
		sameSample=intersect(row.names(cli),row.names(risk))
		risk=risk[sameSample,]
		cli=cli[sameSample,]
		rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
		rocCol=rainbow(ncol(rt)-2)
		aucText=c()
		
		
		pdf(file=outFile,width=6,height=6)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
		  xlab="False positive rate", ylab="True positive rate",
		  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
		abline(0,1)
		
	
		j=1
		for(i in colnames(rt[,3:(ncol(rt)-1)])){
			roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")
			j=j+1
			aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
			lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
		}
		legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
		dev.off()
}
bioROC(riskFile="tcgaRisk.txt",cliFile="xlc.txt",outFile="tcgaROC.pdf")
bioROC(riskFile="geoRisk.txt",cliFile="csc.txt",outFile="geoROC.pdf")


library(rms)
riskFile="tcgaRisk.txt"
cliFile="tcn.txt"
outFile="tcga.Nomogram.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
dd <- datadist(rt)
options(datadist="dd")
f <- cph(Surv(futime, fustat) ~  x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)), 
    lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
pdf(file=outFile,height=7.5,width=11)
plot(nom)
dev.off()