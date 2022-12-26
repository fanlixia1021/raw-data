

library(rms)                
inputFile="input.txt"        


rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)

dd <- datadist(rt)
options(datadist="dd")
f <- cph(Surv(futime, fustat) ~ age+Gender+riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)), 
    lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.5, 0.3,0.1,0.01))  

pdf(file="Nomogram.pdf",height=7,width=8.5)
plot(nom)
dev.off()

#calibration curve
time=3  
f <- cph(Surv(futime, fustat) ~ age+Gender+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=28, B=1000) 
pdf(file="calibration.pdf",height=6,width=7)
plot(cal,xlab="Nomogram-Predicted Probability of 3-Year OS",ylab="Actual 3-Year OS(proportion)",col="red",sub=F)
dev.off()


library(caret)
library(ggDCA)
library(rms)
inputFile="input.txt"        
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)

dd <- datadist(rt)
options(datadist="dd")

m1 <- lrm(status~ANLN,train_data)
m2 <- lrm(status~CENPA+GPR182+BCO2,train_data)
m3 <- lrm(status~ANLN+CENPA+GPR182+BCO2,train_data)d_train <- dca(m1)
ggplot(d_train)
d_train <- dca(m1,m2,m3)
ggplot(d_train)
d_train <- dca(m1,m2,m3,
              model.names =c('ANLN','Other 3 genes','All 4 genes'))
ggplot(d_train)
m1 <- cph(Surv(time,status)~ANLN,rt)
m2 <- cph(Surv(time,status)~CENPA+GPR182+BCO2,rt)
m3 <- cph(Surv(futime, fustat)~age+Gender+riskScore,rt)
d_train <- dca(m3)