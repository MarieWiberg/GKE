#------------------------------------------------
#  R Script for GKE book                       
#------------------------------------------------
# Chapter 9: Examples under the EG design
#------------------------------------------------

## 9.2.1 Preparing the SEPA Data for SNSequate
data("SEPA", package = "SNSequate")

library(equate)
SEPAx <- freqtab(x = SEPA$xscores, scales = 0:50)
SEPAy <- freqtab(x = SEPA$yscores, scales = 0:50)

plot(SEPAx)
plot(SEPAy)

rbind(x=summary(SEPAx),y=summary(SEPAy))

## 9.3 Step 1: Presmoothing

## 9.3.1 Beta4 Models
library(SNSequate) # not as code in the book
beta4nx <- BB.smooth(SEPAx,nparm=4,rel=0)
beta4ny <- BB.smooth(SEPAy,nparm=4,rel=0)

beta4nx

## 9.3.2 Discrete Kernel Estimators
psxB <- discrete.smooth(scores=rep(0:50,SEPAx),kert="bino",h=0.25,x=0:50)
psxD <- discrete.smooth(scores=rep(0:50,SEPAx),kert="dirDU",h=0,x=0:50)
psxT <- discrete.smooth(scores=rep(0:50,SEPAx),kert="triang",h=0.25,x=0:50)
                          
psxB

# 9.4 Step 2: Estimating Score Probabilities

## 9.4.1 Beta Models
rj.b4 <- beta4nx$prob.est
sk.b4 <- beta4ny$prob.est
cbind(0:50, rj.b4, sk.b4)

plot(0:50,as.matrix(SEPAx)/sum(as.matrix(SEPAx)),type="h",ylim=c(0,0.06),
     ylab="Relative Frequency",xlab="Scores")
lines(0:50,rj.b4,type="b",pch=15)

plot(0:50,as.matrix(SEPAy)/sum(as.matrix(SEPAy)),type="h",ylim=c(0,0.06),
     ylab="Relative Frequency",xlab="Scores")
lines(0:50,sk.b4,type="b",pch=15)

## 9.4.2 Discrete Kernel Estimators

psyB <- discrete.smooth(scores=rep(0:50,SEPAy),kert="bino",h=0.25,x=0:50) # not as code in the book

rj.dkb <- psxB$prob.est
sk.dkb <- psyB$prob.est # Note, psyB must be first created as above
cbind(0:50,rj.dkb,sk.dkb)

# In what follows some objects has to be created

rj.dkd <- psxD$prob.est # Not as code in the book
rj.dkt <- psxT$prob.est # Not as code in the book

plot(0:50,as.matrix(SEPAx)/sum(as.matrix(SEPAx)),lwd=2.0,xlab="Scores",
     ylab="Relative Frequency",type="h")
points(0:50,rj.dkb,type="b",pch=0)
points(0:50,rj.dkd,type="b",pch=1)
points(0:50,rj.dkt,type="b",pch=2)
legend("topright",pch=c(16,0,1,2),lty=c(1,1,1,1),
       c("Observed","Binomial","Dirac","Triangular"))

# All four objects below must be created, they are not given as code in the book
psyD <- discrete.smooth(scores=rep(0:50,SEPAy),kert="dirDU",h=0.00,x=0:50)
psyT <- discrete.smooth(scores=rep(0:50,SEPAy),kert="triang",h=0.25,x=0:50)
sk.dkd <- psyD$prob.est
sk.dkt <- psyT$prob.est

plot(0:50,as.matrix(SEPAy)/sum(as.matrix(SEPAy)),lwd=2.0,xlab="Scores",
     ylab="Relative Frequency",type="h")
points(0:50,sk.dkb,type="b",pch=0)
points(0:50,sk.dkd,type="b",pch=1)
points(0:50,sk.dkt,type="b",pch=2)
legend("topright",pch=c(16,0,1,2),lty=c(1,1,1,1),
       c("Observed","Binomial","Dirac","Triangular"))

## 9.5 Step 3: Continuization

## 9.5.1 Bandwidth Selection

rj.obs <- as.matrix(SEPAx)/sum(as.matrix(SEPAx))
sk.obs <- as.matrix(SEPAy)/sum(as.matrix(SEPAy))

hx.obs.ep <- bandwidth(scores=as.vector(SEPAx),kert="epan",
                       design="EG",r=rj.obs)$h
hx.b4.ep <- bandwidth(scores=as.vector(SEPAx),kert="epan",
                      design="EG",r=rj.b4)$h
hx.dkb.ep <- bandwidth(scores=as.vector(SEPAx),kert="epan",
                         design="EG",r=rj.dkb)$h

hx.obs.ad <- bandwidth(scores=as.vector(SEPAx),kert="adap",
                       design="EG",r=rj.obs)$h
hx.b4.ad <- bandwidth(scores=as.vector(SEPAx),kert="adap",
                      design="EG",r=rj.b4)$h
hx.dkb.ad <- bandwidth(scores=as.vector(SEPAx),kert="adap",
                       design="EG",r=rj.dkb)$h

hy.obs.ep <- bandwidth(scores=as.vector(SEPAy),kert="epan",
                         design="EG",r=sk.obs)$h
hy.b4.ep <- bandwidth(scores=as.vector(SEPAy),kert="epan",
                      design="EG",r=sk.b4)$h
hy.dkb.ep <- bandwidth(scores=as.vector(SEPAy),kert="epan",
                       design="EG",r=sk.dkb)$h

hy.obs.ad <- bandwidth(scores=as.vector(SEPAy),kert="adap",
                       design="EG",r=sk.obs)$h
hy.b4.ad <- bandwidth(scores=as.vector(SEPAy),kert="adap",
                      design="EG",r=sk.b4)$h
hy.dkb.ad <- bandwidth(scores=as.vector(SEPAy),kert="adap",
                       design="EG",r=sk.dkb)$h

# Pseudo-Table 9.1 
epa<-c(hx.obs.ep,hx.dkb.ep,hx.b4.ep,hy.obs.ep,hy.dkb.ep,hy.b4.ep)
ada<-c(hx.obs.ad,hx.dkb.ad,hx.b4.ad,hy.obs.ad,hy.dkb.ad,hy.b4.ad)
round(rbind(epa,ada),2)

## 9.6 Step 4: Equating
SEPAmat <- cbind(SEPAx,SEPAy)
eq.epan.b4 <- ker.eq(scores=SEPAmat,hx=hx.b4.ep,hy=hy.b4.ep,kert="epan", 
                     design="EG",r=rj.b4 ,s=sk.b4)
eq.epan.dkb <- ker.eq(scores=SEPAmat,hx=hx.dkb.ep,hy=hy.dkb.ep,kert="epan", 
                      design="EG",r=rj.dkb,s=sk.dkb)

comp.epan <- cbind(Scale=0:50,Epan.b4=eq.epan.b4$eqYx,Epan.dkb=eq.epan.dkb$eqYx)
round(comp.epan,4) 

# Below plot code is not as code in the book
plot(0:50, eq.epan.b4$eqYx, type="l", lwd=2.0,ylab=expression(varphi(x)),
     xlab=expression(x),col=1)
lines(0:50, eq.epan.dkb$eqYx, type="l", lwd=2.0,col=2)
legend("topleft",lty=c(1,1),col=c(1,2),c("Beta4","Binomial Kernel"))

eq.adap.b4 <- ker.eq(scores=SEPAmat,hx=hx.b4.ad,hy=hy.b4.ad,kert="adap",
                     design="EG",r=rj.b4,s=sk.b4,alpha=0.1)
comp.epan.adap <- cbind(Scale=0:50,Epan.b4=eq.epan.b4$eqYx,Adap.b4=eq.adap.b4$eqYx)
round(comp.epan.adap,4)

plot(0:50, eq.epan.b4$eqYx, type="l",ylab=expression(varphi(x)),
     xlab=expression(x),lwd=2.0,col=1)
lines(0:50, eq.adap.b4$eqYx, type="l", lwd=2.0,col=2)
legend("topleft",lty=c(1,1),col=c(1,2),c("Epanechnikov","Adaptive"))

## 9.7 Step 5: Evaluating the Equating Transformation

## 9.7.1 Bootstrap Standard Error of Equating

scores.x <- rep(0:50,SEPAx)
scores.y <- rep(0:50,SEPAy)

Nx <- length(scores.x)
Ny <- length(scores.y)
L <- 500

matboot <- matrix(NA,ncol=51,nrow=L)

set.seed(234)
i <- 0
it <- 0
while(it<L){
i <-i+1
x.b <-sample(scores.x,size=Nx,replace=TRUE)
y.b <-sample(scores.y,size=Ny,replace=TRUE)
SEPAx.b <-freqtab(x.b,scale=0:50)
SEPAy.b <-freqtab(y.b,scale=0:50)
SEPAmat.b<-cbind(SEPAx.b,SEPAy.b)

beta4nx.b<- BB.smooth(SEPAx.b,nparm=4,rel=0)
beta4ny.b<- BB.smooth(SEPAy.b,nparm=4,rel=0)
rj.b4.b <- beta4nx.b$prob.est
sk.b4.b <- beta4ny.b$prob.est

hx.b4.ep.b <- bandwidth(scores=as.vector(SEPAx.b),kert="epan", design="EG",
                        r=rj.b4.b)$h
hy.b4.ep.b <- bandwidth(scores=as.vector(SEPAy.b),kert="epan", design="EG",
                        r=sk.b4.b)$h

eq.epan.b4 <- tryCatch(ker.eq(scores=SEPAmat.b,hx=hx.b4.ep.b,hy=hy.b4.ep.b, 
                              kert="epan",design="EG",r=rj.b4.b,s=sk.b4.b),
                       error=function(e)NULL)
if(!is.null(eq.epan.b4)){
 it=it+1
 matboot[it,]<-eq.epan.b4$eqYx
 }
print(i)
print(it)
}
apply(matboot,2,sd)

# Below plot code is not as code in the book
plot(0:50,apply(matboot,2,sd),type="p",ylab="SEE",xlab="Scores")


## 9.7.2 PRE
data.frame(Epan.b4=round(PREp(eq.epan.b4,10)$preYx,4),
           Epan.dkb=round(PREp(eq.epan.dkb,10)$preYx,4),
           Adap.b4=round(PREp(eq.adap.b4,10)$preYx,4))


## 9.7.3 Freeman-Tukey Residuals
FTr.epan.b4 <-gof(SEPAx,eq.epan.b4$rj*eq.epan.b4$nx,"FT")
FTr.epan.dkb<-gof(SEPAx,eq.epan.dkb$rj*eq.epan.dkb$nx,"FT")
FTr.adap.b4 <-gof(SEPAx,eq.adap.b4$rj*eq.adap.b4$nx,"FT")

plot(0:50,FTr.epan.b4$ft.res,type="l",lty=1,ylab="Residuals",xlab="Scores",
     ylim=c(-3,3))  
lines(0:50,FTr.epan.dkb$ft.res,type="l",lty=2)
lines(0:50,FTr.adap.b4$ft.res,type="l",lty=3)
legend("topleft",lty=c(1,2,3),c("Epan.b4","Epan.dkb","Adap.b4"))
abline(h=0)
