#------------------------------------------------
#  R Script for GKE book                       
#------------------------------------------------
# Chapter 10: Examples under the NEAT design
#------------------------------------------------

## 10.2 ADM Data 
# The ADM data sets are available in the ADM file folder
# at the GitHub repository: https://github.com/MarieWiberg/GKE

## 10.2.1 Preparing the ADM Data for kequate

# Read in the data using one of two alternatives:
#-----------------------------------------------------------------
# 1. Load data directly from github. The links are not in the book.
load(url("https://github.com/MarieWiberg/GKE/raw/main/ADM/ADMX.Rda"))
load(url("https://github.com/MarieWiberg/GKE/raw/main/ADM/ADMY.Rda"))

# 2. Download the ADM data from https://github.com/MarieWiberg/GKE
# Set working directory (where you placed the ADM data)
setwd("C:/Users/....")

# Load data from your working directory
load("ADMX.Rda")
load("ADMY.Rda")
#-------------------------------------------------------------------

verb.xa <- apply(ADMX[,1:40],1,sum)
verb.ya <- apply(ADMY[,1:40],1,sum)
verb.x  <- apply(ADMX[,41:120],1,sum)
verb.y  <- apply(ADMY[,41:120],1,sum)

library(kequate) # not as code in the book

# Obtain bivariate score frequencies for kequate
neat.x <- kefreq(in1 = verb.x, xscores = 0:80,in2 = verb.xa, ascores = 0:40)
neat.y <- kefreq(in1 = verb.y, xscores = 0:80,in2 = verb.ya, ascores = 0:40)

neat.x[1000:1001,1:3]

library(equate) # Install equate if you do not already have it installed

# Obtain bivariate score frequencies for equate
feqX <- freqtab(cbind(verb.x,verb.xa),scales=list(0:80,0:40))
feqY <- freqtab(cbind(verb.y,verb.ya),scales=list(0:80,0:40))

plot(feqX,xlab="Test form X",ylab="Anchor test form A")
plot(feqY,xlab="Test form Y",ylab="Anchor test form A")

## 10.3 Simulated Polytomous Data

Xp <- vector("list", 55)
Yp <- vector("list", 55)
set.seed(7)
for(i in 1:55)
Xp[[i]] <- c(rnorm(1,-0.2,1),rnorm(1,0.4,1),runif(1))
set.seed(9)
for(i in 1:55)
Yp[[i]] <- c(rnorm(1,-0.2,1),rnorm(1,0.4,1),runif(1))

library(ltm)
XsP <- rmvordlogis(1000,Xp,model="gpcm",z.vals = rnorm(1000))
YsP <- rmvordlogis(1000,Yp,model="gpcm",z.vals = rnorm(1000))
colnames(XsP) <- paste("PX",1:55)
colnames(YsP) <- paste("PY",1:55)
XoP <- XsP-1
YoP <- YsP-1

## 10.4  Step 1: Presmoothing

## 10.4.1 Log-Linear Models

NEATvX <- glm(frequency~I(X)+I(X^2)+I(X^3)+I(A)+I(X):I(A),
              family = "poisson", data = neat.x, x = TRUE)
NEATvY <- glm(frequency~I(X)+I(X^2)+I(A)+I(X):I(A), 
              family = "poisson", data = neat.y, x = TRUE)

summary(NEATvX)

## 10.4.1.1 Modeling Complexities

neat.xm <- neat.x

neat.xm$frequency[neat.xm$X==0]  <- 10
neat.xm$frequency[neat.xm$X==35] <- 20
neat.xm$frequency[neat.xm$X==40] <- 20

neat.xm$ix1 <- numeric(length(neat.x$X))
neat.xm$ix2 <- numeric(length(neat.x$X))
neat.xm$ix1[neat.xm$X==0] <- 1
neat.xm$ix2[neat.xm$X %in% c(35,40)] <- 1

niX <- glm(frequency~I(X)+I(X^2)+I(X^3)+I(A)+I(X):I(A),data = neat.xm,
          family = "poisson", x = TRUE)

iX <- glm(frequency~I(X)+I(X^2)+I(X^3)+I(A)+ I(X):I(A)+I(ix1)+I(ix2),
          data = neat.xm,family = "poisson", x = TRUE)

summary(niX) # Not modeling spikes/gaps in the data
summary(iX)  # Modeling the spikes/gaps in the data


## 10.4.1.2 Log-Linear Model Fit

BIC(iX)

NEATvX.2 <- glm(frequency~I(X)+I(X^2)+I(X^3)+I(A)+I(A^2)
                  + I(X):I(A^2), family = "poisson",
                  data = neat.x, x = TRUE)
NEATvX.3 <- glm(frequency~I(X)+I(X^2)+I(A)+I(A^2)
                  + I(X):I(A^2), family = "poisson",
                  data = neat.x, x = TRUE)
NEATvX.4 <- glm(frequency~I(X)+I(X^2)+I(A)+I(A^2)
                  + I(X^2):I(A^2), family = "poisson",
                  data = neat.x, x = TRUE)
aic.com<-AIC(NEATvX,NEATvX.2,NEATvX.3,NEATvX.4)
bic.com<-BIC(NEATvX,NEATvX.2,NEATvX.3,NEATvX.4)
cbind(aic.com,bic.com)

obsNEATx <- matrix(neat.x$freq, nrow=81)/sum(NEATvX$y)
estNEATx <- matrix(NEATvX$fitted.values, nrow=81)/sum(NEATvX$y)
distNEATx <- cdist(est=estNEATx, obs=obsNEATx, 0:80, 0:40)
plot(distNEATx)

# Only the plot is shown in book - not the code
obsNEATx2 <- matrix(neat.x$freq, nrow=81)/sum(NEATvX.2$y)
estNEATx2 <- matrix(NEATvX.2$fitted.values, nrow=81)/sum(NEATvX.2$y)
distNEATx2 <- cdist(est=estNEATx2, obs=obsNEATx2, 0:80, 0:40)
plot(distNEATx2)

## 10.4.2 IRT Models with Binary-Scored Items
library(mirt)
ADMx.2PL<-mirt(data.frame(ADMX),model=1,itemtype="2PL",SE=TRUE)
ADMy.2PL<-mirt(data.frame(ADMY),model=1,itemtype="2PL",SE=TRUE)
head(coef(ADMx.2PL),1)

## 10.4.2.1 IRT Model Fit
admP1 <- as.matrix(ADMX[,c(41:120)])

library(psych) # not as code in the book

scree(admP1,factors=TRUE,pc=FALSE) 
vss(admP1,3) 

head(residuals(ADMx.2PL,"Q3",suppress=0.2),5) 

itemfit(ADMx.2PL) 

head(personfit(ADMx.2PL),2)

# Only plot shown in the book, not the below code
PF <- personfit(ADMx.2PL)
Zh <- PF$Zh
plot(Zh,xlab="Persons")

library(ltm)
admP <- as.matrix(ADMX[,c(41:120,1:40)])
adm2plP <- ltm(admP~z1, IRT.param=TRUE) 
plot(adm2plP,items = c(2,3),legend = TRUE,lty=c(1,2),col=c(1,1),labels=c("item 2","item 3"))

extract.mirt(ADMx.2PL,"logLik")
extract.mirt(ADMx.2PL,"AIC")
extract.mirt(ADMx.2PL,"BIC")

ADMx.R <- mirt(data.frame(ADMX),1,itemtype="Rasch",SE=TRUE)
anova(ADMx.R,ADMx.2PL)

## 10.4.3 IRT Models with Polytomously Scored Items
library(mirt)
X55r <- mirt(XoP,1,itemtype=rep('gpcm',55),SE=T) 
Y55r <- mirt(YoP,1,itemtype=rep('gpcm',55),SE=T)
 
coef(X55r,IRTpars=TRUE,simplify=TRUE)

## 10.4.3.1 IRT Model Fit
XoPm <- as.matrix(XoP[,c(1:40)])
scree(XoPm,factors=TRUE,pc=FALSE) 
vss(XoPm,3)
fa.parallel(XoPm,nfactors=1) 

fa(XoPm,nfactors=2)

X40r <- mirt(XoPm,1,itemtype=rep('gpcm',40),SE=T) 
residuals(X40r,"Q3",suppress=0.2)
 
itemfit(X55r)   
personfit(X55r) 

itemplot(X55r,3) 

## 10.5.1 Step 2: Estimating Score Probabilities

## 10.5.1 Log-Linear Models
### Object defined later in the book 
### glmPSE <- kequate("NEAT_PSE", 0:80, 0:80, NEATvX, NEATvY)
rj <- getScores(glmPSE)$X$r
sk <- getScores(glmPSE)$Y$s
cbind(0:80,rj,sk)

## 10.5.2 IRT Models
### Object defined later in the book
### m2pseSL <- irtose("PSE", ADMx.2PL, ADMy.2PL, 0:80, 0:80,0:40, model= "2pl", eqcoef = "Stocking-Lord")
rj_irt <- getScores(m2pseSL)$X$r
sk_irt <- getScores(m2pseSL)$Y$s
cbind(0:80,rj_irt,sk_irt)

# Marginal plts + presmoothing models (with objects defined later in the book)
library(equate)
marX<-margin(feqX, margin=1)
marY<-margin(feqY, margin=1)

plot(0:80,as.matrix(marX)/sum(as.matrix(marX)),lwd=2.0,
      xlab="Scores",ylab="Relative Frequency",type="h")
points(0:80,rj,type="b",pch=0)
points(0:80,rj_irt,type="b",pch=1)
legend("topright",pch=c(NA,0,1),lty=c(1,1,1),
      c("Observed","log-linear","IRT"))

plot(0:80,as.matrix(marY)/sum(as.matrix(marY)),lwd=2.0,
      xlab="Scores",ylab="Relative Frequency",type="h")
points(0:80,sk,type="b",pch=0)
points(0:80,sk_irt,type="b",pch=1)
legend("topright",pch=c(NA,0,1),lty=c(1,1,1),
      c("Observed","log-linear","IRT"))

## 10.7 Step 4: Equating 

## 10.7.1 Log-Linear Presmoothed Data
glmPSE <- kequate("NEAT_PSE", 0:80, 0:80, NEATvX, NEATvY)
summary(glmPSE)

glmCE <- kequate("NEAT_CE",0:80,0:80,0:40,NEATvX,NEATvY)
glmCEl <- kequate("NEAT_CE",0:80,0:80,0:40,NEATvX,NEATvY,kernel="logistic")
glmCEds <- kequate("NEAT_CE",0:80,0:80,0:40,NEATvX,NEATvY,kernel="uniform",
                   DS=1)
glmCEcv <- kequate("NEAT_CE",0:80,0:80,0:40,NEATvX,NEATvY,kernel="uniform", 
                   CV=1)
glmCEsrt <- kequate("NEAT_CE",0:80,0:80,0:40,NEATvX,NEATvY,kernel="uniform", 
                    altopt=TRUE)
summary(glmCE)

## 10.7.2 IRT Presmoothed Data
m2pseSL <- irtose("PSE", ADMx.2PL, ADMy.2PL, 0:80, 0:80,
                  + 0:40, model= "2pl", eqcoef = "Stocking-Lord")
summary(m2pseSL)

m2ce <- irtose("CE", ADMx.2PL, ADMy.2PL, 0:80, 0:80, 0:40)
summary(m2ce)

# Figure comparing KE and IRT KE transformations

llM  <- getEq(glmPSE) 
irtM <- getEq(m2pseSL)

plot(0:80,llM,lwd=2,xlab="x",type="l",lty=1,ylab=expression(varphi(x)))
lines(0:80,irtM,lty=2, lwd=2)
legend("topleft",lty=c(1,2),c("KE","IRT KE"))

# Figure comparing Gaussian, Logistic and Uniform CE
CEG  <- getEq(glmCE)
CEL  <- getEq(glmCEl) 
glmCEu <- kequate("NEAT_CE",0:80,0:80,0:40,NEATvX,NEATvY,
                 kernel="uniform")
CEU  <- getEq(glmCEu) 

plot(0:80,CEG,lty=1,xlab="x",type="l",ylab=expression(varphi(x)))
lines(0:80,CEL,lty=2, lwd=2)
lines(0:80,CEU,lty=3, lwd=3)
legend("topleft",lty=c(1,2,3),c("Gaussian","Logistic","Uniform"))

# IRT KE with polytomously scored item data
CEP <- irtose("CE", X55r,Y55r,0:80,0:80,0:30, catsX = rep(3,40),
              catsY = rep(3,40), catsA = rep(3,15),model = "GPCM")
summary(CEP) # not in the book


## 10.8 Step 5: Evaluating the Equating Transformation

## 10.8.1 PRE

PREglmPSE <- getPre(glmPSE) # output in Section 10.7.1

## 10.8.2 SEE
SeeADM <- getSee(glmCE)
SeeADM

plot(0:80, SeeADM, xlab="X scores", ylab ="SEE") 

## 10.8.3 Bootstrap Standard Error of Equating
library(ltm)
admP <- as.matrix(ADMX[,c(41:120,1:40)]) 
admQ <- as.matrix(ADMY[,c(41:120,1:40)])
ADMxltm.2PL <- ltm(admP~z1, IRT.param=TRUE)
ADMyltm.2PL <- ltm(admQ~z1, IRT.param=TRUE)

m2ceAN <- irtose("CE", ADMxltm.2PL, ADMyltm.2PL, 0:80, 0:80,
                 0:40, model= "2pl", see="analytical")
ANsee <- getSee(m2ceAN)
m2ceBoot<- irtose("CE", ADMxltm.2PL, ADMyltm.2PL, 0:80, 0:80,
                  0:40, model= "2pl", see="bootstrap", 
                  replications=250) # Takes time to run
                     
Bootsee <- getSee(m2ceBoot)

plot(ANsee,ylim=c(0,0.5),type='l',lty=1,ylab="SEE",xlab="Scores")
lines(Bootsee,col=4,lty=2)
legend(box.lty=0,"topright",inset=0.02,
       lty=c(1,2),col=c(1,4),c("Analytical","Bootstrap"))

# 10.8.4 SEED
SeedADM      <- getSeed(glmCE)
glmCEPSEseed <- genseed(glmCE, glmPSE)
glmCEPSEseed

plot(glmCEPSEseed)

## 10.8.5 MSD, MAD, and RMSD
n <- 81
MSD <- 1/n *sum(getEq(glmCE)-getEq(glmPSE))
MAD <- 1/n *sum(abs(getEq(glmCE)-getEq(glmPSE)))
RMSD <- sqrt(1/n*sum(((getEq(glmCE)-getEq(glmPSE))^2)))

MSD; MAD; RMSD 