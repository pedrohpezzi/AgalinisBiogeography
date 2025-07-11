#script adapted from Valderrama et al. (2022) - https://doi.org/10.3389/fpls.2022.874322
setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/11_HiSSE")

#install.packages('hisse')
library(hisse)
library(diversitree)
library(geiger)
library(phytools)

setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/11_HiSSE")
tree<-read.tree(file='InputTrees/cluster5_nosynonyms_noout.tre')
plot(tree)
lifehistory<-read.csv(file='lifehistory.csv',h=T)
rownames(lifehistory)<-lifehistory$species

lh<-lifehistory$lifehistory
names(lh)<-lifehistory$species
lhist<-lh[names(lh)%in%tree$tip.label]
lhist<-as.factor(lhist)
x<-droplevels(lhist)
x<-x[order(factor(names(x),levels=tree$tip.label))]
x

hd<-cbind(rownames(x),x)
colnames(hd)<-'lifehistory'
hd<-data.frame(hd)
hd$lifehistory<-gsub(1,0,hd$lifehistory)
hd$lifehistory<-gsub(2,1,hd$lifehistory)
lifehd<-data.frame(taxa=rownames(hd),lifehistory=hd$lifehistory)
lifehd

name.check(tree, hd)
name.check(tree, lhist)
is.ultrametric(tree)
f <- c(0.77,0.71)

## Constant rate model
turnover <- c(1,1)
extinction.fraction <- c(1,1)
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0)
print(trans.rates.bisse)
dull.null <- hisse(phy=tree, data=lifehd, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=FALSE,trans.rate=trans.rates.bisse)
dull.null

## BiSSE model
turnover <- c(1,2)
extinction.fraction <- c(1,1)
BiSSE <- hisse(phy=tree, data=lifehd, f=f, turnover=turnover,eps=extinction.fraction, hidden.states=FALSE,trans.rate=trans.rates.bisse)
BiSSE

## HiSSE model
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4)
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)
print(trans.rate.hisse)
hisse <- hisse(phy=tree, data=lifehd, f=f, turnover=turnover,eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate.hisse)
hisse

## CID-2 model
turnover <- c(1, 1, 2, 2)
extinction.fraction <- rep(1, 4)
trans.rate.CID2 <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
print(trans.rate.CID2)
hiCID2 <- hisse(phy=tree, data=lifehd, f=f, turnover=turnover,eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate.CID2)
hiCID2

## CID-4 model
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
trans.rate.hiCID4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
trans.rate.hiCID4 
hiCID4 <- hisse(phy=tree, data=lifehd, f=f, turnover=turnover,eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate.hiCID4)
hiCID4

## Get the best model
model.list <- list(dull.null, BiSSE, hiCID2, hiCID4, hisse)
GetAICWeights(model.list, criterion="AIC")

AIC <- cbind(c(dull.null$loglik,BiSSE$loglik, hiCID2$loglik, hiCID4$loglik, hisse$loglik),
             c(dull.null$AICc,BiSSE$AICc, hiCID2$AICc, hiCID4$AICc, hisse$AICc),
             GetAICWeights(model.list, criterion="AIC"))
AIC
