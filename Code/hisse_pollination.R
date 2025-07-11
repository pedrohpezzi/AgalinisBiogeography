#script adapted from Valderrama et al. (2022) - https://doi.org/10.3389/fpls.2022.874322
setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/11_HiSSE")

#install.packages('hisse')
library(hisse)
library(diversitree)
library(geiger)
library(phytools)

tree<-read.tree(file='InputTrees/concatenated_beast_nosynonyms_noout.tre')
plot(tree)
tips_to_drop <- c("Agalinis_stenantha", "Agalinis_bangii", "Agalinis_reflexidens") 
tree <- drop.tip(tree, tips_to_drop)
plot(tree,no.margin=TRUE,edge.width=2,cex=0.7)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree))
po<-read.csv(file='pollinator_pruned.csv',h=T)
rownames(po)<-po$species

poli<-po$pollinator
names(poli)<-po$species
pol<-poli[names(poli)%in%tree$tip.label]
pol<-as.factor(pol)
x<-droplevels(pol)
x<-x[order(factor(names(x),levels=tree$tip.label))]
x

hd<-cbind(rownames(x),x)
colnames(hd)<-'pollination'
hd<-data.frame(hd)
hd$pollination<-gsub(1,0,hd$pollination)
hd$pollination<-gsub(2,1,hd$pollination)
polhd<-data.frame(taxa=rownames(hd),pollination=hd$pollination)
polhd
name.check(tree, hd)
name.check(tree, po)
is.ultrametric(tree)
f <- c(0.71,0.58)

## Constant rate model
turnover <- c(1,1)
extinction.fraction <- c(1,1)
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0)
print(trans.rates.bisse)
dull.null <- hisse(phy=tree, data=polhd, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=FALSE,trans.rate=trans.rates.bisse)
dull.null

## BiSSE model
turnover <- c(1,2)
extinction.fraction <- c(1,1)
BiSSE <- hisse(phy=tree, data=polhd, f=f, turnover=turnover,eps=extinction.fraction, hidden.states=FALSE,trans.rate=trans.rates.bisse)
BiSSE

## HiSSE model
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4)
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)
print(trans.rate.hisse)
hissepoll <- hisse(phy=tree, data=polhd, f=f, turnover=turnover,eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate.hisse)
hissepoll

## CID-2 model
turnover <- c(1, 1, 2, 2)
extinction.fraction <- rep(1, 4)
trans.rate.CID2 <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
print(trans.rate.CID2)
hiCID2 <- hisse(phy=tree, data=polhd, f=f, turnover=turnover,eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate.CID2)
hiCID2

## CID-4 model
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
trans.rate.hiCID4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
trans.rate.hiCID4 
hiCID4 <- hisse(phy=tree, data=polhd, f=f, turnover=turnover,eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate.hiCID4)
hiCID4

## Get the best model
model.list <- list(dull.null, BiSSE, hiCID2, hiCID4, hissepoll)
GetAICWeights(model.list, criterion="AIC")

AIC <- cbind(c(dull.null$loglik,BiSSE$loglik, hiCID2$loglik, hiCID4$loglik, hissepoll$loglik),
             c(dull.null$AICc,BiSSE$AICc, hiCID2$AICc, hiCID4$AICc, hissepoll$AICc),
             GetAICWeights(model.list, criterion="AIC"))
AIC
