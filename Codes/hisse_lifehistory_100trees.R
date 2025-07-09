library(hisse)
library(diversitree)
library(geiger)
library(phytools)

setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/11_HiSSE")


trees <- read.tree(file='InputTrees/100trees_nosyno_noout.nwk')
are.ultrametric <- sapply(trees, is.ultrametric)
print(are.ultrametric)
lh <- read.csv(file='lifehistory.csv', header=TRUE)
rownames(lh) <- lh$species
lifehistory <- lh$lifehistory
names(lifehistory) <- lh$species

# Loop through each tree
for (i in 1:length(trees)) {
  tree <- trees[[i]]
  
  # Match pollination data to tree tips
  lifehistory<-lh$lifehistory
  names(lifehistory)<-lh$species
  life<-lifehistory[names(lifehistory)%in%tree$tip.label]
  life<-as.factor(life)
  x<-droplevels(life)
  x<-x[order(factor(names(x),levels=tree$tip.label))]
  
  hd<-cbind(rownames(x),x)
  colnames(hd)<-'lifehistory'
  hd<-data.frame(hd)
  hd$lifehistory<-gsub(1,0,hd$lifehistory)
  hd$lifehistory<-gsub(2,1,hd$lifehistory)
  lifehd<-data.frame(taxa=rownames(hd),lifehistory=hd$lifehistory)
  lifehd
  
  
  #define sampling fraction
  f <- c(0.77,0.71)
  
  # Run HiSSE models
  constant <- hisse(phy=tree, data=lifehd, f=f, turnover=c(1,1), eps=c(1,1), hidden.states=FALSE, trans.rate=TransMatMakerHiSSE(hidden.traits=0))
  BiSSE <- hisse(phy=tree, data=lifehd, f=f, turnover=c(1,2), eps=c(1,1), hidden.states=FALSE, trans.rate=TransMatMakerHiSSE(hidden.traits=0))
  hiCID2 <- hisse(phy=tree, data=lifehd, f=f, turnover=c(1,1,2,2), eps=rep(1,4), hidden.states=TRUE, trans.rate=TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE))
  hiCID4 <- hisse(phy=tree, data=lifehd, f=f, turnover=c(1,1,2,2,3,3,4,4), eps=rep(1,8), hidden.states=TRUE, trans.rate=TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE))
  hisse <- hisse(phy=tree, data=lifehd, f=f, turnover=c(1,2,3,4), eps=rep(1,4), hidden.states=TRUE, trans.rate=TransMatMakerHiSSE(hidden.traits=1))

  save(constant, BiSSE, hiCID2, hiCID4, hisse, 
       file = file.path("Results100trees_lifehistory", paste0("Tree", i, "_hisse.Rsave")))
}

