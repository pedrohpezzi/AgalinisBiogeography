library(hisse)
library(diversitree)
library(geiger)
library(phytools)
library(ape)

setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/HiSSE")

# Load the trees
trees <- read.tree(file='InputTrees/100_pruned_trees.newick')
tips_to_drop <- c("Agalinis_stenantha", "Agalinis_bangii", "Agalinis_reflexidens") 
trees_dropped <- lapply(trees, drop.tip, tip = tips_to_drop)
class(trees_dropped) <- "multiPhylo"

# Check if the trees are ultrametric
are.ultrametric <- sapply(trees_dropped, is.ultrametric)
print(are.ultrametric)

# Load the pollination data
po <- read.csv(file='pollinator_pruned.csv', header=TRUE)
rownames(po) <- po$species
poli <- po$pollinator
names(poli) <- po$species

# Loop through each tree
for (i in 1:length(trees_dropped)) {
  tree <- trees_dropped[[i]]
  
  # Match pollination data to tree tips
  pol <- poli[names(poli) %in% tree$tip.label]
  pol <- as.factor(pol)
  x <- droplevels(pol)
  x <- x[order(factor(names(x), levels=tree$tip.label))]

  
  # Create data frame for HiSSE
  hd<-cbind(rownames(x),x)
  colnames(hd)<-'pollination'
  hd<-data.frame(hd)
  hd$pollination<-gsub(1,0,hd$pollination)
  hd$pollination<-gsub(2,1,hd$pollination)
  polhd<-data.frame(taxa=rownames(hd),pollination=hd$pollination)

    #define sampling fraction
  f <- c(0.71,0.58)
  
  # Run HiSSE models
  dull.null <- hisse(phy=tree, data=polhd, f=f, turnover=c(1,1), eps=c(1,1), hidden.states=FALSE, trans.rate=TransMatMakerHiSSE(hidden.traits=0))
  BiSSE <- hisse(phy=tree, data=polhd, f=f, turnover=c(1,2), eps=c(1,1), hidden.states=FALSE, trans.rate=TransMatMakerHiSSE(hidden.traits=0))
  hiCID2 <- hisse(phy=tree, data=polhd, f=f, turnover=c(1,1,2,2), eps=rep(1,4), hidden.states=TRUE, trans.rate=TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE))
  hiCID4 <- hisse(phy=tree, data=polhd, f=f, turnover=c(1,1,2,2,3,3,4,4), eps=rep(1,8), hidden.states=TRUE, trans.rate=TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE))
  hissepoll <- hisse(phy=tree, data=polhd, f=f, turnover=c(1,2,3,4), eps=rep(1,4), hidden.states=TRUE, trans.rate=TransMatMakerHiSSE(hidden.traits=1))
  
  # Save all models to the same Rsave file for the current tree
  save(dull.null, BiSSE, hiCID2, hiCID4, hissepoll, file = paste0("Tree", i, "_hisse.Rsave"))
}
