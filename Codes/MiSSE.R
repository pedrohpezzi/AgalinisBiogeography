# Adapted from Wilson et al. (in prep) from the Beaulieu lab - University of Arkansas

#install.packages("hisse")
library(hisse)
library(parallel)
detectCores()

setwd("/scrfs/storage/ppezzi/misse")

DoOneRun <- function(iteration, phy, f){
  possible.combos <- generateMiSSEGreedyCombinations(max.param=52, turnover.tries=sequence(5), eps.tries=sequence(5), vary.both=TRUE, fixed.eps.tries=NA)
  possible.combos$fixed.eps <- as.numeric(possible.combos$fixed.eps)
  model.set.misse <- MiSSEGreedy(phy=phy, possible.combos=possible.combos, stop.deltaAICc=10, chunk.size=4, f=f, root.type = "herr_als", n.cores=32)
  model.recons.misse <- as.list(1:length(model.set.misse))
  model.set.misse.pruned <- PruneRedundantModels(model.set.misse)
  for (model_index in 1:length(model.set.misse.pruned)) {
    nturnover <- length(unique(model.set.misse.pruned[[model_index]]$turnover))
    neps <- length(unique(model.set.misse.pruned[[model_index]]$eps))
    hisse_recon <- MarginReconMiSSE(phy=model.set.misse.pruned[[model_index]]$phy, f=f, hidden.states=nturnover, pars=model.set.misse.pruned[[model_index]]$solution, AIC=model.set.misse.pruned[[model_index]]$AIC, root.type="herr_als", n.cores = 32)
    model.recons.misse[[model_index]] <- hisse_recon
  }
  tip.rates.misse <- GetModelAveRates(model.recons.misse[[model_index]], type = c("tips"))
  tip.rates.misse <- tip.rates.misse[order(tip.rates.misse$taxon),]
  
  save(model.set.misse, model.set.misse.pruned, model.recons.misse, tip.rates.misse, file=paste("Agalinis_tree", iteration, ".Rsave", sep=""))
}


DoAllTrees <- function(phy, f){
  ntrees <- length(phy)
  for(tree.index in 1:ntrees){
    cat("Running tree number:", tree.index, "\n")
    tmp <- DoOneRun(tree.index, phy[[tree.index]], f=f)
  }
}

#100 trees
trees <- read.tree("trees.trees")
f <- 0.72
DoAllTrees(phy=trees, f=f)
