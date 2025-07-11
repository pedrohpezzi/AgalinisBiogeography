library(hisse)
library(dplyr)
library(ape)

setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/12_MuHiSSE")

trees <- read.tree(file='InputTrees/100trees_nosyno_noout.nwk')
tips_to_drop <- c("Agalinis_stenantha", "Agalinis_bangii", "Agalinis_reflexidens")
trees_dropped <- lapply(trees, drop.tip, tip = tips_to_drop)
class(trees_dropped) <- "multiPhylo"
poll_lifehist <- read.csv("pollinator_lifehistory_pruned.csv")
row.names(poll_lifehist) <- poll_lifehist$species

# Format data
poll_lifehist <- poll_lifehist %>%
  mutate(
    pollinator = ifelse(pollinator == "bee", 0, ifelse(pollinator == "bird", 1, NA)),
    lifehistory = ifelse(lifehistory == "annual", 0, ifelse(lifehistory == "perennial", 1, NA))
  )

# Loop through trees
for (i in 1:length(trees_dropped)) {
  tree <- trees_dropped[[i]]
  tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% poll_lifehist$species])
  
  data <- poll_lifehist[tree$tip.label, ]
  f <- c(0.77, 0.56, 1.0, 0.58)
  
  ## Models
  trans0 <- TransMatMakerMuHiSSE(hidden.traits = 0)
  trans1_null <- TransMatMakerMuHiSSE(hidden.traits = 1, make.null = TRUE)
  trans1 <- TransMatMakerMuHiSSE(hidden.traits = 1)
  constant <- MuHiSSE(phy = tree, data = data, f = f, turnover = rep(1, 4), eps = rep(1, 4), hidden.states = FALSE, trans.rate = trans0)
  nohidden1 <- MuHiSSE(phy = tree, data = data, f = f, turnover = c(1,2,1,2), eps = rep(1, 4), hidden.states = FALSE, trans.rate = trans0)
  nohidden2 <- MuHiSSE(phy = tree, data = data, f = f, turnover = c(1,1,2,2), eps = rep(1, 4), hidden.states = FALSE, trans.rate = trans0)
  nohidden3 <- MuHiSSE(phy = tree, data = data, f = f, turnover = c(1,2,3,4), eps = rep(1, 4), hidden.states = FALSE, trans.rate = trans0)
  muhisse1 <- MuHiSSE(phy = tree, data = data, f = f, turnover = c(1,1,2,2,3,3,4,4), eps = rep(1, 8), hidden.states = TRUE, trans.rate = trans1_null)
  muhisse2 <- MuHiSSE(phy = tree, data = data, f = f, turnover = c(1,2,1,2,3,4,3,4), eps = rep(1, 8), hidden.states = TRUE, trans.rate = trans1_null)
  muhisse3 <- MuHiSSE(phy = tree, data = data, f = f, turnover = c(1,2,3,4,5,6,7,8), eps = rep(1, 8), hidden.states = TRUE, trans.rate = trans1)
  null8 <- MuHiSSE(phy = tree, data = data, f = f, turnover = c(1,1,1,1,2,2,2,2), eps = rep(1, 8), hidden.states = TRUE, trans.rate = trans1_null)

  # Save results
  save(constant, nohidden1, nohidden2, nohidden3,
       muhisse1, muhisse2, muhisse3, null8,
       file = paste0("Tree", i, "_muhisse.Rsave"))
}
