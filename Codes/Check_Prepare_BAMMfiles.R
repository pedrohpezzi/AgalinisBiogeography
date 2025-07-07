library(ape)
library(phytools)

setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/08_BAMM")
tree <- read.tree("TREE_INPUT.tre")
is.ultrametric(tree)
plot(tree)
tree <- force.ultrametric(tree)
is.ultrametric(tree)
is.binary(tree)
# Now to check min branch length:
min(tree$edge.length)
plot(tree)

#drop outgroup tips
tips <- c("Aureolaria_pedicularia", "Brachystigma_wrightii", "Dasistoma_macrophylla")
tree2 <- drop.tip(tree, tips)
plot(tree2)

write.tree(tree2, file = "TREE_OUTPUT.tre")
######################################################################################################

tree <- read.tree("TREE_OUTPUT.tre")

#set BAMM priors
#install.packages("BAMMtools")
library(BAMMtools)

setBAMMpriors(
  tree,
  total.taxa = 72,
  traits = NULL,
  outfile = "priors.txt",
  Nmax = 1000,
  suppressWarning = FALSE
)
