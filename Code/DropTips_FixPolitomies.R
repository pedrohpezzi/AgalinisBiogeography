setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/BioGeoBEARS")

library(ape)
library(RRphylo)

# load my data
trfn = "tree.nwk"
tr = read.tree(trfn)
tr
plot(tr, cex=0.5)
axisPhylo() # plots timescale
tr$edge.length
sum(tr$edge.length == 0)
TF = tr$edge.length == 0
TF
nums = (1:length(tr$edge.length))[TF]
nums
tr$edge.length[nums]

# drop tips - species that are synonyms
tips <- c("Agalinis_acuta", "Agalinis_paupercula")
pruned_tree <- drop.tip(tr, tips)
plot(pruned_tree)
write.tree(pruned_tree, file = "tree_nosynonyms.nwk")

# fix politomy for r8s and congruification trees
tr2 = fix.poly(pruned_tree, type="resolve")
is.binary(tr2)
tr2$edge.length
sum(tr2$edge.length == 0)
TF = tr2$edge.length == 0
TF
nums = (1:length(tr2$edge.length))[TF]
nums
tr2$edge.length[nums]
tr2$edge.length<-ifelse(test=tr2$edge.length==0, yes=0.00001, no=tr2$edge.length)
write.tree(tr2, file = "tree_nosynonyms_bifurcating.nwk")
plot(tr2)
