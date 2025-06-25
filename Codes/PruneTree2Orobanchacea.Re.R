library("ape")
library("phytools")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("readr")

#http://blog.phytools.org/2014/11/pruning-trees-to-one-member-per-genus.html
# Script used to prune the Lamiales tree from Fonseca et al. (2021) down to only Orobanchaceae species.  
# The file `OrobanchaceaeSpList.txt` is in the `Data/congruification` folder.

setwd("C:/Users/pedro/OneDrive/Desktop/Agalinis/Analyses/Congruification")

tree <- read.tree("Fonseca2021_SuppMaterial/1-s2.0-S1055790321002207-mmc1/Lamiales_dated.nwk.tre")
ggtree(tree)
ggtree(tree, layout='circular')

splist <- read_lines("InputData_Orobanchaceae/OrobanchaceaeSpList.txt")

pruned.tree<-drop.tip(tree,tree$tip.label[-match(splist, tree$tip.label)])
plot(pruned.tree)
is.rooted(pruned.tree)

write.tree(pruned.tree, file = "Fonseca_OrobanPauloPruned.tre")
