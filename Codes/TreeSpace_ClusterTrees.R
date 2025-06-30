#library(devtools)
#install_github("thibautjombart/treespace")
#install.packages("igraph", configure.args = "--disable-glpk")
library("treespace")
library("adegenet")
library("adegraphics")
library("rgl")

setwd("/storage/ppezzi/BEAST_Agalinis/clustertrees")

trees <- read.nexus("BEAST_concatenated_runs_combined.trees")

# choose random 75k trees - 25% of trees
subsamples <- trees[sample(1:length(trees),75000)]

# save the subsamples file 
write.nexus(subsamples, file ="75_subsamples.tre")

# use treespace to find and project the distances:
Dscape <- treespace(trees, nf=3)

# simple plot 
plotGrovesD3(Dscape$pco)

# find clusters:
BEASTGroves <- findGroves(Dscape, nclust=10)

# plot results
plotGrovesD3(BEASTGroves, legend_width=50, col_lab="Cluster")
# plot axes 2 and 3. This helps to show why, for example, clusters X and Y have been identified as separate, despite them appearing to overlap when viewing axes 1 and 2.
plotGrovesD3(BEASTGroves, xax=2, yax=3, legend_width=50, col_lab="Cluster")

# find median trees for the 10 clusters identified earlier:
res <- medTree(trees, BEASTGroves$groups)
names(res)

# get the first median of each
med.trees <- lapply(res, function(e) ladderize(e$trees[[1]]))

# plot trees
par(mfrow=c(2,2))
for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i),cex=0.5)

# compare median trees from clusters x and y - change numbers to compare different trees
plotTreeDiff(med.trees[[1]],med.trees[[2]], use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)

# loop through each tree and save it to a separate file
for (i in 1:length(med.trees)) {
  write.nexus(med.trees[[i]], file = paste0("median_tree_cluster_", i, ".nex"))
}
# save all trees in one file
write.nexus(med.trees, file = "median_trees_all_clusters.nex")
