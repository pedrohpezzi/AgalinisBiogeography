library(ape)
setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/100randomtrees")

# Load files form BEAST
mcmc_file <- "TreeCombined_Sigma1.0.trees"
mcmc_trees <- read.nexus(mcmc_file)
num_trees <- length(mcmc_trees)

# select trees at least 1000 positions apart
select_trees <- function(trees, num_selected, min_distance) {
  selected_indices <- c()
  while (length(selected_indices) < num_selected) {
    candidate <- sample(seq_along(trees), 1)
    if (all(abs(candidate - selected_indices) >= min_distance)) {
      selected_indices <- c(selected_indices, candidate)
    }
  }
  trees[selected_indices]
}

selected_trees <- select_trees(mcmc_trees, 100, 1000)

# write each to a newick file
for (i in seq_along(selected_trees)) {
  write.tree(selected_trees[[i]], file = paste0("tree_", i, ".nwk"))
}

write.tree(selected_trees, file = "100trees.nwk")