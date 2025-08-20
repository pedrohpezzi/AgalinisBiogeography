# http://phylo.wikidot.com/biogeobears

###########################################################
# Run Biogeobears on multiple trees and store the output 
###########################################################

###########################################################
#Load libraries and additional source code
###########################################################
library(optimx)         
library(FD)       
library(snow)     
library(parallel)
library(BioGeoBEARS)
library(stringr)

source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")

calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

############################################################
# Universal setup, input file path
############################################################
setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/07_BioGeoBEARS_multipletrees")
trfn = "concatenated_tree_nosynonyms.nwk"
geogfn = "GeoData.txt"
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
max_range_size = 3

newick_fns = scan("treefiles.txt", what="character", sep=NULL)
geog_fns = scan("geodatafiles.txt", what="character", sep=NULL)
resfns = scan("outputnames.txt", what="character", sep=NULL)


############################################################
# Run DEC+j on multiple trees
############################################################
### Set BioGeoBEARS object - DEC+j
# Run DEC+J
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE         # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();
# if "GenSA", use Generalized Simulated Annealing, which seems better on high-dimensional
# problems (5+ parameters), but seems to sometimes fail to optimize on simple problems
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = res$outputs@params_table["d","est"]
estart = res$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


# Run the analysis, on the master tree
# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "BGB_originaltree.Rdata"
if (runslow)
    {
    master_tree_res = bears_optim_run(BioGeoBEARS_run_object)
    master_tree_res

    save(master_tree_res, file=resfn)
    resDEC = master_tree_res
    } else {
    # Loads to "master_tree_res"
    load(resfn)
    resDEC = master_tree_res
    }


# multiple tree settings
#' @BioGeoBEARS_run_object Set up the inference model you want to run on 
#' each tree, like you would for a normal single-tree run.
#' @newick_fns A list of the Newick files (e.g., you should extract some trees
#' from a BEAST NEXUS MCMC output)
#' @model_name The name you would like added to output filenames
#' @geog_fns A list of corresponding geography files (by default, these are just .geog instead of .newick)
#' @resfns A list of results filenames for .Rdata files, either for saving BioGeoBEARS analyses on each tree, or 
#' for loading previously-saved runs (by default, these are just _model_name.Rdata instead of .newick)
#' @run_or_collect If you want to run BioGeoBEARS on each tree (slow), use "run". If you just want to 
#' collect the results over all the trees, use "collect".  For both, pick "both".
#' @start_treenum Default 1. Change if you want to skip some trees
#' @end_treenum Default is length(newick_fns). Change if you want to run a subset of trees.
#' @runslow If FALSE, old, saved .Rdata results are loaded via load(). Default TRUE generates new .Rdata files and 
#' saves them via save()

runslow = TRUE
if (runslow)
	{
	res = run_bears_optim_on_multiple_trees(BioGeoBEARS_run_object, newick_fns, model_name="AgalinisDECj", geog_fns, resfns, run_or_collect="run", start_treenum=1, end_treenum=length(newick_fns), runslow=TRUE, plot_params=FALSE)
	
	# Save res to "res_over_multiple_trees" 
	res_over_multiple_trees = res
	save(res_over_multiple_trees, file="res_over_multiple_trees.Rdata")
	} else {
	# Loads to "res_over_multiple_trees"
	load(file="res_over_multiple_trees.Rdata")
	} # END if (runslow)



# Now, collect the saved results from each individual run
res2 = run_bears_optim_on_multiple_trees(BioGeoBEARS_run_object, newick_fns, model_name="AgalinisDECj", geog_fns, resfns, run_or_collect="collect", start_treenum=1, end_treenum=length(newick_fns), runslow=TRUE, plot_params=FALSE)
names(res2)

# This has the matrix of state probabilities at nodes
# for ALL of the trees you specified, in one huge matrix. 
dim(res2$state_probs_at_nodes_across_all_trees)
# The rightmost column is text, and contains the alphabetical, 
# comma-delimited list of tipnames (OTU names) descending from 
# that node or corner.
#
# These tipnames are used to identify "common" nodes between
# individual trees and the master tree. (See debates about
# the meaning of "common nodes", below.)
#
# The function summarize_stateprobs_on_master_tree(), below, is 
# what will go through this matrix, extract the probabilities
# for common nodes, and produce a matrix appropriate for
# plotting with standard functions.

# This matrix has probabilities for corners (branch bottoms), 
# instead of nodes.
dim(res2$state_probs_at_corners_across_all_trees)

# This has the averaged parameter values
res2$optim_results_mean


#################################################################
# To average the probabilities across all trees
#################################################################

# To plot with the standard BioGeoBEARS plot function you will need to 
# - run a standard analysis on the master tree, and 
# - then copy these average state probabilities into that res object 
#   (into the state probabilities at nodes, and state probabilities at branch bottoms)
# - then copy res2$optim_results_mean over the standard optim results


# Let's assume this tree is the master tree
best_trfn = "concatenated_tree_nosynonyms.nwk"
tr = read.tree(best_trfn)
master_OTUs = tr$tip.label


stateprobs_list = summarize_stateprobs_on_master_tree(master_tree_fn=best_trfn, state_probs_at_nodes_across_all_trees=res2$state_probs_at_nodes_across_all_trees, state_probs_at_corners_across_all_trees=res2$state_probs_at_corners_across_all_trees, plotflag=FALSE)
stateprobs_list

#####################################################
# Get the single-tree master tree analysis
#####################################################
# Load the tipranges
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Load master-tree analysis to "master_tree_res"
load(resfn)
resDECj_averaged = master_tree_res



#####################################################
# Input the averaged probabilities
#####################################################

# At nodes
resDECj_averaged$ML_marginal_prob_each_state_at_branch_top_AT_node = stateprobs_list$stateprobs_nodes_master_tr

# At corners
resDECj_averaged$ML_marginal_prob_each_state_at_branch_bottom_below_node = stateprobs_list$stateprobs_corners_master_tr

# Average parameter inferences
resDECj_averaged$optim_result = res2$optim_results_mean


#####################################################
# PLEASE NOTE: As you can see, 
# SOME NODES HAVE NO STATES OR PIE CHARTS. This is 
# because those nodes were not found in any of 
# the 3 trees sampled in example, above.
#
# In the probabilities matrix, these will be rows 
# that have "NA" instead of numbers.
#
# (Also see discussion of the ambiguity of what 
#  "same node" even means, in this context.)
#####################################################




#######################################################
# PDF plots
#######################################################
pdffn = "AverageOf100trees_on_bestBEASTtree.pdf"
pdf(pdffn, width=6, height=18)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DECj, averaged over 100 trees, plotted on best BEAST tree"

# Setup
results_object = resDECj_averaged
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it




####################################
# Probabilities of states/ranges at each node
####################################
# To get the probabilities of each state/range at each node:
# What you want, if "res" is your results object, is:
resDECj_averaged$ML_marginal_prob_each_state_at_branch_top_AT_node

# In this table:
# - columns are states/ranges
# - rows are nodes, in APE order (tips, then root, then internal)

#  You can see the node numbers in the same APE order with:
trtable = prt(tr, printflag=FALSE)
head(trtable)
tail(trtable)

# You can plot APE node labels with:
plot(tr, cex = 0.7)
axisPhylo()
nodelabels(cex = 0.5, frame = "n")
tiplabels(1:length(tr$tip.label), cex = 0.5, frame = "n")

# Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
max_range_size =3
areas = getareas_from_tipranges_object(tipranges)

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based))
{    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

# Look at the ranges list
ranges_list

# Make the node numbers the row names
# Make the range_list the column names
range_probabilities = as.data.frame(resDECj_averaged$ML_marginal_prob_each_state_at_branch_top_AT_node)
row.names(range_probabilities) = trtable$node
names(range_probabilities) = ranges_list

# Look at the table (first six rows)
head(range_probabilities)

# Write the table to a tab-delimited text file (for Excel etc.)
write.table(range_probabilities, file="range_probabilities_100averaged.txt", quote=FALSE, sep="\t")

# Look at the file
moref("range_probabilities_100averaged.txt")








###############################################################
# Nick's advice/rant on the pseudo-Bayesian approach
###############################################################
# 
# After you've done all this, you'll have a plot of probabilities of states, 
# averaged over many trees. Congratulations.
# 
# However, I bet it will probably look pretty much like a plot of one run on the single ML or MCC tree.
# This is one of several reasons why I am skeptical of the popularity of "Bayes DIVA", "Bayes Lagrange", etc. 
# approaches. 
#
# SEVERAL POTENTIAL PROBLEMS WITH THE (PSEUDO) "BAYES" APPROACHES 
# TO HISTORICAL BIOGEOGRAPHY (Bayes DIVA, Bayes-Lagrange) 
#
# 1. First, they aren't really fully Bayesian -- more like pseudo-Bayesian.  
#
# 2. Second, the "Bayes" approach is sometimes done to mollify reviewers who heard somewhere 
# it was a good idea, probably because they think it helps "account for uncertainty." 
# However, The REALLY BIG uncertainty, is in things like:
#
# 2a. Model misspecification. For example, Bayes-DIVA is still giving DIVA results, and 
#     DIVA may be a poorly-fitting model on many datasets.  Running it over 2000 trees 
#     or whatever doesn't fix any of that, despite the elaborate and impressive-sounding
#     computational effort.
#
# 2b. Missing species (due to incomplete sampling, or extinction) will be missing in any
#     analysis, and this will not be fixed by running over many trees.
# 
# 2c. The "single most probable state" plots hide the uncertainty, because they just plot
#     the single-most probable estimated ancestral state at each node anyway. You need 
#     to do stochastic mapping to real feel, and characterize, the uncertainty.
#
# 3.  A third problem is that "averaging" "node probabilities" across trees is a weird 
#     operation. What is the definition of "same node" between two different trees?
#
# Consider these two clades, from two sampled trees:
# tree 1: ((A,B),C)
# tree 2: ((A,B),D)
#
# Is the common ancestor of (A,B) the "same" between the two trees?  The analysis above
# assumes "yes".
#
# (And, I think, all Bayes-DIVA/Bayes-Lagrange approaches assume this, except perhaps 
# the one I implemented in Wood, Matzke et al. (2013), Sys. Bio., which had more options.)
#
# Also, consider the "corner" at the bottom of the branch of the (A,B) node. Is *it* the same between
# the two trees?  The BioGeoBEARS program above assumes "yes".
#
# How about:
# tree 1: (((A,B),C),D)
# tree 2: (((A,C),B),D)
#
# Is the node that is common ancestor to ((A,B),C) "the same" as the node that is common ancestor to ((A,C),B)?
# Sort of, but not really. But, BioGeoBEARS and presumably the other programs above assume "yes".
#
# You can imagine getting "better matches" by requiring more constraints (same daughter tips, and same sister),
# but this will reduce matches between trees, and still doesn't really solve the problem, because unless
# trees agree perfectly in topology, there will always be mismatches of this sort.  And because biogeographical
# models have sudden cladogenetic transitions in addition to the continuous-time anagenetic ones, topological
# differences may have more significance than they would under typical, purely continuous-time models.
# 
#
# WHAT TO DO INSTEAD
#
# After thinking about it, my current conclusion is that the pseudo-Bayesian approaches
# need some serious critical thinking before being employed.
#
# Rather than running hundreds or thousands of trees and presenting dubious "average",
#
# * It would be better to spend time *thinking* about and discussing all of the sources
# of uncertainty, above, rather than spending days running a model on a bunch of (typically) very 
# similar trees that (typically) all say about the same thing.
#
# * It would be better to just run the analysis on a few different trees, and put those
#   plots in the paper.
#
# * If you are really concerned about some specific events, run a number of different trees, and then
#   do stochastic mapping on each one. The summary counts of events would be one way of summarizing
#   the uncertainty in event histories.  You could also identify some "common" nodes by some
#   explicitly-stated criterion, and count up the ranges inhabited and the counts of events at those nodes.
#
# The best argument for doing the pseudo-Bayesian approach is if you have a particular 
# node you are interested in, and the geography estimated at that node is strongly 
# affected by nearby uncertainty in the tree.  In that case,
# running a model over a bunch of trees can be useful.  
#
# I may be overly skeptical. But I have yet to see a case where the pseudo-Bayesian approach
# gave any real additional insight. It's only likely if your tree is very uncertain, but
# typically people get enough data to avoid that. Very often, I guess, people just want to
# be able to say, "we ran it on 100 trees and got the same result".  But really, if your
# tree is decently estimated, you actually knew that going in, from first principles.
#
# The pseudo-Bayesian approach might be useful if your dating is highly uncertain (as is common), 
# and you are doing a time-stratified analysis where uncertain nodes are near the time-strata.
# This would be the case I would focus on.
