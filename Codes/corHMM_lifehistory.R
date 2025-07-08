# Adapted from Sazatornil et al. 2023
# https://zenodo.org/records/6507023

#library(devtools)
#install_github("thej022214/corHMM")

setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/10_corHMM")
library(corHMM)
library(phytools)
library(svglite)

# Read tree
tree <- read.tree("InputTrees/concatenated_beast_nosynonyms_noout.tre")
plot(tree,no.margin=TRUE,edge.width=2,cex=0.7)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree))

# Load lifehistory and clade data from Kriebel et al.
lifehistory <- read.csv("lifehistory.csv")
row.names(lifehistory) <- lifehistory$species

# Erase species not in tree
lifehistory <- lifehistory[tree$tip.label, ]

# Recode lifehistory syndromes 
# 1 code for annual
# 2 code for perennial
# 1&2 indicates ambiguity
lifehistory$lifehistory <- gsub("annual", "1", lifehistory$lifehistory)
lifehistory$lifehistory <- gsub("perennial", "2", lifehistory$lifehistory)
lifehistory$lifehistory <- as.factor(lifehistory$lifehistory)

# Final base
b0 <- lifehistory[,c("species","lifehistory")]


#####################################################
# 2. ANCESTRAL STATE RECONSTRUCTION ON the MCC tree #
#####################################################

####################
# 2.1 Build models
####################

# model names code for:
# SMM. Simple Markov model
# HMM. Hidden Markov model
# ARD. All rates different
# SYM. Symmetric rates

# 2.1.1 SMM_ARD
SMM_ARD <- corHMM(tree, b0, model ="ARD", 
                  node.states = "marginal", rate.cat = 1, 
                  root.p = "yang",  get.tip.states = TRUE)
SMM_ARD #summary

# 2.1.2 SMM_SYM
SMM_SYM <- corHMM(tree, b0, model ="SYM",
                  node.states = "marginal", rate.cat = 1, 
                  root.p = "yang",  get.tip.states = TRUE)
SMM_SYM # summary

# 2.1.3 HMM_ARD_ARD_CTARD
# rate transitions ARD ARD
# rate category transitions ARD
HMM_ARD_ARD_CTARD <- corHMM(tree, b0, model ="ARD",
                            node.states = "marginal", rate.cat = 2, 
                            root.p = "yang", get.tip.states = TRUE)
HMM_ARD_ARD_CTARD #summary

# 2.1.4 HMM_ARD_SYM_CTARD
# rate transitions ARD SYM
# rate category transitions ARD 
RateCat1 <- getStateMat4Dat(b0, model = "ARD", dual = T)$rate.mat # R1
RateCat2 <- getStateMat4Dat(b0)$rate.mat # R2
RateCat2 <- equateStateMatPars(RateCat2, c(1:4))
RateClassMat <- getRateCatMat(2) 
StateMats <- list(RateCat1, RateCat2)
FM_ard_sym_ctard <- getFullMat(StateMats, RateClassMat)
HMM_ARD_SYM_CTARD <- corHMM(phy = tree, data = b0, rate.cat = 2,
                            rate.mat = FM_ard_sym_ctard, node.states = "marginal", 
                            root.p = "yang", get.tip.states = TRUE)
HMM_ARD_SYM_CTARD # summary

# 2.1.5 HMM_SYM_SYM_CTARD
# rate transitions SYM SYM
# rate category transitions ARD 
RateCat1 <-getStateMat4Dat(b0)$rate.mat
RateCat1 <-equateStateMatPars(RateCat1, c(1,2)) 
RateCat2 <-getStateMat4Dat(b0)$rate.mat
RateCat2 <-equateStateMatPars(RateCat2, c(1,2))
StateMats <-list(RateCat1, RateCat2)
RateClassMat <-getRateCatMat(2)              
FM_sym_sym_ctard <-getFullMat(StateMats, RateClassMat)
HMM_SYM_SYM_CTARD <- corHMM(phy = tree, data = b0, rate.cat = 2,
                            rate.mat = FM_sym_sym_ctard, node.states = "marginal", 
                            root.p = "yang", get.tip.states = TRUE)
HMM_SYM_SYM_CTARD # summary

# 2.1.6 HMM_ARD_ARD_CTSYM
# rate transitions ARD ARD
# rate category transitions SYM
RateCat1 <-getStateMat4Dat(b0)$rate.mat
RateCat2 <-getStateMat4Dat(b0)$rate.mat
StateMats <-list(RateCat1, RateCat2)
RateClassMat <-getRateCatMat(2) 
RateClassMat <- equateStateMatPars(RateClassMat, c(1,2))
FM_ard_ard_ctsym <-getFullMat(StateMats, RateClassMat)
HMM_ARD_ARD_CTSYM <- corHMM(tree, b0, rate.mat = FM_ard_ard_ctsym,
                            node.states = "marginal", rate.cat = 2,
                            root.p = "yang", get.tip.states = TRUE)
HMM_ARD_ARD_CTSYM  # summary

# 2.1.7 HMM_ARD_SYM_CTSYM
# rate transitions ARD SYM
# rate category transitions SYM
RateCat1 <-getStateMat4Dat(b0)$rate.mat
RateCat2 <-getStateMat4Dat(b0)$rate.mat
RateCat2 <- equateStateMatPars(RateCat2, c(1,2))
StateMats <-list(RateCat1, RateCat2)
RateClassMat <-getRateCatMat(2) 
RateClassMat <- equateStateMatPars(RateClassMat, c(1,2))
FM_ard_sym_ctsym <-getFullMat(StateMats, RateClassMat)
HMM_ARD_SYM_CTSYM <- corHMM(tree, b0, rate.mat = FM_ard_sym_ctsym,
                            node.states = "marginal", rate.cat = 2, 
                            root.p = "yang", get.tip.states = TRUE)
HMM_ARD_SYM_CTSYM  # summary 

# 2.1.8 HMM_SYM_SYM_CTSYM
# rate transitions SYM SYM
# rate category transitions SYM
RateCat1 <- getStateMat4Dat(b0)$rate.mat
RateCat1 <- equateStateMatPars(RateCat1, c(1,2))  
RateCat2 <- getStateMat4Dat(b0)$rate.mat
RateCat2 <- equateStateMatPars(RateCat2, c(1,2))  
StateMats <- list(RateCat1, RateCat2)
RateClassMat <- getRateCatMat(2) 
RateClassMat <- equateStateMatPars(RateClassMat, c(1,2))
FM_sym_sym_ctsym <- getFullMat(StateMats, RateClassMat)
HMM_SYM_SYM_CTSYM <- corHMM(tree, b0, rate.mat = FM_sym_sym_ctsym,
                            node.states = "marginal", rate.cat = 2,
                            root.p = "yang", get.tip.states = TRUE)
HMM_SYM_SYM_CTSYM  # summary 

#----------------------------------------------------
save(FM_ard_sym_ctard, FM_sym_sym_ctard, 
    FM_ard_ard_ctsym, FM_ard_sym_ctsym, 
    FM_sym_sym_ctsym, file = "beast_concatenated_start_matrices_lifehistory.RData")
#----------------------------------------------------

#####################################
# 2.2 Model comparison (MCC tree)
#####################################

obj<- list(SMM_ARD, SMM_SYM, 
           HMM_ARD_ARD_CTARD, HMM_ARD_SYM_CTARD, 
           HMM_SYM_SYM_CTARD, HMM_ARD_ARD_CTSYM,
           HMM_ARD_SYM_CTSYM, HMM_SYM_SYM_CTSYM)

names(obj) <- c("SMM_ARD", "SMM_SYM", 
                "HMM_ARD_ARD_CTARD", "HMM_ARD_SYM_CTARD", 
                "HMM_SYM_SYM_CTARD", "HMM_ARD_ARD_CTSYM",
                "HMM_ARD_SYM_CTSYM", "HMM_SYM_SYM_CTSYM")

# Obtain AIC weights for model averaging
# This part of the routine was modified from Boyko & Beaulieu 2021
# https://doi.org/10.1111/2041-210X.13534 

AICcs <- unlist(lapply(obj, function(x) x$AICc))
AICwt <- exp(-0.5 * AICcs - min(AICcs))/sum(exp(-0.5 * AICcs - min(AICcs)))
res <- matrix(0, dim(obj[[1]]$states)[1], dim(obj[[1]]$states)[2])
for(i in 1:length(obj)){
  States <- colnames(obj[[i]]$solution)
  if(length(grep("R2", States)) == 0){
    ASR_i <- obj[[i]]$states[,grep("R1", States)]
  }else{
    ASR_i <- obj[[i]]$states[,grep("R1", States)] + obj[[i]]$states[,grep("R2", States)]
  }
  res <- res + (ASR_i * AICwt[i])
}
colnames(res) <- c("annual", "perennial")

# Define the nodes of interest
nodes <- seq_len(nrow(obj[[1]]$states))

# Function to extract annual probability for a given node
extr_annual <- function(x, node){
  annual <- ifelse(length(x$states[node, ]) == 2, 
                x$states[node, 1],
                sum(x$states[node, c(1,3)]))
  annual
}

# Function to extract perennial probability for a given node
extr_perennial <- function(x, node){
  perennial <- ifelse(length(x$states[node, ]) == 2, 
                 x$states[node, 2],
                 sum(x$states[node, c(2,4)]))
  perennial
}

# Initialize empty data frame to store results
table1 <- data.frame(model = names(obj), 
                     AIC_c = round(AICcs, 5),
                     weights = round(AICwt, 5))

# Loop through each node and calculate probabilities
for (node in nodes) {
  # Calculate probabilities for each model
  annual_probs <- sapply(obj, extr_annual, node = node)
  perennial_probs <- sapply(obj, extr_perennial, node = node)
  
  # Add columns for annual and perennial probabilities for each node
  table1[paste0("annual_node_", node)] <- round(annual_probs, 5)
  table1[paste0("perennial_node_", node)] <- round(perennial_probs, 5)
}

# Calculate model-averaged probabilities
averaged_probs <- sapply(nodes, function(node) {
  annual_col <- paste0("annual_node_", node)
  perennial_col <- paste0("perennial_node_", node)
  average_annual <- sum(table1[[annual_col]] * table1$weights)
  average_perennial <- sum(table1[[perennial_col]] * table1$weights)
  c(average_annual = round(average_annual, 5),
    average_perennial = round(average_perennial, 5))
})

# Display the averaged probabilities with node headers
colnames(averaged_probs) <- paste0("Node_", nodes)
averaged_probs

# Transpose the averaged_probs matrix for easier manipulation
transposed_probs <- t(averaged_probs)

# Convert to a data frame
states_df <- as.data.frame(transposed_probs)

# Rename columns to match the desired format
colnames(states_df) <- c("(annual)", "(perennial)")

# Optionally, you can remove the row names if you want clean indexing
rownames(states_df) <- NULL

# Display the final object
states_df

# Save SVG
svglite(filename = "beast_concatenated_lifehistory.svg", width = 8, height = 8)
plotRECON(tree, likelihoods = states_df, cex = 0.8,  pie.cex = 0.4,
          piecolors = c("#0072B1", "#B4E5A2"), show.tip.label = TRUE)
title("BEAST MCC Concatenated Tree")
colors <- as.factor(b0$lifehistory)
colors <- factor(colors, labels = c("#0072B1", "#B4E5A2"))
names(colors) <- b0$species
colors <- as.character(colors)
tiplabels(pch = 19, cex = 1, col = colors)
dev.off()

# Save PDF
pdf(file = "beast_concatenated_lifehistory.pdf", width = 8, height = 8)
plotRECON(tree, likelihoods = states_df, 
          cex = 0.8, pie.cex = 0.4,
          piecolors = c("#0072B1", "#B4E5A2"), 
          show.tip.label = TRUE)
title("BEAST MCC Concatenated Tree")
colors <- as.factor(b0$lifehistory)
colors <- factor(colors, labels = c("#0072B1", "#B4E5A2"))
names(colors) <- b0$species
colors <- as.character(colors)
tiplabels(pch = 19, cex = 1, col = colors)
dev.off()

#write tables
write.table(table1, file = "beast_concatenated_lifehistory_fullresults_AIC_probs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(averaged_probs, file = "beast_concatenated_lifehistory_node32_averaged_probs_table.txt", sep = "\t", quote = FALSE, row.names = TRUE)
