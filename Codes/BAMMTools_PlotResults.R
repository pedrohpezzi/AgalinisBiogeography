#script adapted from https://github.com/ryanafolk/astragalus_niche_biogeo/blob/main/BAMM/bammtools_diversification.R

library("phytools")
library("ape")
library("BAMMtools")
library("openxlsx")
library("coda")

# Set WD
setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/08_BAMM")

# Check convergence
mcmcout <- read.csv("Results/concatenated_beast/concatenated_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
# Burn-in 20% of samples
burnstart <- floor(0.20 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation)

# check the effective sample sizes of the log-likelihood 
# and the number of shift events present in each sample, should be larger than 200
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

# The number of macroevolutionary rate regimes on our phylogenetic tree
post_probs <- table(postburn$N_shifts) / nrow(postburn)
# Names(post_probs)
post_probs

# Compute the posterior odds ratio for (say) two models 
post_probs["1"] / post_probs["2"] # How much more posterior probability is in 1 shifts than 2

## Summarize the posterior distribution of the number of shifts using summary methods
tree = read.tree("Results/concatenated_beast/concatenated_beast_nosynonyms_noout.tre")
edata <- getEventData(tree, eventdata = "Results/concatenated_beast/concatenated_event_data.txt", burnin=0.2)
shift_probs <- summary(edata)
shift_probs

############
# Plot diversification tree
best_diversification <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
msc.set <- maximumShiftCredibility(edata, maximize='product')
msc.config <- subsetEventData(edata, index = msc.set$sampleindex)
pdf("concatenated_diversification.pdf", width=5, height=5)
plot.bammdata(best_diversification, lwd = 0.5, method ='phylogram', labels=TRUE, spex = "netdiv", cex=0.5, logcolor = TRUE, breaksmethod = "jenks", legend=TRUE)
addBAMMshifts(best_diversification, method ='phylogram', cex=0.4)
dev.off()

# Plot rate through time
ratematrix <- getRateThroughTimeMatrix(edata) # Calculating ahead of time avoids repeating calculations to adjust figure; still need to recalculate for different nodes
# plotRateThroughTime(ratematrix,intervalCol="red", avgCol="red")
pdf("concatenated_rate_through_time.pdf", width=6, height=5)
plotRateThroughTime(ratematrix,intervalCol="skyblue", avgCol="skyblue3",ratetype="netdiv", ylim = c(0, 1))
dev.off()

# Bayes Factor
bfmat <- computeBayesFactors(postburn, expectedNumberOfShifts=1, burnin=0.2)
bfmat

# Evolutionary rates:
allrates <- getCladeRates(edata)
allrates

tiprates <- getTipRates(edata)
tiprates

# Save data to .Rsave file
save(edata, 
     allrates,
     best_diversification,
     tiprates,
     mcmcout,
     postburn,
     ratematrix,
     shift_probs,
     msc.set,
     bfmat,
     file="concatenated_BAMMFiles.Rsave")

# Get tip diversification rates in a excel file format
lambda_avg <- tiprates[["lambda.avg"]]
mu_avg <- tiprates[["mu.avg"]]
net_diversification <- lambda_avg - mu_avg
df <- data.frame(Species = names(lambda_avg), 
                 Lambda_Avg = lambda_avg, 
                 Mu_Avg = mu_avg, 
                 Net_Diversification = net_diversification)
write.xlsx(df, file = "tiprates_concatenated.xlsx")
