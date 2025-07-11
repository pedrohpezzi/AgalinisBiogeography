##############################
setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/12_MuHiSSE")

library(hisse)
library(dplyr)

tree<-read.tree(file='concatenated_beast_nosynonyms_noout.tre')
plot(tree)
tips_to_drop <- c("Agalinis_stenantha", "Agalinis_bangii", "Agalinis_reflexidens") 
tree <- drop.tip(tree, tips_to_drop)
plot(tree,no.margin=TRUE,edge.width=2,cex=0.7)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree))

poll_lifehist<- read.csv( "pollinator_lifehistory_pruned.csv")
tree<- drop.tip(tree, tree$tip.label[!tree$tip.label%in%poll_lifehist$species]) #######
row.names(poll_lifehist) <- poll_lifehist$species

poll_lifehist <- poll_lifehist %>%
  mutate(
    pollinator = ifelse(pollinator == "bee", 0, 
                        ifelse(pollinator == "bird", 1, NA)),
    lifehistory = ifelse(lifehistory == "annual", 0, 
                         ifelse(lifehistory == "perennial", 1, NA))
  )

##00 bee annual
##01 bee perennial
##10 hummingbird annual 100%
##11 hummingbird perennial

f = c(0.77,0.56,1.0,0.58)

# 00 bee_annual；01 bee_perennial；10 bird_annual；11 bird_perennial

############################################
turnover <- c(1,1,1,1) 
extinction.fraction <- rep(1, 4) 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
constant <- MuHiSSE(phy=tree, data=poll_lifehist, f=f, turnover=turnover,
                      eps=extinction.fraction, hidden.states=FALSE,
                      trans.rate=trans.rate)
constant 

###########
turnover <- c(1,2,1,2)
extinction.fraction <- rep(1, 4) 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)

trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
nohidden1 <- MuHiSSE(phy=tree, data=poll_lifehist, f=f, turnover=turnover,
                      eps=extinction.fraction, hidden.states=FALSE,
                      trans.rate=trans.rate)
nohidden1

###########
turnover <- c(1,1,2,2)
extinction.fraction <- rep(1, 4) 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)

trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
nohidden2 <- MuHiSSE(phy=tree, data=poll_lifehist, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=FALSE,
                     trans.rate=trans.rate)
nohidden2

###########
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4) 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)

trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
nohidden3 <- MuHiSSE(phy=tree, data=poll_lifehist, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=FALSE,
                     trans.rate=trans.rate)
nohidden3

#####################################################################
#00A, 01A, 10A, 11A, 00B, 01B, 10B, 11B
turnover <- c(1,1,2,2,3,3,4,4)
extinction.fraction <- rep(1, 8) 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1, make.null=TRUE)
muhisse1<- MuHiSSE(phy=tree, data=poll_lifehist, f=f, turnover=turnover,
                    eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate)
muhisse1

#######################
turnover <- c(1,2,1,2,3,4,3,4)
extinction.fraction <- rep(1, 8) 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1, make.null=TRUE)
muhisse2 <- MuHiSSE(phy=tree, data=poll_lifehist, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate)
muhisse2

######################
turnover <- c(1,2,3,4,5,6,7,8)
extinction.fraction <- rep(1, 8) 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1)

muhisse3 <- MuHiSSE(phy=tree, data=poll_lifehist, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate)
muhisse3

########################################################################
turnover <- c(1,1,1,1,2,2,2,2)
extinction.fraction <- rep(1, 8) 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1, make.null=TRUE)
null8 <- MuHiSSE(phy=tree, data=poll_lifehist, f=f, turnover=turnover,
                  eps=extinction.fraction, hidden.states=TRUE,trans.rate=trans.rate)
null8

## Select the best model

model.list <- list(constant, nohidden1, nohidden2, nohidden3, muhisse1, muhisse2, muhisse3, null8)
GetAICWeights(model.list, criterion="AIC")

AIC <- cbind(c(constant$loglik,nohidden1$loglik, nohidden2$loglik, nohidden3$loglik, muhisse1$loglik, muhisse2$loglik, muhisse3$loglik, null8$loglik),
             c(constant$AICc,nohidden1$AICc, nohidden2$AICc, nohidden3$AICc, muhisse1$AICc, muhisse2$AICc, muhisse3$AICc, null8$AICc),
             GetAICWeights(model.list, criterion="AIC"))

