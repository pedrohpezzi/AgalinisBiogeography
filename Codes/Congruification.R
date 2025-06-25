#install.packages("geiger")
#install.packages("ape")
library(geiger)
library(ape)

setwd("C:/Users/pedro/OneDrive - University of Arkansas/Agalinis/Analyses/02_Congruification")

Sys.getenv("PATH")
Sys.setenv(PATH=paste("C:/Users/pedro/OneDrive/Desktop/Softwares/PATHd8/PATHd8", Sys.getenv("PATH"), sep=":"))
Sys.getenv("PATH")

#Run
reference <- read.tree("InputData_Orobanchaceae/Fonseca_OrobanPauloPruned.tre")
plot(reference)
is.rooted(reference)

target <- read.tree("InputData_Orobanchaceae/AgalinisRooted_MrBayes.tre")
plot(target)
is.ultrametric(ultra_target)
is.rooted(target)

tax <- read.table("InputData_Orobanchaceae/TaxaEquivalency.txt")

congru_tree <- congruify.phylo(reference, target, tax, tol = 0, scale="PATHd8", ncores=4)

print(congru_tree)

cal=congru_tree$calibrations

write.treePL(target, cal, opts=list(smooth=0.1, nthreads=2, opt=1, optad=1, thorough=TRUE))