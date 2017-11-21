## Run genetic_msa by reading in BAYES input files (.bse and .mix)
rm(list = ls())
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/Spring Troll 2010-2017/")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
# dir.create("genetic_msa")

# Get objects/ define arguments
LocusControl <- dget(file = "Objects/LocusControl.txt")
loci <- dget(file = "Objects/GAPSLoci_reordered.txt")
basename <- "GAPS357pops13loci"
mixname <- "D114Troll_2010"
groupnames <- dget(file = "Objects/GroupNames26.txt")
groups <- dget(file = "Objects/GroupVec26RG_357.txt")
ngroups <- ifelse(length(groupnames) == max(groups), max(groups), stop("groupnames and groups are not equal!!!"))
group_prior <- rep(1 / max(groups), max(groups))
nchains = 5
group_inits <- MultiChainInits.GCL(npops = ngroups, nchains = nchains, prop = 0.9)
nits = 4e4
nburn = 2e4
thin = 10
q_out = FALSE
level = 0.1
nalleles <- LocusControl$nalleles[loci]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create directory for mixture
if(!dir.exists("genetic_msa")) {stop("Need to create 'genetic_msa' directory!!!")}
dir.create(paste0("genetic_msa/", mixname))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in .bse
bayes.bse <- as.matrix(read.table(file = paste0("BAYES/Baseline/", basename, ".bse")))
npops <- max(bayes.bse[, 1])
nloci <- max(bayes.bse[, 2])

# Check if nloci from the BAYES .bse file is the same as 'loci'
if(nloci != length(loci)) {stop("Length of 'loci' is not equal to the number of loci in the BAYES .bse file!!!")}

# Create base.txt
base <- t(sapply(seq(npops), function(i) {
  Reduce('c', 
         sapply(seq(nloci), function(l) {
           row.il <- which(bayes.bse[, 1] == i & bayes.bse[, 2] == l)
           bayes.bse[row.il, seq(nalleles[l]) + 3]
         } )
  )
} ))

# Write base.txt
write.table(x = base, file = paste0("genetic_msa/", mixname, "/base.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in .mix
test.mix <- as.matrix(read.table(paste0("BAYES/Mixture/", mixname, ".mix")) )
nind <- nrow(test.mix)
if(nloci != ncol(test.mix)) {stop("Length of 'loci' is not equal to the number of loci in the BAYES .mix file!!!")}

# Create mix.txt
bayes.mix <- as.matrix(read.table("BAYES/Mixture/D101102Troll_2010.mix", colClasses = rep("character", nloci))) 
mix <- t(sapply(
  apply(bayes.mix, 1, function(ind) {paste(ind, collapse = '')
  } ),
  function(ind.c) {as.numeric(unlist(strsplit(x = ind.c, split = "")))
  } ) 
)

# Write mix.txt
write.table(x = mix, file = paste0("genetic_msa/", mixname, "/mix.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Start genetic_msa
source("C:/Users/krshedd/Documents/R/GCL-R-Scripts/genetic_msa.R")

orig_dir <- getwd()
setwd(paste0("genetic_msa/", mixname))

genetic_msa_Estimates <- genetic_msa(nalleles = nalleles, groupnames = groupnames, groups = groups, group_prior = group_prior, group_inits = group_inits, 
                                     nchains = nchains, nits = nits, nburn = nburn, thin = thin, q_out = q_out, level = level)

dput(x = genetic_msa_Estimates, file = "mixname_genetic_msa_Estimates.txt")

setwd(orig_dir)

