#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### SEAK Chinook Spring Troll 2010-2017 ####
# Kyle Shedd Tue Nov 07 14:08:51 2017
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
# The goal of this script is to revisit Chinook salmon mixtures from the SEAK
# commercial spring troll harvests from 2010-2017 looking at D14 using the GAPS3.0
# baseline containing 357 populations in 26 reporting groups characterized by 
# 13 uSATs. All mixtures are to be analyzed with the program BAYES.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Specific Objectives ####
# This script will:
# 1) Import mixture data
# 2) Add attribute data
# 3) Define spatio-temporal strata
# 4) Perform a data QC on mixtures
# 5) Prepare BAYES input files
# 6) Summarize BAYES results
# 7) Generate plots and tables of results

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/Spring Troll 2010-2017")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
username <- "krshedd"
password <- "********"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
Spring2010Mixtures <- c("KSPRING10H", "KSPRING10J", "KSPRING10K", "KSPRING10P", "KSPRING10S", "KSPRING10W")
Spring2011Mixtures <- c("KSPRING11J", "KSPRING11K", "KSPRING11P", "KSPRING11S", "KSPRING11W")
Spring2012Mixtures <- c("KSPRING12J", "KSPRING12K", "KSPRING12P", "KSPRING12S", "KSPRING12W")  # Stikine and Taku directed fishery samples never extracted "KTROL12SR" "KTROL12TR"
Spring2013Mixtures <- c("KSPRING13J", "KSPRING13K", "KSPRING13P", "KSPRING13S", "KSPRING13W", "KSPRING13Y")
Spring2014Mixtures <- c("KSPRING14C", "KSPRING14J", "KSPRING14K", "KSPRING14P", "KSPRING14S", "KSPRING14W", "KSPRING14Y")
Spring2015Mixtures <- c("KSPRING15C", "KSPRING15J", "KSPRING15K", "KSPRING15P", "KSPRING15S", "KSPRING15W", "KSPRING15Y")
Spring2016Mixtures <- c("KTROL16SP")  # "KTROL16D8" not used, no extractions
Spring2017Mixtures <- c("KTROL17SP")

## Pull genotypes
LOKI2R_GAPS.GCL(sillyvec = unlist(sapply(objects(pattern = "Spring"), get)), username = username, password = password)


## Save unaltered .gcls
# dir.create("Raw genotypes")
# dir.create("Raw genotypes/OriginalCollections")
invisible(sapply(unlist(sapply(objects(pattern = "Spring"), get)), function(silly) {dput(x = get(paste0(silly, ".gcl")), file = paste0("Raw genotypes/OriginalCollections/" , silly, ".txt"))} )); beep(8)

# dir.create("Objects")
dput(x = LocusControl, file = "Objects/LocusControl.txt")
invisible(sapply(objects(pattern = "Mixtures"), function(mix) {dput(x = get(mix), file = paste0("Objects/", mix, ".txt"))}))
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GAPSLoci_reordered.txt", to = "Objects")
GAPSLoci_reordered <- dget(file = "Objects/GAPSLoci_reordered.txt")

dimnames(KTROL16SP.gcl$counts)[[2]]
GAPSLoci_reordered

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pool into a single silly per year
PoolCollections.GCL(collections = Spring2010Mixtures, loci = GAPSLoci_reordered, newname = "KTROL10SP")
PoolCollections.GCL(collections = Spring2011Mixtures, loci = GAPSLoci_reordered, newname = "KTROL11SP")
PoolCollections.GCL(collections = Spring2012Mixtures, loci = GAPSLoci_reordered, newname = "KTROL12SP")
PoolCollections.GCL(collections = Spring2013Mixtures, loci = GAPSLoci_reordered, newname = "KTROL13SP")
PoolCollections.GCL(collections = Spring2014Mixtures, loci = GAPSLoci_reordered, newname = "KTROL14SP")
PoolCollections.GCL(collections = Spring2015Mixtures, loci = GAPSLoci_reordered, newname = "KTROL15SP")
PoolCollections.GCL(collections = Spring2016Mixtures, loci = GAPSLoci_reordered, newname = "KTROL16SP")
PoolCollections.GCL(collections = Spring2017Mixtures, loci = GAPSLoci_reordered, newname = "KTROL17SP")

dimnames(KTROL16SP.gcl$counts)[[2]]

sapply(paste0("KTROL", 10:17, "SP"), function(silly) {get(paste0(silly, ".gcl"))$n} )
# KTROL10SP KTROL11SP KTROL12SP KTROL13SP KTROL14SP KTROL15SP KTROL16SP KTROL17SP 
# 1106      1260      1072      1427      1162      1105      1114      1010 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Change FK_FISH_ID to the back end of SillySource
str(KTROL10SP.gcl$attributes$FK_FISH_ID)
str(KTROL16SP.gcl$attributes$FK_FISH_ID)

KTROL10SP.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL10SP.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL11SP.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL11SP.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL12SP.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL12SP.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL13SP.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL13SP.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL14SP.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL14SP.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL15SP.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL15SP.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL16SP.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL16SP.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL17SP.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL17SP.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )


# dir.create("Raw genotypes/PooledCollections")
invisible(sapply(paste0("KTROL", 10:17, "SP"), function(silly) {dput(x = get(paste0(silly, ".gcl")), file = paste0("Raw genotypes/PooledCollections/" , silly, ".txt"))} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/Spring Troll 2010-2017/")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
SEAKobjects <- list.files(path = "Objects", recursive = FALSE)
# SEAKobjects <- SEAKobjects[-which(SEAKobjects == "Vials" | SEAKobjects == "OLD_BAD_LOCUSCONTROL")]
SEAKobjects

invisible(sapply(SEAKobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

## Get un-altered mixtures
invisible(sapply(paste0("KTROL", 10:17, "SP"), function(silly) {assign(x = paste0(silly, ".gcl"), value = dget(file = paste0("Raw genotypes/PooledCollections/", silly, ".txt")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pair with metadata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pair with district data from Anne
# NOTE that fish from 113-95 and 113-97 got moved to 114-21
require(xlsx)

spring_troll.df <- read.xlsx(file = "2010-2017 Spring troll asl by district.xlsx", sheetName = "D114com", startRow = 2)
str(spring_troll.df)

ids <- sapply(paste0("KTROL", 10:17, "SP"), function(silly) {get(paste0(silly, ".gcl"))$attributes$FK_FISH_ID} )
str(ids)

# Are we missing metadata for fish we have genotyped?
table(unlist(ids) %in% spring_troll.df$Dna.Specimen.No)
# FALSE  TRUE 
# 252    9004 

# Which years are missing metadata
table(sapply(names(unlist(ids)[!unlist(ids) %in% spring_troll.df$Dna.Specimen.No]), function(id) {
  unlist(strsplit(x = unlist(strsplit(x = id, split = "KTROL"))[2], split = "SP"))[1]
} ))
# 10  11  12  13  15  16 
#  1   2   6   1   1 241

# Paste the missing fish into ASL .csv to see what project they are from
writeClipboard(as.character(unlist(ids)[!unlist(ids) %in% spring_troll.df$Dna.Specimen.No]))
# All of the 2016 fish are from "District 108 Spring Troll" project
# E-mailed Anne to see if it is safe to assume that they were all caught in 108


# Match data up by year and look at district breakdowns
for(yr in 10:17){
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  match.yr <- match(my.gcl$attributes$FK_FISH_ID, spring_troll.df$Dna.Specimen.No)
  # table(spring_troll.df$Year[match.yr])
  my.gcl$attributes$District <- spring_troll.df$District.[match.yr]
  my.gcl$attributes$StatWeek <- spring_troll.df$Stat.Week[match.yr]
  
  assign(x = paste0("match.20", yr), value = match.yr)
  assign(x = paste0("KTROL", yr, "SP.gcl"), value = my.gcl)
}

sapply(as.character(10:17), function(yr) {
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  addmargins(table(my.gcl$attributes$District, useNA = "always"))
} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

Spring10_17_Strata <- paste0("KTROL", 10:17, "SP")

Spring10_17_Strata_SampleSizes <- matrix(data = NA, nrow = length(Spring10_17_Strata), ncol = 4, 
                                       dimnames = list(Spring10_17_Strata, c("Genotyped", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_Spring10_17_Strata_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = Spring10_17_Strata, loci = GAPSLoci_reordered)
min(Original_Spring10_17_Strata_SampleSizebyLocus)  ## 991
apply(Original_Spring10_17_Strata_SampleSizebyLocus, 1, min) / apply(Original_Spring10_17_Strata_SampleSizebyLocus, 1, max)  ## Good, 0.928

Original_Spring10_17_Strata_PercentbyLocus <- apply(Original_Spring10_17_Strata_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_Spring10_17_Strata_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_Spring10_17_Strata_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares

#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_Spring10_17_Strata_ColSize <- sapply(paste0(Spring10_17_Strata, ".gcl"), function(x) get(x)$n)
Spring10_17_Strata_SampleSizes[, "Genotyped"] <- Original_Spring10_17_Strata_ColSize

### Missing
## Remove individuals with >20% missing data
Spring10_17_Strata_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = Spring10_17_Strata, proportion = 0.8)
dput(x = Spring10_17_Strata_MissLoci, file = "Objects/Spring10_17_Strata_MissLoci.txt")

## Get number of individuals per silly after removing missing loci individuals
ColSize_Spring10_17_Strata_PostMissLoci <- sapply(paste0(Spring10_17_Strata, ".gcl"), function(x) get(x)$n)
Spring10_17_Strata_SampleSizes[, "Missing"] <- Original_Spring10_17_Strata_ColSize - ColSize_Spring10_17_Strata_PostMissLoci

### Duplicate
## Check within collections for duplicate individuals.
Spring10_17_Strata_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = Spring10_17_Strata, loci = GAPSLoci_reordered, quantile = NULL, minproportion = 0.95)
Spring10_17_Strata_DuplicateCheckReportSummary <- sapply(Spring10_17_Strata, function(x) Spring10_17_Strata_DuplicateCheck95MinProportion[[x]]$report)
Spring10_17_Strata_DuplicateCheckReportSummary
dput(x = Spring10_17_Strata_DuplicateCheckReportSummary, file = "Objects/Spring10_17_Strata_DuplicateCheckReportSummary.txt")

## Remove duplicate individuals
Spring10_17_Strata_RemovedDups <- RemoveDups.GCL(Spring10_17_Strata_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_Spring10_17_Strata_PostDuplicate <- sapply(paste0(Spring10_17_Strata, ".gcl"), function(x) get(x)$n)
Spring10_17_Strata_SampleSizes[, "Duplicate"] <- ColSize_Spring10_17_Strata_PostMissLoci-ColSize_Spring10_17_Strata_PostDuplicate

### Final
Spring10_17_Strata_SampleSizes[, "Final"] <- ColSize_Spring10_17_Strata_PostDuplicate
Spring10_17_Strata_SampleSizes

dput(x = Spring10_17_Strata_SampleSizes, file = "Objects/Spring10_17_Strata_SampleSizes.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save PostQC .gcl's as back-up:
# dir.create("Raw genotypes/PooledCollections_PostQC")
invisible(sapply(Spring10_17_Strata, function(silly) {
  dput(x = get(paste(silly, ".gcl", sep = '')), file = paste0("Raw genotypes/PooledCollections_PostQC/" , silly, ".txt"))
} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Tables by District ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pivots to show number of fish by district by year that we have genotypes for

t(sapply(as.character(10:17), function(yr) {
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  addmargins(table(my.gcl$attributes$District, useNA = "always"))
} ))

# Anne's data only has quadrant projects
table(spring_troll.df$Year, spring_troll.df$Project)

# Read in OceanAK ASL data to add District info to D108 and D111 projects
spring_troll_oceanAK.df <- read.csv(file = "Harvest - Detailed ASL Samples 2010-2017.csv")
str(spring_troll_oceanAK.df)

# Which years?
table(spring_troll_oceanAK.df$ï..Year, spring_troll_oceanAK.df$Project)  # 2012 and 2016

# Fish IDs
ids_D108 <- spring_troll_oceanAK.df$Dna.Specimen.No[spring_troll_oceanAK.df$Project == "District 108 Spring Troll"]
ids_D111 <- spring_troll_oceanAK.df$Dna.Specimen.No[spring_troll_oceanAK.df$Project == "District 111 Spring Troll"]


table(KTROL12SP.gcl$attributes$FK_FISH_ID %in% ids_D108)
table(KTROL12SP.gcl$attributes$FK_FISH_ID %in% ids_D111)
table(KTROL16SP.gcl$attributes$FK_FISH_ID %in% ids_D108)

# Create character vector of district
for(yr in 10:17){
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  my.gcl$attributes$District.chr <- as.character(my.gcl$attributes$District)

  assign(x = paste0("KTROL", yr, "SP.gcl"), value = my.gcl)
}

# Add 108 data from 2016
table(KTROL16SP.gcl$attributes$District)
table(KTROL16SP.gcl$attributes$District.chr)

match.oceanAK.2016 <- match(KTROL16SP.gcl$attributes$FK_FISH_ID[KTROL16SP.gcl$attributes$FK_FISH_ID %in% ids_D108], 
                            spring_troll_oceanAK.df$Dna.Specimen.No)
KTROL16SP.gcl$attributes$District.chr[KTROL16SP.gcl$attributes$FK_FISH_ID %in% ids_D108] <- as.character(spring_troll_oceanAK.df$District[match.oceanAK.2016])

table(KTROL16SP.gcl$attributes$District.chr)

levels(KTROL10SP.gcl$attributes$District)


# Create new factor with all districts
for(yr in 10:17){
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  my.gcl$attributes$District.fac <- factor(x = my.gcl$attributes$District.chr, levels = c(" ", as.character(101:115), "183"))
  
  assign(x = paste0("KTROL", yr, "SP.gcl"), value = my.gcl)
}


# Pivot of years by district
addmargins(t(sapply(as.character(10:17), function(yr) {
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  table(my.gcl$attributes$District.fac, useNA = "always")
} )))

KTROL17SP.gcl$attributes$FK_FISH_ID[is.na(KTROL17SP.gcl$attributes$District.fac)]

# Add column in Anne's sheet to denote which samples have been genotyped
ids.genotyped <- sapply(as.character(10:17), function(yr) {
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  my.gcl$attributes$FK_FISH_ID
})

str(ids.genotyped)
table(spring_troll.df$Dna.Specimen.No %in% unlist(ids.genotyped))

spring_troll.df$Genotyped <- spring_troll.df$Dna.Specimen.No %in% unlist(ids.genotyped)
addmargins(table(spring_troll.df$Year, spring_troll.df$District., spring_troll.df$Genotyped))

# options(java.parameters = "-Xmx100g")
write.xlsx(x = spring_troll.df, file = "2010-2017 Spring troll asl by district.xlsx", sheetName = "D114com_genotyped", append = TRUE)
write.table(x = spring_troll.df, file = "2010-2017 Spring troll asl by district.txt", row.names = FALSE)

# save.image("SpringTroll2010-2017.RData")
