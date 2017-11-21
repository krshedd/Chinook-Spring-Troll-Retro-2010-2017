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
# NOTE that fish from 113-95 and 113-97 are going to be included in D114 along with 112-65

# require(xlsx)
# spring_troll.df <- read.xlsx(file = "2010-2017 Spring troll asl by district_gh_ks.xlsx", sheetName = "original data", startRow = 1)
# str(spring_troll.df)

spring_troll.df <- read.table(file = "2010-2017 Spring troll asl by district_original_data.txt", sep = "\t", header = TRUE)
str(spring_troll.df)
table(spring_troll.df$Year, spring_troll.df$District)  # all samples

spring_troll.df$Sub.District.char <- sapply(as.character(spring_troll.df$Sub.District), function(i) {if(!is.na(i) & nchar(i) == 1) {paste0(0, i)} else {i} } )

spring_troll.df$Stat.Area <- paste0(spring_troll.df$District, spring_troll.df$Sub.District.char)
spring_troll.df$Stat.Area[is.na(spring_troll.df$District)] <- NA
table(spring_troll.df$Year, spring_troll.df$Stat.Area)  # all samples


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
  my.gcl$attributes$StatWeek <- spring_troll.df$Stat.Week[match.yr]
  my.gcl$attributes$Port <- spring_troll.df$Port.Code[match.yr]
  my.gcl$attributes$Quadrant <- spring_troll.df$Quadrant[match.yr]
  my.gcl$attributes$District <- spring_troll.df$District[match.yr]
  my.gcl$attributes$SubDistrict <- spring_troll.df$Sub.District.char[match.yr]
  my.gcl$attributes$StatArea <- spring_troll.df$Stat.Area[match.yr]
  my.gcl$attributes$Age <- spring_troll.df$Age.European[match.yr]
  my.gcl$attributes$LengthType <- spring_troll.df$Length.Type[match.yr]
  my.gcl$attributes$Length <- spring_troll.df$Length.Millimeters[match.yr]
  
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
dput(x = Spring10_17_Strata, file = "Objects/Spring10_17_Strata.txt")

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
invisible(sapply(paste0("KTROL", 10:17, "SP"), function(silly) {assign(x = paste0(silly, ".gcl"), value = dget(file = paste0("Raw genotypes/PooledCollections_PostQC/", silly, ".txt")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


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
table(unlist(ids.genotyped) %in% spring_troll.df$Dna.Specimen.No)  # all but District 108 fish
table(spring_troll.df$Dna.Specimen.No %in% unlist(ids.genotyped))  # not all fish from District 171-174 were genotyped (and passed data QC)

spring_troll.df$Genotyped <- spring_troll.df$Dna.Specimen.No %in% unlist(ids.genotyped)
addmargins(table(spring_troll.df$Year, spring_troll.df$District, spring_troll.df$Genotyped))

# options(java.parameters = "-Xmx100g")
# write.xlsx(x = spring_troll.df, file = "2010-2017 Spring troll asl by district.xlsx", sheetName = "D114com_genotyped", append = TRUE)
write.table(x = spring_troll.df, file = "2010-2017 Spring troll asl by district_original_data_genotyped.txt", row.names = FALSE)

# save.image("SpringTroll2010-2017.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Variable for Mixture ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(KTROL15SP.gcl$attributes$District.fac %in% 101:102)
table(KTROL15SP.gcl$attributes$District.fac %in% 103)
table(KTROL15SP.gcl$attributes$District.fac %in% 106:108 & KTROL15SP.gcl$attributes$StatArea != 10643)
table(KTROL15SP.gcl$attributes$District.fac %in% c(109:110, 112) & KTROL15SP.gcl$attributes$StatArea != 11265)
table(KTROL15SP.gcl$attributes$District.fac %in% 113)
table(KTROL15SP.gcl$attributes$District.fac %in% 114 | KTROL15SP.gcl$attributes$StatArea %in% c(11265, 11395, 11397))
table(KTROL15SP.gcl$attributes$District.fac %in% 183)


# Add new factor with mixtures
for(yr in 10:17){
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))

  my.gcl$attributes$Mixture <- NA
  my.gcl$attributes$Mixture[my.gcl$attributes$District.fac %in% 101:102] <- "101/102"
  my.gcl$attributes$Mixture[my.gcl$attributes$District.fac %in% 103] <- "103"
  my.gcl$attributes$Mixture[my.gcl$attributes$District.fac %in% 106:108 & my.gcl$attributes$StatArea != 10643] <- "106/107/108"
  my.gcl$attributes$Mixture[my.gcl$attributes$District.fac %in% c(109:110, 112) & my.gcl$attributes$StatArea != 11265] <- "109/110/112"
  my.gcl$attributes$Mixture[my.gcl$attributes$District.fac %in% 113] <- "113"
  my.gcl$attributes$Mixture[my.gcl$attributes$District.fac %in% 114 | my.gcl$attributes$StatArea %in% c(11265, 11395, 11397)] <- "114"
  my.gcl$attributes$Mixture[my.gcl$attributes$District.fac %in% 183] <- "183"
  
  my.gcl$attributes$Mixture <- factor(x = my.gcl$attributes$Mixture, levels = c("101/102", "103", "106/107/108", "109/110/112", "113", "114", "183"))
  
  assign(x = paste0("KTROL", yr, "SP.gcl"), value = my.gcl)
}


# Pivot of years by mixtures with NA
addmargins(t(sapply(as.character(10:17), function(yr) {
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  table(my.gcl$attributes$Mixture, useNA = "always")
} )))


# Pivot of years by mixtures without NA
addmargins(t(sapply(as.character(10:17), function(yr) {
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  table(my.gcl$attributes$Mixture)
} )))

#     101/102 103 106/107/108 109/110/112  113 114  183  Sum
# 10      148   0          10          16  296 242    0  712
# 11      152   0         193          70  459 111    0  985
# 12      128   0         257         132  301 168    0  986
# 13      133   0         105          97  251 117  497 1200
# 14      142  96         126          94  209  69  377 1113
# 15      109 100         156         112  177  24  316  994
# 16       51 102         135         140  104  78   95  705
# 17       50 101         125          90  187  61    0  614
# Sum     913 399        1107         751 1984 870 1285 7309


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save PostQC_Metadata .gcl's as back-up:
# dir.create("Raw genotypes/PooledCollections_PostQC_Metadata")
invisible(sapply(Spring10_17_Strata, function(silly) {
  dput(x = get(paste(silly, ".gcl", sep = '')), file = paste0("Raw genotypes/PooledCollections_PostQC_Metadata/" , silly, ".txt"))
} )); beep(8)



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
invisible(sapply(paste0("KTROL", 10:17, "SP"), function(silly) {assign(x = paste0(silly, ".gcl"), value = dget(file = paste0("Raw genotypes/PooledCollections_PostQC_Metadata/", silly, ".txt")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Change rownames for scores and counts to FK_FISH_ID ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## All of the .gcl functions that rely on IDs look in the rownames for scores,
## FK_FISH_ID and rownames(scores) + rownames(counts) need to match!!!

for(yr in 10:17){
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  
  rownames(my.gcl$scores) <- rownames(my.gcl$counts) <- as.character(my.gcl$attributes$FK_FISH_ID)
  
  assign(x = paste0("KTROL", yr, "SP.gcl"), value = my.gcl)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save PostQC_Metadata .gcl's as back-up:
# dir.create("Raw genotypes/PooledCollections_PostQC_Metadata_Rename")
invisible(sapply(Spring10_17_Strata, function(silly) {
  dput(x = get(paste(silly, ".gcl", sep = '')), file = paste0("Raw genotypes/PooledCollections_PostQC_Metadata_Rename/" , silly, ".txt"))
} )); beep(8)




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
invisible(sapply(paste0("KTROL", 10:17, "SP"), function(silly) {assign(x = paste0(silly, ".gcl"), value = dget(file = paste0("Raw genotypes/PooledCollections_PostQC_Metadata_Rename/", silly, ".txt")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mixtures <- levels(KTROL10SP.gcl$attributes$Mixture)
dput(x = mixtures, file = "Objects/mixtures.txt")

mixtures.names <- setNames(object = c("D101102Troll", "D103Troll", "D106107108Troll", "D109110112Troll", "D113Troll", "D114Troll", "D183Troll"), 
                           nm = mixtures)
dput(x = mixtures.names, file = "Objects/mixtures.names.txt")

# dir.create("BAYES")
# dir.create("BAYES/Mixture")

# Loop over years and mixtures
for(yr in 10:17){
  my.gcl <- get(paste0("KTROL", yr, "SP.gcl"))
  
  for(mix in mixtures) {
    IDs <- list("my" = na.omit(AttributesToIDs.GCL(silly = "my", attribute = "Mixture", matching = mix)))
    if(length(IDs[["my"]])) {
      invisible(CreateMixture.GCL(sillys = "my", loci = GAPSLoci_reordered, IDs = IDs, 
                                  mixname = paste0(mixtures.names[mix], "_20", yr), 
                                  dir = "BAYES/Mixture/", type = "BAYES", PT = FALSE))
    }  # if IDS
  }  # mixture within year
  
}  # year


# Double check mixture files
sum(sapply(list.files(path = "BAYES/Mixture/", full.names = TRUE), function(fle) {nrow(read.table(file = fle, header = FALSE))} ))

# Remove mixtures we do not intend to run!
setwd("BAYES/Mixture/")
unlink(x = c("D106107108Troll_2010.mix"))
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/Spring Troll 2010-2017")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("BAYES/Baseline")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/BAYES/Baseline/GAPS357pops13loci.bse", to = "BAYES/Baseline/")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupVec26RG_357.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupNames26.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/SEAKPops357.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GAPS357PopsInits.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/WASSIPSockeyeSeeds.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/mixfortran.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/bayesfortran_357.txt", to = "Objects")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("BAYES/Control")

SEAKobjects <- list.files(path = "Objects", recursive = FALSE)
invisible(sapply(SEAKobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); rm(SEAKobjects); beep(2)


# Flat Pop prior
GAPS357PopFlatPrior <- Prior.GCL(groupvec = 1:357, groupweights = rep(1/357, 357), minval = 0.001)
dput(x = GAPS357PopFlatPrior, file = "Objects/GAPS357PopFlatPrior.txt")

# Dump Control Files
GroupVec26RG_357 <- dget(file = "Objects/GroupVec26RG_357.txt")
SEAKPops357 <- dget(file = "Objects/SEAKPops357.txt")


all.mixtures <- sapply(list.files(path = "BAYES/Mixture"), function(mix) {unlist(strsplit(x = mix, split = ".mix"))[1]}, USE.NAMES = FALSE)
dput(x = all.mixtures, file = "Objects/all.mixtures.txt")

sapply(all.mixtures, function(Mix) {
  invisible(CreateControlFile.GCL(sillyvec = SEAKPops357, loci = GAPSLoci_reordered, mixname = Mix, basename = "GAPS357pops13loci", suffix = "", nreps = 40000, nchains = 5,
                                  groupvec = GroupVec26RG_357, priorvec = GAPS357PopFlatPrior, initmat = GAPS357PopsInits, dir = "BAYES/Control",
                                  seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = mixfortran, basefortran = bayesfortran_357, switches = "F T F T T T F"))
} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("BAYES/Output")
sapply(all.mixtures, function(Mix) {
  invisible(dir.create(path = paste0("BAYES/Output/", Mix)))
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Compare BAYES and genetic_msa ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BAYES_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:26, groupnames = groupnames, maindir = "BAYES/Output/", mixvec = mixname,
                                                prior = '', ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

genetic_msa_Estimates

require(gplots)
barplot2(height = rbind(BAYES_Estimates$D114Troll_2010[, "mean"], genetic_msa_Estimates$mean), beside = TRUE)  # very bad
legend("topright", legend = c("BAYES", "genetic_msa"), fill = c("red", "yellow"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize BAYES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Spring10_17_26RG_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = 1:26, groupnames = GroupNames26, maindir = "BAYES/Output", mixvec = all.mixtures,
                               prior = '', ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)
# dir.create("Estimates objects")
dput(x = Spring10_17_26RG_Estimates, file = "Estimates objects/Spring10_17_26RG_Estimates.txt")
sapply(Spring10_17_26RG_Estimates, function(mix) {table(mix[, "GR"] > 1.2)})

file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupNames4.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupNames4Pub.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupVec4.txt", to = "Objects")

Spring10_17_4RG_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec4, groupnames = GroupNames4Pub, maindir = "BAYES/Output", mixvec = all.mixtures,
                               prior = '', ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)
dput(x = Spring10_17_4RG_Estimates, file = "Estimates objects/Spring10_17_4RG_Estimates.txt")
Spring10_17_4RG_Estimates <- dget(file = "Estimates objects/Spring10_17_4RG_Estimates.txt")

sapply(Spring10_17_4RG_Estimates, function(mix) {table(mix[, "GR"] > 1.2)})
sapply(Spring10_17_4RG_Estimates, function(mix) {mix[c("Alaska", "TBR"), "mean"]})

sapply(Spring10_17_4RG_Estimates, function(mix) {mix[, "mean"]} )
sapply(Spring10_17_4RG_Estimates, function(mix) {mix[, "sd"] / mix[, "mean"]} )  # CV
sapply(Spring10_17_4RG_Estimates, function(mix) {sum(mix[, "sd"] / mix[, "mean"] < 0.20)} )  # how many RGs with CV < 20%
sapply(Spring10_17_4RG_Estimates, function(mix) {mix[, "95%"] - mix[, "5%"]} )  # 90% CI range


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Read in Harvest and Sample Size Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(xlsx)
harvest.df <- read.xlsx(file = "2010-2017 Spring troll asl by district_gh_ks.xlsx", sheetName = "CE000522", startRow = 23, header = TRUE)
str(harvest.df)

harvest.df$Mixture <- NA
harvest.df$Mixture[harvest.df$District %in% 101:102] <- "101/102"
harvest.df$Mixture[harvest.df$District %in% 103] <- "103"
harvest.df$Mixture[harvest.df$District %in% 106:108 & harvest.df$Area.Value != 10643] <- "106/107/108"
harvest.df$Mixture[harvest.df$District %in% c(109:110, 112) & harvest.df$Area.Value != 11265] <- "109/110/112"
harvest.df$Mixture[harvest.df$District %in% 113] <- "113"
harvest.df$Mixture[harvest.df$District %in% 114 | harvest.df$Area.Value %in% c(11265, 11395, 11397)] <- "114"
harvest.df$Mixture[harvest.df$District %in% 183] <- "183"

harvest.df$Mixture <- factor(x = harvest.df$Mixture, levels = c("101/102", "103", "106/107/108", "109/110/112", "113", "114", "183"))

dput(x = harvest.df, file = "Objects/harvest.df.txt")

require(reshape)
harvest_mix.df <- aggregate(N.Catch ~ Year + Mixture, data = harvest.df, sum)
str(harvest_mix.df)
cast(harvest_mix.df, Year ~ Mixture, value = "N.Catch")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sample sizes
all.mixtures.samplesize <- sapply(all.mixtures, function(mix) {dim(read.table(file = paste0("BAYES/Mixture/", mix, ".mix")))[1]} )
dput(x = all.mixtures.samplesize, file = "Objects/all.mixtures.samplesize.txt")

mixtures.names
mixtures.names2 <- names(mixtures.names)
names(mixtures.names2) <- mixtures.names

mixtures.df <- as.data.frame(t(sapply(all.mixtures, function(mix) {unlist(strsplit(x = mix, split = "_"))} )), stringsAsFactors = FALSE)
names(mixtures.df) <- c("Mixname", "Year")
mixtures.df$Year <- as.numeric(mixtures.df$Year)
mixtures.df$Full.Mixname <- all.mixtures
mixtures.df$Mixture <- factor(x = mixtures.names2[mixtures.df$Mixname], levels = levels(harvest_mix.df$Mixture))
mixtures.df$n <- all.mixtures.samplesize

dput(x = mixtures.df, file = "Objects/mixtures.df.txt")


str(mixtures.df)

all.mixtures.n100 <- names(which(all.mixtures.samplesize >= 100))
dput(x = all.mixtures.n100, file = "Objects/all.mixtures.n100.txt")

# Subset data for n >= 100
Spring10_17_4RG_Estimates_n100 <- Spring10_17_4RG_Estimates[all.mixtures.n100]
dput(x = Spring10_17_4RG_Estimates_n100, file = "Estimates objects/Spring10_17_4RG_Estimates_n100.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge harvest and sample size
spring.estimates.df <- merge(x = harvest_mix.df, y = mixtures.df, by = c("Mixture", "Year"), all = FALSE)
spring.estimates.df$Alaska.mean.p <- sapply(Spring10_17_4RG_Estimates, function(mix) {mix["Alaska", "mean"]})
spring.estimates.df$TBR.mean.p <- sapply(Spring10_17_4RG_Estimates, function(mix) {mix["TBR", "mean"]})
spring.estimates.df$Alaska.mean.C <- spring.estimates.df$Alaska.mean.p * spring.estimates.df$N.Catch
spring.estimates.df$TBR.mean.C <- spring.estimates.df$TBR.mean.p * spring.estimates.df$N.Catch

round(cast(spring.estimates.df, Year ~ Mixture, value = "Alaska.mean.C"))
round(cast(spring.estimates.df, Year ~ Mixture, value = "TBR.mean.C"))
round(cast(spring.estimates.df, Year ~ Mixture, value = "n"))

dput(x = spring.estimates.df, file = "Objects/spring.estimates.df.txt")

harvest <- setNames(object = spring.estimates.df$N.Catch, nm = spring.estimates.df$Full.Mixname)
dput(x = harvest, file = "Objects/harvest.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Subset to only mixtures with >= 100 fish
spring.estimates.n100.df <- spring.estimates.df
spring.estimates.n100.df[spring.estimates.n100.df$n < 100, c("Alaska.mean.p", "Alaska.mean.C", "TBR.mean.p", "TBR.mean.C")] <- NA

# Heatmap of total Catch
require(lattice)
new.colors <- colorRampPalette(c("white", "darkgreen"))
data.mat <- as.matrix(round(cast(spring.estimates.n100.df, Year ~ Mixture, value = "N.Catch")))
# data.mat[is.na(data.mat)] <- 0
levelplot(data.mat, 
          col.regions = new.colors, 
          at = seq(from = 0, to = max(data.mat, na.rm = TRUE), length.out = 100), 
          main = "Total Catch", xlab = "Year", ylab = "District Area", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill",
          panel = function(...) {
            panel.fill("black")
            panel.levelplot(...)}
)  # aspect = "iso" will make squares


# Heatmap of mean Alaska
require(lattice)
new.colors <- colorRampPalette(c("white", "darkblue"))
data.mat <- as.matrix(cast(spring.estimates.n100.df, Year ~ Mixture, value = "Alaska.mean.p")) * 100
levelplot(data.mat, 
          col.regions = new.colors, 
          at = seq(from = 0, to = 100, length.out = 100), 
          main = "Mean Alaska %", xlab = "Year", ylab = "District Area", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill", 
          panel = function(...) {
            panel.fill("black")
            panel.levelplot(...)}
)  # aspect = "iso" will make squares

# Heatmap of mean TBR %
require(lattice)
new.colors <- colorRampPalette(c("white", "darkblue"))
data.mat <- as.matrix(cast(spring.estimates.n100.df, Year ~ Mixture, value = "TBR.mean.p")) * 100
levelplot(data.mat, 
          col.regions = new.colors, 
          at = seq(from = 0, to = 100, length.out = 100), 
          main = "Mean TBR %", xlab = "Year", ylab = "District Area", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill", 
          panel = function(...) {
            panel.fill("black")
            panel.levelplot(...)}
)  # aspect = "iso" will make squares


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create 4RG Summary Tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# dir.create("Estimates tables")
require(xlsx)

EstimatesStats <- Spring10_17_4RG_Estimates_n100
SampSizes <- all.mixtures.samplesize
HarvestVec <- harvest
PubNames <- setNames(object = paste("Spring Troll", spring.estimates.df$Year, "District(s)", spring.estimates.df$Mixture),
                     nm = spring.estimates.df$Full.Mixname)

for(mix in all.mixtures.n100) {
  
  TableX <- matrix(data = "", nrow = 7, ncol = 7)
  TableX[1, 1] <- paste0(PubNames[mix], " (n=", SampSizes[mix], ", catch=", formatC(x = HarvestVec[mix], digits = 0, big.mark = ",", format = "f"), ")")
  TableX[2, 6] <- "90% CI"
  TableX[3, 2:7] <- c("Reporting Group", "Mean", "SD", "Median", "5%", "95%")
  TableX[4:7, 1] <- 1:4
  TableX[4:7, 2] <- rownames(EstimatesStats[[mix]])
  TableX[4:7, 3:7] <- formatC(x = EstimatesStats[[mix]][, c("mean", "sd", "median", "5%", "95%")], digits = 3, format = "f")
  
  write.xlsx(x = TableX, file = "Estimates tables/SpringTroll2017_4RG_Estimates.xlsx",
             col.names = FALSE, row.names = FALSE, sheetName = paste(mix, " 4RG"), append = TRUE)
  
}


# save.image("SpringTroll2010-2017.RData")