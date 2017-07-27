## SNP selection script
# this script selects data for logistic regression analysis
# input: WTCCC2 patient's data for 196 SNPs found in the intersection
# some data is missing therefore we select best SNPs for which not more than 5% of data is missing 
# we remove patients for which more than 5% of data is missing for the SNPs previously selected
#
# Copyright: 12/2014 Emilia Wieczorek ewieczorek@berkeley.edu

rm(list=ls())
library(plyr)
library(ggplot2)

# filenames and paths -----------------------------------------------------

WTCCC2_filename <- 'outputs/WTCC2.csv' # patient data for 196 SNPs in the final intersection
SNPtable_filename <- 'outputs/finalIntersection.csv'


# missing data handling ---------------------------------------------------

# read in the data
WTCCC2 <- read.csv(WTCCC2_filename, stringsAsFactors = FALSE)
SNPtable <- read.csv(SNPtable_filename, stringsAsFactors = FALSE)
SNPtable$X <- NULL

# check how much data is missing
missing_rate <- sapply(WTCCC2[8:ncol(WTCCC2)], function(x) mean(is.na(x)))

# add missing_rate info to SNP table
oldNames <- names(WTCCC2[8:ncol(WTCCC2)])
newNames <- sapply(oldNames, function(x) sub("_.*", "", x))
indexing <- match(SNPtable$Proxy, newNames)
MissingRate <- unname(missing_rate[indexing])
SNPtable <- cbind(SNPtable, MissingRate)


# Have a look at the missing rate / rsquared relation
ggplot(SNPtable) + geom_point(aes(x = MissingRate, y= RSquared, group = rsIDforJoseph, colour = rsIDforJoseph)) +
  geom_line(aes(x = MissingRate, y= RSquared, group = rsIDforJoseph, colour = rsIDforJoseph)) +
  scale_x_log10() +
  theme_bw()

# First missing rate filtering
threshold <- 0.05
SNPtableFiltered <- SNPtable[SNPtable$MissingRate<threshold,]
regions <- unique(SNPtable$rsIDforJoseph)
removedRegions <- regions[! regions %in% unique(SNPtableFiltered$rsIDforJoseph)]
cat("### Warning: With a threshold of", threshold, "on the missing rate, the", length(removedRegions), "following regions have been removed:\n", removedRegions, "\n" )

# from the SNPs found in the final intersection select the "best one"
bestSNPtable <- ddply(SNPtableFiltered, ~ rsIDforJoseph, function(dfi) {
  num_proxies = nrow(dfi)
  ## Best Rsquared within top missing rate
  bestMissing <- dfi[ dfi$MissingRate == min(dfi$MissingRate),c("Proxy", "MissingRate", "RSquared", "Order")]
  if(nrow(bestMissing) > 1) {
    bestMissing <- bestMissing[which.max(bestMissing$RSquared),]
  }
  names(bestMissing) <- paste0("BMR_", names(bestMissing))
  ## Best missing within top RSquared (SNP table is ordered by RSquared)
  bestRSquared <- dfi[dfi$RSquared == max(dfi$RSquared),c("Proxy", "MissingRate", "RSquared", "Order")]
  if(nrow(bestRSquared) > 1) {
    bestRSquared <- bestRSquared[which.min(bestRSquared$MissingRate),]
  }
  names(bestRSquared) <- paste0("BR_", names(bestRSquared))
  return(data.frame(num_proxies, bestRSquared, bestMissing))
})

ggplot(bestSNPtable) + theme_bw() +
  geom_density(aes(x = BR_MissingRate), fill = "red", alpha = 0.4) +
  geom_density(aes(x = BMR_MissingRate), fill = "blue", alpha = 0.4)

ggplot(bestSNPtable) + theme_bw() +
  geom_density(aes(x = BR_RSquared), fill = "red", alpha = 0.4) +
  geom_density(aes(x = BMR_RSquared), fill = "blue", alpha = 0.4)

ggplot(bestSNPtable) + theme_bw() +
  geom_point(aes(x = BR_MissingRate, y = BR_RSquared), fill = "red", alpha = 0.4, shape = 21) +
  geom_point(aes(x = BMR_MissingRate, y = BMR_RSquared), fill = "blue", alpha = 0.4, shape = 21) 
  

cols_to_keep <- match(bestSNPtable$BR_Proxy,newNames)
# remember to move index by 7, first seven cols of WTCCC2 are: X,IID,FID,PAT,MAT,SEX,PHENOTYPE
cols_to_keep <- cols_to_keep + 7
WTCCC2_reduced <- WTCCC2[,c(2:7,cols_to_keep)]

# keep only those subjects that have no missing data for bestSNPs
WTCCC2_final <- na.omit(WTCCC2_reduced)

# write to file
write.csv(file = file.path("outputs", "SelectedRegionsData.csv"), bestSNPtable)
write.csv(file = file.path("outputs", 'WTCCC2_nonmissing.csv'), WTCCC2_final)

