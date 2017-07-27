# MSGB computation script for 23andMe genome file
# input: 23andMe file
# output: MSGB score
#
# Copyright 12/2014 Emilia Wieczorek ewieczorek@berkeley.edu

# Intialize --------------------------------------------------------------

rm(list = ls())
source('/Users/Emilka/Documents/Research/23andMe-SNPs/R_analysis/code/04_MSGBfunction.r')

args <- commandArgs(trailingOnly = TRUE)
stopifnot (length(args) == 1)
file_23andMe <- args[1]


# Prepare data ------------------------------------------------------------

# read it in
genotypeData <- read.table(file_23andMe, stringsAsFactors = FALSE)
# get only the information that we need
genotypeData <- genotypeData[,c(1,2)]
colnames(genotypeData) <- c("rsID", "genotype")
# keeping only rsID for which we have coefficients
genotypeDataReduced <- genotypeData[genotypeData$rsID %in% coeff$rsID, ]
# retrieve the allele corresponding to the coefficients
alleleOfCoeffs <- coeff$allele[match(genotypeDataReduced$rsID, coeff$rsID)]
# count the mAF
genotype <- mapply(function(genotypeAlleles, alleleOfCoeff) 2-sum(genotypeAlleles %in% alleleOfCoeff), 
                   strsplit(genotypeDataReduced$genotype, split = ""), alleleOfCoeffs)

names(genotype) <- paste(genotypeDataReduced$rsID, alleleOfCoeffs, sep = "_")

MSGBscore <- MSGB(genotype)
cat(MSGBscore)

