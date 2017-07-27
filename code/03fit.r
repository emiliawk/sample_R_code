# logistic regression script
# input: WTCCC2 patients data for selected SNPs out of the 196 SNPs first found in the intersection of all platforms
#
# Copyright: 12/2014 Emilia Wieczorek ewieczorek@berkeley.edu


# intializing packages and constants --------------------------------------

library(ggplot2)

rm(list = ls())
WTCCC2 <- read.csv('outputs/WTCCC2_nonmissing.csv', row.names = 1) 
bestSNPtable <- read.csv('outputs/SelectedRegionsData.csv', stringsAsFactors = F) [-1]

# functions ---------------------------------------------------------------

get_rsID <- function(listOfNames) {
  temp <- strsplit(listOfNames, "_")
  rsID <- sapply(temp, "[", 1)
  return(rsID)
}

get_allele <- function(listOfNames) {
  temp <- strsplit(listOfNames, "_")
  allele <- sapply(temp, "[", 2)
  return(allele)
}


# Prepare the risk alleles ------------------------------------------------

meanDoseControls <- apply(WTCCC2[ WTCCC2$PHENOTYPE == 1 ,7:ncol(WTCCC2)], 2, mean)
meanDoseCases <- apply(WTCCC2[ WTCCC2$PHENOTYPE == 2 ,7:ncol(WTCCC2)], 2, mean)  

qplot(meanDoseCases - meanDoseControls, geom = "density")

isRiskAllele <- as.logical(sign(meanDoseCases - meanDoseControls)+1)
WTCCC2[, (7:ncol(WTCCC2))[!isRiskAllele]] <- 2 - WTCCC2[, (7:ncol(WTCCC2))[!isRiskAllele]]

# fit ---------------------------------------------------------------------

WTCCC2$PHENOTYPE <- factor(WTCCC2$PHENOTYPE)
f <- as.formula(paste('PHENOTYPE ~', paste(colnames(WTCCC2[7:ncol(WTCCC2)]), collapse='+')))
# fit the model
fit <- glm(formula=f, data = WTCCC2, family = "binomial", maxit = 30)
# summary
summary(fit)
## odds ratios only
coeff <- summary(fit)$coefficients


# join data ---------------------------------------------------------------

# Transform coeff table:
coeff <- coeff[-1,] # remove intercept
coeffFinal <- data.frame(rsID = get_rsID(row.names(coeff)),
                         allele = get_allele(row.names(coeff)),
                         isRiskAllele = isRiskAllele,
                         coeff = coeff[,"Estimate"],
                         OR = exp(coeff[,"Estimate"]),
                         z.val = coeff[,"z value"],
                         p.val = coeff[,"Pr(>|z|)"],
                         stringsAsFactors = F)


# merge the two tables
# View(merge(bestSNPtable, coeffFinal, by.x = "BR_Proxy", by.y = "rsID"))
index <- match(bestSNPtable$BR_Proxy, coeffFinal$rsID) 
stopifnot(all(diff(index) == 1)) # Just to check that these are 1:X
bestSNPtableFinal <- cbind(bestSNPtable, coeffFinal)
  
ggplot(bestSNPtableFinal) + theme_bw() +
  geom_density(aes(x = log10(p.val)), binwidth = 0.1) 

ggplot(bestSNPtableFinal) + theme_bw() +
  geom_point(aes(x = log10(p.val), y=BR_RSquared)) +
  geom_smooth(aes(x = log10(p.val), y=BR_RSquared), method = lm) +
  scale_x_continuous(limits = c(-10,0))

## Cutoff on the p-values:
threshold.pval <- 0.01
b.index <- bestSNPtableFinal$p.val < threshold.pval
cat("## Removing", sum(!b.index), "Regions that have weak coefficients (p value >", threshold.pval, "): \n")
print(bestSNPtableFinal[!b.index, c("rsIDforJoseph", "rsID", "allele")])
outCoeffs <- bestSNPtableFinal[b.index, c("rsIDforJoseph", "rsID", "allele", "isRiskAllele", "OR", "p.val")]

#write to file
write.csv(file = 'outputs/SelectedRegionsDataWithCoeffs.csv', bestSNPtableFinal)
write.csv(file = 'outputs/MSGBCoeffsWithInfo.csv', outCoeffs)
write.csv(file = 'outputs/MSGBCoeffs.csv', outCoeffs[c("rsID", "allele", "isRiskAllele", "OR")])
