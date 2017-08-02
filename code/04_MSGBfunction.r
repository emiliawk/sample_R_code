## MSGB Computation Script
# Reads in the coefficient file and defines the function to compute msgb from a "genotype" named vector 
# that specifies the dosage found in the patient.
#
# Copyright 12/2014 Emilia Wieczorek ewieczorek@berkeley.edu

# Read in the coefficient file --------------------------------------------

coeff_filename <- '/Users/Emilka/Documents/Research/23andMe-SNPs/R_analysis/outputs/MSGBCoeffs.csv'
coeff <- read.csv(coeff_filename, stringsAsFactors = FALSE, row.names = 1)

# the function to compute MSGB --------------------------------------------

MSGB <- function(genotype) {
  # MSGB computes the MSGB corresponding to the gnotype passed as argument,
  # based on SNPs that are both present in the 23andMe data and the WTCC2.
  # - "genotype" is a named vector of minor allele counts with rsIDs as row names
  
  indexOfCoeffInGenotype <- match(rownames(coeff), names(genotype))
  
  if (any(is.na(indexOfCoeffInGenotype))) {
    stop("ARG_ERROR: the genotype provided doesn't have matching rsID_allele as names.")
  }
  
  multiplier <- log(coeff$OR)
  orderedGenotype <- genotype[indexOfCoeffInGenotype]
  burden <- ifelse( coeff$isRiskAllele, orderedGenotype, 2-orderedGenotype) * multiplier
  
  # compute and return
  return(sum(burden))
}


# test --------------------------------------------------------------------

if (test != FALSE) {
  coeff_old <- coeff
  coeff <- data.frame(row.names = c("rs1234_A", "rs2345_G"),
                      OR = c(0.8, 1.2),
                      isRiskAllele = c(FALSE, TRUE))
  #it should work
  MSGB(c(rs2345_G = 2, rs1234_A = 1))
  
  #it should not work
  tryCatch(
    MSGB(c(rs2345_A = 2, rs1234_A = 1)),
    error = function(e) {
      cat("## We got an expected error:", e$message, "\n")
    })
  
  coeff <- coeff_old
}
