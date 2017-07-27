## MSGB computation for the 20K patients in the WTCCC2 file
#
# Copyright 12/2014 Emilia Wieczorek ewieczorek@berkeley.edu


# initializing packages and contents --------------------------------------

rm(list = ls())
library(ggplot2)
library(ROCR)
source('./code/04_MSGBfunction.r')

# filename
WTCCC2_filename <- './outputs/WTCCC2_nonmissing.csv'
# read in the data
WTCCC2 <- read.csv(WTCCC2_filename, stringsAsFactors = FALSE, row.names=1)


# MSGB computation for 20K WTCCC2 patients --------------------------------

numPatients <- nrow(WTCCC2)
MSGBscores <- data.frame(phenotype = numeric(numPatients), MSGBscore = numeric(numPatients))
genotype <- rep(NA, ncol(WTCCC2)-6)
cat("[")
for (patient in 1:numPatients) {
  if (patient %% (numPatients%/%50) == 0) cat("-")
  genotype <- as.numeric(WTCCC2[patient,7:ncol(WTCCC2)]) # [1:length(genotype)] 
  names(genotype) <- names(WTCCC2)[7:ncol(WTCCC2)]
  score <- MSGB(genotype)
  MSGBscores$phenotype[patient] <- WTCCC2[patient, 6]
  MSGBscores$MSGBscore[patient] <- score
}
cat("]\n")


# plots -------------------------------------------------------------------

population <- factor(MSGBscores$phenotype, labels = c("Controls", "Cases"))
ggsave("outputs/MSGB_pop.pdf",
       w = 8, h = 6,
       qplot(x = MSGBscores$MSGBscore, 
             group = population, 
             color = population,
             fill = population, 
             alpha = I(0.4),
             geom = "density") + 
         theme_bw() + 
         labs(x = "MSGB", y="")
)


MSGBpred <- prediction(MSGBscores$MSGBscore, MSGBscores$phenotype)
MSGBROC <- performance(MSGBpred, "tpr", "fpr")
performance(MSGBpred, "auc")@y.values[[1]]
plot(MSGBROC)
