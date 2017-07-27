# File path & names -------------------------------------------------------

WTCCC2SNPsFolder <- "data/WTCCC2 data/SNP_lists"
MSGB_SNPsandProxies <- "data/SNPsProxy.txt"
platformSNPs23andMe <- "data/SNPPlatformdata.Rdata"


# Extraction Functions ----------------------------------------------------

getWTCCC2snps <- function(SNPfileFolder) {
  # gets SNPs available on the WTCC2 platform, there is one file for each chromosome
  
  # chromosome list
  filenames <- file.path(SNPfileFolder, paste('chr.',1:22,'.All.bim', sep = ''))
  SNPs <- c()
  
  for (filename in filenames) {
    tmp <- read.table(filename, stringsAsFactors=FALSE)
    # subset to only include data for the current chromosome (this ammounts to removing first few rows)
    data <- subset(tmp, tmp$V1 %in% as.character(1:22))
    # add the SNPs to the list of all SNPs available on WTCCC2 platform
    SNPs <- c(SNPs, data$V2)
  }
  return(SNPs)
}

get23andMeSNPs <- function(filename) {
  # load the 23anaMe platform data
  load(filename)
  # drop the Lumigenix since we are not using it
  compare_rsIDs$Lumigenix <- NULL
  # cols now are "rsID" "V2_1" "V2_2" "V2_3" "V2_4" "V2_5" "V3"   "V3_1" "V4"  
  # omit those that have at least one 'NA' (get the intersection of snps available on all 23andMe platforms)
  data <- na.omit(compare_rsIDs)
  # get a list of SNPs available on all 23andMe platforms (convert from factor to character)
  SNPs <- as.character(data$rsID)
  return(SNPs)
}

liftoverMask <- function(liftover_filename ) {
  liftoverData <- read.csv(liftover_filename, stringsAsFactors=FALSE)
  # get data only for those that match
  data <- subset(liftoverData, liftoverData$Match == 'YES')
  # now get the SNPs
  SNPs <- data$rsID
  return(SNPs)
}

getProxySNPs <- function(filename) {
  SNPtable <- read.table(filename, stringsAsFactors = FALSE, header = T)
  # get a smaller subset of the data, only those that are of interest
  data <- subset(SNPtable, select=c(rsIDforJoseph,Proxy,Distance,RSquared,DPrime,Order,MAF,Coordinate_HG17, Coordinate_HG18, Chromosome))
  return(data)
}



# Script ------------------------------------------------------------------

# get the list of SNPs
SNPs_23andMe <- get23andMeSNPs(platformSNPs23andMe)
#SNPs_liftover <- liftoverMask(liftover_filename)
SNPs_WTCCC2 <- getWTCCC2snps(WTCCC2SNPsFolder)
# this is a data frame 
SNPs_proxy <- getProxySNPs(MSGB_SNPsandProxies)

# get the intersections
liftover_23andMe_intersection <- SNPs_23andMe #intersect(SNPs_23andMe, SNPs_liftover)
# get SNPs available on all 23andMe platforms AND WTCCC2 platform
WTCCC2_23andMe_intersection <- intersect(liftover_23andMe_intersection, SNPs_WTCCC2)
# finally get the intersection with the proxy SNPs
# get the list of proxy SNPs
proxies <- SNPs_proxy$Proxy
# get the intersection of all
WTCCC2_23andMe_proxy_intersection <- intersect(WTCCC2_23andMe_intersection, proxies)

# kepp only the data for SNPs in the final intersection
SNPtableFinal <- subset(SNPs_proxy, SNPs_proxy$Proxy %in% WTCCC2_23andMe_proxy_intersection)
# write to file
write.csv(file = file.path("outputs", 'finalIntersection.csv'), SNPtableFinal)

# some info to look at
cat("## number of SNPs available on 23andMe platforms:", length(SNPs_23andMe), "\n")
# print(length(liftover_23andMe_intersection)) # number of 23andMe SNPs after doing the liftover
cat("## number of SNPs available on 23andMe and WTCCC2 platforms:", length(WTCCC2_23andMe_intersection), "\n") 
cat("## number of SNPs available in the intersection with proxy SNPs:", length( WTCCC2_23andMe_proxy_intersection), "\n")
cat("## number of unique SNPs:",length(unique(SNPtableFinal$rsIDforJoseph)), "\n")
cat("## number of original SNPs that are associated with MS:",nrow(subset(SNPtableFinal, Order == 1)), "\n")

library(plyr)
SNPtableReduced <- ddply(SNPtableFinal, ~ rsIDforJoseph, function(dfi) {
  return(dfi[dfi$Order == min(dfi$Order),])
})

print(sum(SNPtableReduced$RSquared == 1))
print(sum(SNPtableReduced$DPrime == 1))

library(ggplot2)
ggplot(SNPtableFinal) + geom_line(aes(x = Order, y= RSquared, group = rsIDforJoseph, colour = rsIDforJoseph)) +
  theme_bw()

# create a file with a list of SNPs for PLINK
write(file = "outputs/SNPsForPLINK.txt", SNPtableFinal$Proxy)
