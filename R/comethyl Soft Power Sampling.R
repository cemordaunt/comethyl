# comethyl Soft Power Sampling ---------------------------------------------
# Charles Mordaunt
# 2/4/20

# Setup ####
setwd("/share/lasallelab/Yihui/MARBLES_WGBS_fastq/cytosine_reports")
.libPaths("/share/lasallelab/Charles/comethylated/R")
sapply(c("scales", "openxlsx", "rlist", "tidyverse", "ggdendro", "cowplot", "bsseq", "dmrseq", "WGCNA", "sva"), require, 
       character.only = TRUE)
options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
disableWGCNAThreads()
set.seed(5)

# Functions ####
getSoftPower <- function(meth, powerVector = 1:20, corType = c("pearson", "bicor"), maxPOutliers = 0.1, 
                         RsquaredCut = 0.8, blockSize = 40000, gcInterval = blockSize - 1, verbose = TRUE){
        corType <- match.arg(corType)
        if(verbose){
                message("[getSoftPower] Analyzing scale-free topology with ", corType, 
                        " correlation to determine best soft-thresholding power")
                verboseNum <- 10
        } else {
                verboseNum <- 0
        }
        if(corType == "pearson"){
                sft <- pickSoftThreshold(meth, RsquaredCut = RsquaredCut, powerVector = powerVector, networkType = "signed",
                                         moreNetworkConcepts = TRUE, corFnc = "cor", blockSize = blockSize, 
                                         gcInterval = gcInterval, verbose = verboseNum)
        } else {
                if(corType == "bicor"){
                        sft <- pickSoftThreshold(meth, RsquaredCut = RsquaredCut, powerVector = powerVector, 
                                                 networkType = "signed", moreNetworkConcepts = TRUE, corFnc = "bicor",
                                                 corOptions = list(maxPOutliers = maxPOutliers), blockSize = blockSize, 
                                                 gcInterval = gcInterval, verbose = verboseNum)
                } else {
                        stop("[getSoftPower] Error: corType must be either pearson or bicor")
                }
        }
        if(is.na(sft$powerEstimate)){
                fit <- -sign(sft$fitIndices[,"slope"]) * sft$fitIndices[,"SFT.R.sq"]
                sft$powerEstimate <- sft$fitIndices$Power[fit == max(fit)]
        }
        if(verbose){
                message("[getSoftPower] At soft power threshold = ", sft$powerEstimate, 
                        ", fit = ", round(sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == sft$powerEstimate], 3), 
                        " and mean connectivity = ", round(sft$fitIndices$mean.k.[sft$fitIndices$Power == sft$powerEstimate], 1))
        }
        return(sft)
}

# Estimate Soft Power from Different Region Subsets ####
methAdj <- readRDS("WGCNA_MARBLES_placenta/Adjusted_Region_Methylation_Male.rds")
dim(methAdj) # 62 311678
fitIndices <- NULL
totalRegions <- ncol(methAdj)

for(i in seq(40000, to = 320000, by = 40000)){
        i <- min(i, totalRegions)
        message("Sampling ", i, " regions from a total of ", totalRegions)
        subRegions <- sample(totalRegions, size = i)
        methAdjSub <- methAdj[,subRegions]
        sft <- getSoftPower(methAdjSub, powerVector = 12:19, blockSize = 40000)
        fitIndices <- rbind(fitIndices, cbind("sample" = i, sft$fitIndices))
}

# Running on Epigenerate