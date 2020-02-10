# comethyl Soft Power Sampling ---------------------------------------------
# Charles Mordaunt
# 2/4/20

# Setup ####
setwd("/share/lasallelab/Yihui/MARBLES_WGBS_fastq/cytosine_reports")
.libPaths("/share/lasallelab/Charles/comethyl/R")
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
        fitIndices <- rbind(fitIndices, cbind("Regions" = i, sft$fitIndices))
}

# Analysis ####
# Load Data
setwd("~/Documents/Programming/MARBLES Placenta WGCNA")
fitIndices <- read.csv("Soft Power Data.csv", stringsAsFactors = FALSE)
table(fitIndices$slope < 0) # All TRUE, don't need to sign
fitIndices <- reshape2::melt(fitIndices, id.vars = c("Regions", "Power"))

# Plot fitIndices
gg <- ggplot(data = fitIndices) +
        geom_line(aes(x = Power, y = value, group = Regions, color = Regions)) +
        facet_wrap(vars(variable), nrow = 3, ncol = 3, scales = "free_y", strip.position = "left") +
        xlab("Soft Power Threshold") +
        scale_x_continuous(breaks = breaks_pretty(n = 4), expand = expand_scale(mult = c(0.05, 0.03))) +
        scale_y_continuous(breaks = breaks_pretty(n = 4)) +
        scale_color_gradient("Sampled\nRegions", breaks = breaks_pretty(n = 4)) +
        theme_bw(base_size = 24) +
        theme(axis.text = element_text(size = 14, color = "black"),
              axis.ticks = element_line(size = 1.25, color = "black"), axis.title.x = element_text(size = 18),
              axis.title.y = element_blank(), legend.background = element_blank(), legend.position = c(1.07, 0.89), 
              legend.title = element_text(size = 18), legend.text = element_text(size = 14), 
              panel.border = element_rect(color = "black", size = 1.25), 
              panel.grid = element_blank(), panel.spacing.x = unit(0.3, "lines"),
              panel.spacing.y = unit(0.8, "lines"), plot.margin = unit(c(1,7,0.7,0.2), "lines"),
              strip.background = element_blank(), strip.placement = "outside", 
              strip.switch.pad.wrap = unit(0, "lines"), strip.text = element_text(size = 18))
ggsave(filename = "Soft Power Sampling.pdf", plot = gg, dpi = 600, width = 13, height = 9, units = "in")

# Estimated Soft Power Threshold
fitTest <- subset(fitIndices, variable == "SFT.R.sq" & value > 0.8)
tapply(fitTest$Power, INDEX = fitTest$Regions, FUN = min)

# 40000  80000 120000 160000 200000 
#    18     18     16     17     16 

#' Conclusions
#' Fit differs somewhat with different number of regions sampled, but it's not that different
#' Estimated soft power threshold is within 1 or 2 of that with all regions
#' Connectivity measures decrease linearly with decreasing regions, possibly due to lower likelyhood of high correlations detected
#' Density and heterogeneity are identical, while centralization is different but has the same trend
#' Sampling can probably be used instead of analyzing entire set, but still need to try this with greater number of regions
#' 240k regions needs up to 500 GB of RAM






