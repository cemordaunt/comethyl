# Comethylated Analysis Pipeline ------------------------------------------
# Charles Mordaunt

# Load Packages ####
.libPaths("/share/lasallelab/Charles/comethylated/R")
sapply(c("tidyverse", "bsseq", "dmrseq", "WGCNA", "openxlsx"), require, character.only = TRUE)

# Functions ####
getCpGs <- function(colData, path = getwd(), pattern = "*CpG_report.txt.gz", 
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), BPPARAM = MulticoreParam(10), 
                    cov = 2, perSample = 0.75, save = TRUE, file = "Filtered_BSseq.rds", verbose = TRUE){
        if(verbose){
                message("[getCpGs] Loading CpG-level data")
        }
        files <- list.files(path, pattern = pattern) %>% .[pmatch(rownames(colData), table = .)]
        loci <- bsseq:::.readBismarkAsFWGRanges(files[1]) %>% chrSelectBSseq(seqnames = chroms)
        bs <- read.bismark(files, loci = loci, colData = colData, BPPARAM = BPPARAM, verbose = verbose)
        if(verbose){
                message("[getCpGs] Filtering CpG-level data for loci with at least ", cov, " reads in at least ", 
                        perSample * 100, "% of samples")
        }
        covSample <- (getCoverage(bs) >= cov) %>% DelayedMatrixStats::rowSums2()
        bs <- bs[covSample >= (perSample * ncol(bs)),]
        if(verbose){
                message("[getCpGs] Final BSseq Object:")
                print(bs)
        }
        if(save){
                if(verbose){
                        message("[getCpGs] Saving BSseq object as ", file)
                }
                saveRDS(bs, file = file)
        }
        if(verbose){
                message("[getCpGs] Complete!")
        }
        return(bs)
}

getRegions <- function(bs, maxGap = 150, n = 3, covMin = 10, methSD = 0.05, save = TRUE, file = "Filtered_Regions.txt",
                       verbose = TRUE){
        if(verbose){
                message("[getRegions] Calling regions when at least ", n, " CpGs are no more than ", maxGap, " bases apart")
        }
        regions <- bsseq:::regionFinder3(x = as.integer(rep(1, length(bs))), chr = as.character(seqnames(bs)), 
                                         positions = start(bs), maxGap = maxGap, verbose = FALSE)[["up"]]
        regions <- regions[regions$n >= n,]
        regions$chr <- as.character(regions$chr)
        regions$width <- regions$end - regions$start
        if(verbose){
                message("[getRegions] Filtering regions for at least ", covMin, 
                        " reads in all samples and methylation SD of at least ", methSD * 100, "%")
        }
        cov <- getCoverage(bs, regions = regions[,c("chr", "start", "end")], what = "perRegionTotal")
        regions$covMean <- DelayedMatrixStats::rowMeans2(cov)
        regions$covSD <- DelayedMatrixStats::rowSds(cov)
        regions$covMin <- DelayedArray::rowMins(cov)
        regions <- regions[regions$covMin >= covMin,]
        meth <- getMeth(bs, regions = regions[,c("chr", "start", "end")], type = "raw", what = "perRegion")
        regions$methMean <- DelayedMatrixStats::rowMeans2(meth)
        regions$methSD <- DelayedMatrixStats::rowSds(meth)
        regions <- regions[regions$methSD >= methSD,]
        regions$RegionID <- paste("Region", 1:nrow(regions), sep = "_")
        regions$width <- regions$end - regions$start
        regions <- regions[,c("RegionID", "chr", "start", "end", "width", "n", "covMin", "covMean", "covSD", "methMean", 
                              "methSD")]
        if(save){
                if(verbose){
                        message("[getRegions] Saving regions as ", file)
                }
                write.table(regions, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
        }
        if(verbose){
                message("[getRegions] Complete!")
        }
        return(regions)
}

# Read and Filter Bismark CpG Reports ###
colData <- read.xlsx("sample_info.xlsx", colNames = TRUE, rowNames = TRUE)
colData <- colData[4:9,] # subset for testing
bs <- getCpGs(colData)

# Call and Filter Regions ####
regions <- getRegions(bs)


