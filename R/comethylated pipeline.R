# Comethylated Analysis Pipeline ------------------------------------------
# Charles Mordaunt
# 1/6/19

# Load Packages ####
.libPaths("/share/lasallelab/Charles/comethylated/R")
sapply(c("tidyverse", "bsseq", "dmrseq", "WGCNA", "openxlsx"), require, character.only = TRUE)

# Functions ####
getCpGs <- function(colData, path = getwd(), pattern = "*CpG_report.txt.gz", 
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), BPPARAM = MulticoreParam(10), 
                    cov = 2, perSample = 0.75, saveBSseq = TRUE, outfile = "Filtered_BSseq.rds", verbose = TRUE){
        if(verbose){
                message("[getCpGs] Loading CpG-level data")
        }
        files <- list.files(path, pattern = pattern) %>% .[pmatch(rownames(colData), table = .)]
        loci <- bsseq:::.readBismarkAsFWGRanges(files[1]) %>% chrSelectBSseq(seqnames = chroms)
        bs <- read.bismark(files, loci = loci, colData = colData, BPPARAM = BPPARAM, verbose = verbose)
        if(verbose){
                message("[getCpGs] Filtering CpG-level data for loci with at least ", cov, " reads in ", perSample * 100, 
                        "% of samples")
        }
        covSample <- (getCoverage(bs) >= cov) %>% DelayedMatrixStats::rowSums2()
        bs <- bs[covSample >= (perSample * ncol(bs)),]
        if(verbose){
                message("[getCpGs] Final BSseq Object:")
                print(bs)
        }
        if(saveBSseq){
                if(verbose){
                        message("[getCpGs] Saving BSseq Object as ", outfile)
                }
                saveRDS(bs, file = outfile)
        }
        if(verbose){
                message("[getCpGs] Complete!")
        }
        return(bs)
}

# Read and Filter Bismark CpG Reports ###
colData <- read.xlsx("sample_info.xlsx", colNames = TRUE, rowNames = TRUE)
colData <- colData[4:9,] # subset for testing
bs <- getCpGs(colData)

# Call and Filter Regions ####
regions <- bsseq:::regionFinder3(x = as.integer(rep(1, length(bs))), chr = as.character(seqnames(bs)), 
                                 positions = start(bs), maxGap = 150, verbose = FALSE)[["up"]] %>%
        subset(n >= 3, select = c("chr", "start", "end", "n"))
regions$chr <- as.character(regions$chr)
regions$width <- regions$end - regions$start
regions$minCov <- getCoverage(bs, regions = regions[,c("chr", "start", "end")], what = "perRegionTotal") %>% 
        DelayedArray::rowMins()
regions$SDs <- getMeth(bs, regions = regions[,c("chr", "start", "end")], type = "raw", what = "perRegion") %>% 
        DelayedMatrixStats::rowSds()
regions <- subset(regions, minCov >= 10 & SDs >= 0.05)
regions$RegionID <- paste("Region", 1:nrow(regions), sep = "_")
write.table(regions, file = "WGCNA_males/Males WGCNA minCov 10 SD 5 Filtered regions.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)




