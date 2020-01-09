# Comethylated Analysis Pipeline ------------------------------------------
# Charles Mordaunt

# Load Packages ####
.libPaths("/share/lasallelab/Charles/comethylated/R")
sapply(c("scales", "openxlsx", "tidyverse", "bsseq", "dmrseq", "WGCNA", "sva"), require, character.only = TRUE)

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

plotRegionStats <- function(regions, bins = 100, histCol = "#132B43", lineCol = "red", nBreaks = 4, save = TRUE, 
                            file = "Region_Plots.pdf", width = 11, height = 8.5, verbose = TRUE){
        if(verbose){
                message("[plotRegionStats] Plotting histograms of region statistics")
        }
        regions <- reshape2::melt(regions[,c("RegionID", "width", "n", "covMin", "covMean", "methMean", "methSD")], 
                                  id.vars = "RegionID")
        medians <- aggregate(value ~ variable, data = regions, FUN = median)
        gg <- ggplot(data = regions)
        gg <- gg +
                geom_histogram(aes(x = value), bins = bins, fill = histCol, color = histCol) +
                geom_vline(data = medians, aes(xintercept = value), color = lineCol) +
                facet_wrap(vars(variable), nrow = 2, ncol = 3, scales = "free") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(expand = expand_scale(mult = c(0.008, 0.05))) +
                theme_bw(base_size = 24) +
                theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
                      legend.position = "none", strip.text.x = element_text(size = 16), 
                      axis.ticks.x = element_line(size = 1.25, color = "black"), axis.ticks.y = element_blank(), 
                      axis.title = element_blank(), strip.background = element_blank(), 
                      plot.margin = unit(c(0,1,1,0.4), "lines"), panel.spacing.y = unit(0, "lines"), 
                      axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_blank())
        if(save){
                if(verbose){
                        message("[plotRegionStats] Saving plots as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        if(verbose){
                message("[plotRegionStats] Complete!")
        }
        return(gg)
}

plotSDstats <- function(regions, bins = 100, nBreaks = 4, legend.position = c(1.085,0.938), save = TRUE, 
                        file = "SD_Plots.pdf", width = 8.5, height = 8.5, verbose = TRUE){
        if(verbose){
                message("[plotSDstats] Plotting methylation SD vs region statistics")
        }
        regions <- reshape2::melt(regions[,c("RegionID", "n", "covMin", "covMean", "methMean", "methSD")], 
                                  id.vars = c("RegionID", "methSD"))
        gg <- ggplot(data = regions)
        gg <- gg +
                geom_bin2d(aes(x = value, y = methSD, color = ..count..), bins = bins) +
                facet_wrap(vars(variable), nrow = 2, ncol = 2, scales = "free_x") +
                scale_fill_continuous(name = "Count", trans = "log10") +
                scale_color_continuous(guide = FALSE, trans = "log10") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks), expand = expand_scale(mult = c(0.0062, 0.05))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks), expand = expand_scale(mult = c(0.006, 0.05))) +
                theme_bw(base_size = 24) +
                theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
                      legend.position = legend.position, legend.background = element_blank(),
                      legend.title = element_text(size = 14), legend.text = element_text(size = 12),
                      strip.text.x = element_text(size = 16), 
                      axis.ticks = element_line(size = 1.25, color = "black"), axis.title.x = element_blank(),
                      axis.title.y = element_text(size = 14, color = "black"), 
                      strip.background = element_blank(), plot.margin = unit(c(0.5,6,1,1), "lines"), 
                      panel.spacing.y = unit(0, "lines"), axis.text = element_text(size = 12, color = "black"))
        if(save){
                if(verbose){
                        message("[plotSDstats] Saving plots as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        if(verbose){
                message("[plotSDstats] Complete!")
        }
        return(gg)
}

getRegionMeth <- function(regions, bs, type = "raw", save = TRUE, file = "Region_Methylation.rds", verbose = TRUE){
        if(verbose){
                message("[getRegionMeth] Calculating region methylation from BSseq object")
        }
        meth <- getMeth(bs, regions = regions[,c("chr", "start", "end")], type = type, what = "perRegion")
        rownames(meth) <- regions$RegionID
        if(save){
                if(verbose){
                        message("[getRegionMeth] Saving region methylation as ", file)
                }
                saveRDS(meth, file = file)
        }
        if(verbose){
                message("[getRegionMeth] Complete!")
        }
        return(meth)
}

adjustRegionMeth <- function(meth, mod = matrix(1, nrow = ncol(meth), ncol = 1), save = TRUE, 
                             file = "Adjusted_Region_Methylation.rds", verbose = TRUE){
        if(verbose){
                message("[adjustRegionMeth] Determining number of principal components to adjust for")
        }
        n.pc <- num.sv(meth, mod = mod, seed = 5)
        if(verbose){
                message("[adjustRegionMeth] Adjusting region methylation for the top ", n.pc, " principal components")
        }
        methAdj <- sva_network(meth, n.pc = n.pc) %>% t()
        if(save){
                if(verbose){
                        message("[adjustRegionMeth] Saving adjusted region methylation as ", file)
                }
                saveRDS(methAdj, file = file)
        }
        if(verbose){
                message("[adjustRegionMeth] Complete!")
        }
        return(methAdj)
}

# Read and Filter Bismark CpG Reports ####
colData <- read.xlsx("sample_info.xlsx", colNames = TRUE, rowNames = TRUE)
colData <- colData[4:9,] # subset for testing
bs <- getCpGs(colData)

# Call and Filter Regions ####
regions <- getRegions(bs)
plotRegionStats(regions)
plotSDstats(regions)

# Adjust Methylation Data for PCs ####
meth <- getRegionMeth(regions, bs = bs)
mod <- model.matrix(~1, data = pData(bs))
methAdj <- adjustRegionMeth(meth, mod = mod)

sampleTree <- (1 - cor(t(methAdj))) %>% as.dist() %>% hclust(method = "average")
pdf(file = "Sample_Dendrogram.pdf", height = 5, width = 11)
par(mar = c(0, 5, 1, 0))
plot(sampleTree, main = "", sub = "", xlab = "", cex = 0.8, cex.lab = 1, cex.axis = 1)
invisible(dev.off())







