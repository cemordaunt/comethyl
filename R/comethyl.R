# Comethyl Analysis Pipeline ------------------------------------------
# Charles Mordaunt

# Load Packages ####
.libPaths("/share/lasallelab/Charles/comethylated/R")
sapply(c("scales", "openxlsx", "tidyverse", "ggdendro", "bsseq", "dmrseq", "WGCNA", "sva"), require, character.only = TRUE)

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
                        message("[getCpGs] Saving file as ", file)
                }
                saveRDS(bs, file = file)
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
        regions$covMean <- DelayedMatrixStats::rowMeans2(cov) %>% round(digits = 5)
        regions$covSD <- DelayedMatrixStats::rowSds(cov) %>% round(digits = 5)
        regions$covMin <- DelayedArray::rowMins(cov)
        regions <- regions[regions$covMin >= covMin,]
        meth <- getMeth(bs, regions = regions[,c("chr", "start", "end")], type = "raw", what = "perRegion")
        regions$methMean <- DelayedMatrixStats::rowMeans2(meth, na.rm = TRUE) %>% round(digits = 5)
        regions$methSD <- DelayedMatrixStats::rowSds(meth, na.rm = TRUE) %>% round(digits = 5)
        regions <- regions[regions$methSD >= methSD,]
        regions$RegionID <- paste("Region", 1:nrow(regions), sep = "_")
        regions$width <- regions$end - regions$start
        regions <- regions[,c("RegionID", "chr", "start", "end", "width", "n", "covMin", "covMean", "covSD", "methMean", 
                              "methSD")]
        if(save){
                if(verbose){
                        message("[getRegions] Saving file as ", file)
                }
                write.table(regions, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
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
                facet_wrap(vars(variable), nrow = 2, ncol = 3, scales = "free", strip.position = "bottom") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(expand = expand_scale(mult = c(0.008, 0.05))) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_blank(),
                      axis.ticks.x = element_line(size = 1.25, color = "black"), axis.ticks.y = element_blank(), 
                      axis.title = element_blank(), legend.position = "none", 
                      panel.border = element_rect(color = "black", size = 1.25), 
                      panel.grid = element_blank(), panel.spacing.x = unit(0.6, "lines"),
                      panel.spacing.y = unit(0.8, "lines"), plot.margin = unit(c(1,1,0.3,0.6), "lines"),
                      strip.background = element_blank(), strip.placement = "outside", 
                      strip.switch.pad.wrap = unit(0, "lines"), strip.text.x = element_text(size = 16))
        if(save){
                if(verbose){
                        message("[plotRegionStats] Saving plots as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

plotSDstats <- function(regions, bins = 100, nBreaks = 4, legend.position = c(1.09,0.9), save = TRUE, 
                        file = "SD_Plots.pdf", width = 8.5, height = 8.5, verbose = TRUE){
        if(verbose){
                message("[plotSDstats] Plotting methylation SD vs region statistics")
        }
        regions <- reshape2::melt(regions[,c("RegionID", "n", "covMin", "covMean", "methMean", "methSD")], 
                                  id.vars = c("RegionID", "methSD"))
        gg <- ggplot(data = regions)
        gg <- gg +
                geom_bin2d(aes(x = value, y = methSD, color = ..count..), bins = bins) +
                facet_wrap(vars(variable), nrow = 2, ncol = 2, scales = "free_x", strip.position = "bottom") +
                scale_fill_continuous(name = "Count", trans = "log10") +
                scale_color_continuous(guide = FALSE, trans = "log10") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks), expand = expand_scale(mult = c(0.0062, 0.05))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks), expand = expand_scale(mult = c(0.006, 0.05))) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 12, color = "black"), 
                      axis.ticks = element_line(size = 1.25, color = "black"),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = 16, color = "black"),
                      legend.background = element_blank(), legend.position = legend.position, 
                      legend.text = element_text(size = 12), legend.title = element_text(size = 16), 
                      panel.border = element_rect(color = "black", size = 1.25), panel.grid = element_blank(), 
                      panel.spacing.x = unit(1, "lines"), panel.spacing.y = unit(0.8, "lines"), 
                      plot.margin = unit(c(1,6,0.3,1), "lines"), strip.background = element_blank(), 
                      strip.placement = "outside", strip.switch.pad.wrap = unit(0, "lines"), 
                      strip.text.x = element_text(size = 16))
        if(save){
                if(verbose){
                        message("[plotSDstats] Saving plots as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
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
                        message("[getRegionMeth] Saving file as ", file)
                }
                saveRDS(meth, file = file)
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
                        message("[adjustRegionMeth] Saving file as ", file)
                }
                saveRDS(methAdj, file = file)
        }
        return(methAdj)
}

getDendro <- function(x, transpose = FALSE, distance = c("euclidean", "pearson", "bicor"), maxPOutliers = 0.1, 
                      verbose = TRUE){
        if(transpose){
                if(verbose){
                        message("[getDendro] Transposing data")
                }
                x <- t(x)
        }
        distance <- match.arg(distance)
        if(distance == "euclidean"){
                if(verbose){
                        message("[getDendro] Clustering with euclidean distance")
                        dist <- dist(x)
                } else {
                        if(distance == "pearson"){
                                if(verbose){
                                        message("[getDendro] Clustering with pearson correlation as the distance")
                                }
                                dist <- (1 - cor(x)) %>% as.dist()
                        } else {
                                if(distance == "bicor"){
                                        if(verbose){
                                                message("[getDendro] Clustering with bicor correlation as the distance")
                                        }
                                        dist <- (1 - bicor(x, maxPOutliers = maxPOutliers)) %>% as.dist()
                                } else {
                                        stop("[getDendro] Error: Distance must be either euclidean, pearson, or bicor")
                                }
                        }
                }
        }
        dendro <- hclust(dist, method = "average")
        return(dendro)
}

plotDendro <- function(dendro, label = TRUE, labelSize = 2.5, expandX = c(2,2), expandY = c(7,2), nBreaks = 4, save = TRUE,
                       file = "Dendrogram.pdf", width = 11, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotDendro] Plotting dendrogram")
        }
        dendroPlot <- dendro_data(dendro)
        fix <- dendroPlot$segments$yend == 0
        dendroPlot$segments$yend[fix] <- dendroPlot$segments$y[fix] - max(dendroPlot$segments$y) * 0.05
        dendroPlot$labels$y <- dendroPlot$segments$yend[fix] - max(dendroPlot$segments$y) * 0.01
        gg <- ggplot()
        gg <- gg +
                geom_segment(data = dendroPlot$segments, aes(x = x, y = y, xend = xend, yend = yend), lwd = 0.3) +
                scale_x_continuous(expand = expand_scale(add = expandX)) +
                scale_y_continuous(expand = expand_scale(add = expandY), breaks = breaks_pretty(n = nBreaks)) +
                ylab("Height") +
                theme_dendro() +
                theme(plot.margin = unit(c(1,1,0,1), "lines"), 
                      panel.background = element_rect(color = "black", size = 1.1),
                      axis.ticks.y = element_line(), axis.text.y = element_text(size = 12, color = "black"), 
                      axis.title.y = element_text(size = 16, angle = 90, vjust = 2))
        if(label){
                gg <- gg +
                        geom_text(data = dendroPlot$labels, aes(x = x, y = y, label = label), angle = 90, hjust = 1, 
                                  size = labelSize)
        }
        if(save){
                if(verbose){
                        message("[plotDendro] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

getSoftPower <- function(meth, powerVector = 1:20, corType = c("pearson", "bicor"), maxPOutliers = 0.1, 
                         RsquaredCut = 0.8, verbose = TRUE){
        corType <- match.arg(corType)
        if(verbose){
                message("[getSoftPower] Analyzing scale-free topology with ", corType, 
                        " correlation to determine best soft-thresholding power")
                verboseNum <- 10
        } else {
                verboseNum <- 0
        }
        if(corType == "pearson"){
                sft <- pickSoftThreshold(meth, RsquaredCut = RsquaredCut, powerVector = powerVector, 
                                         networkType = "signed", corFnc = "cor", blockSize = 40000, verbose = verboseNum)
        } else {
                if(corType == "bicor"){
                        sft <- pickSoftThreshold(meth, RsquaredCut = RsquaredCut, powerVector = powerVector, 
                                                 networkType = "signed", corFnc = "bicor", 
                                                 corOptions = list(maxPOutliers = maxPOutliers), blockSize = 40000, 
                                                 verbose = verboseNum)
                } else {
                        stop("[getSoftPower] Error: corType must be either pearson or bicor")
                }
        }
        if(is.na(sft$powerEstimate)){
                sft$powerEstimate <- sft$fitIndices$Power[sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq)]
        }
        if(verbose){
                message("[getSoftPower] At soft power threshold = ", sft$powerEstimate, 
                        ", fit = ", round(sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == sft$powerEstimate], 3), 
                        " and mean connectivity = ", round(sft$fitIndices$mean.k.[sft$fitIndices$Power == sft$powerEstimate], 1))
        }
        return(sft)
}

plotSoftPower <- function(sft, pointCol = "#132B43", lineCol = "red", nBreaks = 4, save = TRUE, 
                          file = "Soft_Power_Plots.pdf", width = 8.5, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotSoftPower] Plotting scale-free topology fit and mean connectivity by soft power threshold")
        }
        fitIndices <- data.frame(power = sft$fitIndices$Power, 
                                 fit = -sign(sft$fitIndices[,"slope"]) * sft$fitIndices[,"SFT.R.sq"],
                                 log10_meanConnectivity = log10(sft$fitIndices$mean.k.),
                                 powerEstimate = sft$powerEstimate) %>% 
                reshape2::melt(id.vars = c("power", "powerEstimate"))
        gg <- ggplot(data = fitIndices)
        gg <- gg +
                geom_vline(aes(xintercept = powerEstimate), color = lineCol) +
                geom_text(aes(x = powerEstimate, y = 0, label = powerEstimate), color = lineCol, nudge_x = -1) +
                geom_point(aes(x = power, y = value), color = pointCol, size = 1.2) +
                facet_wrap(vars(variable), nrow = 1, ncol = 2, scales = "free_y", strip.position = "left") +
                xlab("Soft Power Threshold") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                expand_limits(x = 0, y = c(0,1)) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 12, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"), axis.title.x = element_text(size = 16),
                      axis.title.y = element_blank(), legend.position = "none", 
                      panel.border = element_rect(color = "black", size = 1.25), 
                      panel.grid = element_blank(), panel.spacing.x = unit(0.3, "lines"),
                      panel.spacing.y = unit(0.8, "lines"), plot.margin = unit(c(1,1,0.7,0.2), "lines"),
                      strip.background = element_blank(), strip.placement = "outside", 
                      strip.switch.pad.wrap = unit(0, "lines"), strip.text.x = element_text(size = 16))
        if(save){
                if(verbose){
                        message("[plotSoftPower] Saving plots as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

getModules <- function(methAdj, power = NULL, corType = c("pearson", "bicor"), deepSplit = 4, minModuleSize = 10, 
                       mergeCutHeight = 0.1, nThreads = 6, save = TRUE, file = "Modules.rds", verbose = TRUE){
        if(is.null(power)){
                stop("[getModules] You must select a soft power threshold")
        }
        if(verbose){
                message("[getModules] Constructing network and detecting modules in blocks")
                verboseNum <- 10
        } else {
                verboseNum <- 0
        }
        corType <- match.arg(corType)
        modules <- blockwiseModules(methAdj, checkMissingData = FALSE, maxBlockSize = 40000, corType = corType, 
                                    power = power, networkType = "signed", TOMtype = "signed", deepSplit = deepSplit, 
                                    minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, nThreads = nThreads,
                                    verbose = verboseNum)
        if(save){
                if(verbose){
                        message("[getModules] Saving modules as ", file)
                }
                saveRDS(modules, file = file)
        }
        return(modules)
}

plotRegionDendro <- function(modules, save = TRUE, file = "Region_Dendrograms.pdf", width = 11, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotRegionDendro] Plotting region dendrograms and modules for each block")
        }
        blockColors <- lapply(modules$blockGenes, function(x) modules$colors[x])
        blockNames <- paste("Block ", 1:length(modules$dendrograms), " (", sapply(modules$blockGenes, length), " regions)", sep = "")
        if(save){
                if(verbose){
                        message("[plotRegionDendro] Saving plot as ", file)
                }
                pdf(file = file, width = width, height = height)
        }
        invisible(mapply(FUN = plotDendroAndColors, dendro = modules$dendrograms, colors = blockColors, main = blockNames, 
                         MoreArgs = list(groupLabels = "Modules", dendroLabels = FALSE, marAll = c(1.5,5,3,1.5), saveMar = FALSE, 
                                         cex.lab = 1.2, cex.colorLabels = 1.2, autoColorHeight = FALSE, lwd = 0.8, 
                                         colorHeight = 0.15, cex.axis = 1, frame.plot = TRUE)))
        invisible(dev.off())
}

.regionFilterTotals <- function(regions, covMin, methSD){
        regions <- regions[regions$covMin >= covMin & regions$methSD >= methSD,]
        totals <- c("covMin" = covMin, "methSD" = methSD, "totalRegions_K" = nrow(regions)/10^3, 
                    "totalWidth_Mb" = sum(regions$width)/10^6, "totalN_M" = sum(regions$n)/10^6)
        return(totals)
}

getRegionFilterTotals <- function(regions, covMin = rep(seq(0,20,2), each = 11), methSD = rep(seq(0,0.1,0.01), 11),
                                  save = TRUE, file = "Region_Filter_Totals.txt", verbose = TRUE){
        if(verbose){
                message("[getRegionFilterTotals] Calculating region totals at specified covMin and methSD cutoffs")
        }
        regionFilterTotals <- mapply(FUN = .regionFilterTotals, covMin = covMin, methSD = methSD, 
                                     MoreArgs = list(regions = regions)) %>% t() %>% as.data.frame()
        if(save){
                if(verbose){
                        message("[getRegionFilterTotals] Saving file as ", file)
                        write.table(regionFilterTotals, file = file, sep = "\t", row.names = FALSE)
                }
        }
        return(regionFilterTotals)
}

plotRegionFilterTotals <- function(regionFilterTotals, nBreaks = 4, legend.position = c(1.11,0.925), save = TRUE, 
                                   file = "Region_Filter_Totals.pdf", width = 11, height = 11, verbose = TRUE){
        if(verbose){
                message("[plotRegionFilterTotals] Plotting region filter totals")
        }
        regionFilterTotals <- reshape2::melt(regionFilterTotals, id.vars = c("covMin", "methSD"))
        regionFilterTotals$variable <- as.character(regionFilterTotals$variable) %>% 
                str_replace_all(pattern = c("totalRegions_K" = "Total Regions (Thousands)", "totalWidth_Mb" = "Total Width (Mb)", 
                                            "totalN_M" = "Total CpGs (Millions)")) %>%
                factor(levels = c("Total Regions (Thousands)", "Total Width (Mb)", "Total CpGs (Millions)"))
        gg <- ggplot(data = regionFilterTotals)
        gg <- gg +
                geom_line(aes(x = covMin, y = value, group = methSD, color = methSD)) +
                geom_text(data = subset(regionFilterTotals, covMin == min(covMin)), 
                          aes(x = covMin, y = value, group = methSD, color = methSD, label = methSD), 
                          size = 4.5, check_overlap = TRUE, nudge_x = -0.2, hjust = 1) +
                facet_wrap(vars(variable), nrow = 3, ncol = 1, scales = "free_y", strip.position = "left") +
                xlab("Minimum Coverage") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks), expand = expand_scale(mult = c(0.07, 0.03))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_color_gradient("Minimum SD", breaks = breaks_pretty(n = nBreaks - 1)) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 14, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"), axis.title.x = element_text(size = 18),
                      axis.title.y = element_blank(), legend.background = element_blank(), legend.position = legend.position, 
                      legend.title = element_text(size = 18), legend.text = element_text(size = 14), 
                      panel.border = element_rect(color = "black", size = 1.25), 
                      panel.grid = element_blank(), panel.spacing.x = unit(0.3, "lines"),
                      panel.spacing.y = unit(0.8, "lines"), plot.margin = unit(c(1,9,0.7,0.2), "lines"),
                      strip.background = element_blank(), strip.placement = "outside", 
                      strip.switch.pad.wrap = unit(0, "lines"), strip.text = element_text(size = 18))
        if(save){
                if(verbose){
                        message("[plotRegionFilterTotals] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

# Set Global Options ####
options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
enableWGCNAThreads(nThreads = 6)

# Read and Filter Bismark CpG Reports ####
colData <- read.xlsx("sample_info.xlsx", rowNames = TRUE)
bs <- getCpGs(colData)

# Call Regions without Filtering ####
regions_unf <- getRegions(bs, covMin = 0, methSD = 0, file = "Unfiltered_Regions.txt")
plotRegionStats(regions_unf, file = "Unfiltered_Region_Plots.pdf")
plotSDstats(regions_unf, file = "Unfiltered_SD_Plots.pdf")

# Examine Region Totals at Different Cutoffs ####
regionFilterTotals <- getRegionFilterTotals(regions_unf)
plotRegionFilterTotals(regionFilterTotals)

# Call Regions with Filtering ####
regions <- getRegions(bs)
plotRegionStats(regions)
plotSDstats(regions)

# Adjust Methylation Data for PCs ####
meth <- getRegionMeth(regions, bs = bs)
mod <- model.matrix(~1, data = pData(bs))
methAdj <- adjustRegionMeth(meth, mod = mod)
dendro <- getDendro(methAdj)
plotDendro(dendro, file = "Sample_Dendrogram.pdf")

# Select Soft Power Threshold ####
sft <- getSoftPower(methAdj)
plotSoftPower(sft)

# Get Comethylation Modules ####
modules <- getModules(methAdj, power = sft$powerEstimate)
plotRegionDendro(modules)
