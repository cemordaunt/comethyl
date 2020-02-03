# Comethyl Analysis Pipeline ------------------------------------------
# Charles Mordaunt

# Load Packages ####
.libPaths("/share/lasallelab/Charles/comethylated/R")
sapply(c("scales", "openxlsx", "rlist", "tidyverse", "ggdendro", "cowplot", "bsseq", "dmrseq", "WGCNA", "sva"), require, 
       character.only = TRUE)

# Functions ####
getCpGs <- function(colData, path = getwd(), pattern = "*CpG_report.txt.gz", 
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), BPPARAM = MulticoreParam(10), 
                    save = TRUE, file = "Unfiltered_BSseq.rds", verbose = TRUE){
        if(verbose){
                message("[getCpGs] Loading CpG-level data")
        }
        files <- list.files(path, pattern = pattern) %>% .[pmatch(rownames(colData), table = .)]
        loci <- bsseq:::.readBismarkAsFWGRanges(files[1]) %>% chrSelectBSseq(seqnames = chroms)
        bs <- read.bismark(files, loci = loci, colData = colData, BPPARAM = BPPARAM, verbose = verbose)
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

.nCpGsBySample <- function(covSample, nSample){
        n <- length(covSample[covSample >= nSample])
        return(n)
}

.nCpGsByCov <- function(bsCov, cov, nSample){
        covSample <- (bsCov >= cov) %>% DelayedMatrixStats::rowSums2()
        nCpGs <- sapply(nSample, FUN = .nCpGsBySample, covSample = covSample)
        return(nCpGs)
}

getCpGtotals <- function(bs, cov = seq(0,10,1), perSample = seq(0.5,1,0.05), save = TRUE, file = "CpG_Totals.txt",
                         verbose = TRUE){
        if(verbose){
                message("[getCpGtotals] Calculating CpG totals at specified cov and perSample cutoffs")
        }
        nSample <- (perSample * ncol(bs)) %>% ceiling()
        bsCov <- getCoverage(bs)
        nCpGs <- sapply(cov, FUN = .nCpGsByCov, bsCov = bsCov, nSample = nSample)
        nCpGs <- as.integer(nCpGs)
        perCpGs <- (nCpGs / length(bs)) %>% round(digits = 5)
        CpGtotals <- data.frame(cov = rep(cov, each = length(perSample)), perSample = rep(perSample, times = length(cov)), 
                                nSample = rep(nSample, times = length(cov)), nCpGs_M = nCpGs / 10^6, perCpGs = perCpGs)
        if(save){
                if(verbose){
                        message("[getCpGtotals] Saving file as ", file)
                }
                write.table(CpGtotals, file = file, sep = "\t", row.names = FALSE)
        }
        return(CpGtotals)
}

plotCpGtotals <- function(CpGtotals, nBreaks = 4, legend.position = c(1.08,0.73), save = TRUE, 
                          file = "CpG_Totals.pdf", width = 11, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotCpGtotals] Plotting CpG totals")
        }
        CpGtotals$perSample <- CpGtotals$perSample * 100
        gg <- ggplot(data = CpGtotals)
        gg <- gg +
                geom_line(aes(x = perSample, y = nCpGs_M, group = cov, color = cov)) +
                geom_text(data = subset(CpGtotals, perSample == min(perSample)), 
                          aes(x = perSample, y = nCpGs_M, group = cov, color = cov, label = cov), 
                          size = 4.5, check_overlap = TRUE, nudge_x = -0.5, hjust = 1) +
                xlab("Samples (%) Cutoff") +
                ylab("Total CpGs (Millions)") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks + 1), expand = expand_scale(mult = c(0.05, 0.03))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_color_gradient("Coverage\nCutoff", breaks = breaks_pretty(n = nBreaks - 1)) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 14, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"), axis.title = element_text(size = 18),
                      legend.background = element_blank(), legend.position = legend.position, 
                      legend.title = element_text(size = 18), legend.text = element_text(size = 14), 
                      panel.border = element_rect(color = "black", size = 1.25), 
                      panel.grid = element_blank(), plot.margin = unit(c(1,7,0.7,0.7), "lines"))
        if(save){
                if(verbose){
                        message("[plotCpGtotals] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

filterCpGs <- function(bs, cov = 2, perSample = 0.75, save = TRUE, file = "Filtered_BSseq.rds", verbose = TRUE){
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

getRegions <- function(bs, maxGap = 150, n = 3, save = TRUE, file = "Unfiltered_Regions.txt", verbose = TRUE){
        if(verbose){
                message("[getRegions] Calling regions when at least ", n, " CpGs are no more than ", maxGap, " bases apart")
        }
        regions <- bsseq:::regionFinder3(x = as.integer(rep(1, length(bs))), chr = as.character(seqnames(bs)), 
                                         positions = start(bs), maxGap = maxGap, verbose = FALSE)[["up"]]
        regions <- regions[regions$n >= n,]
        regions$RegionID <- paste("Region", 1:nrow(regions), sep = "_")
        regions$chr <- as.character(regions$chr)
        if(verbose){
                message("[getRegions] Calculating region statistics")
        }
        regions$width <- regions$end - regions$start
        cov <- getCoverage(bs, regions = regions[,c("chr", "start", "end")], what = "perRegionTotal")
        regions$covMin <- DelayedArray::rowMins(cov)
        regions$covMean <- DelayedMatrixStats::rowMeans2(cov) %>% round(digits = 5)
        regions$covSD <- DelayedMatrixStats::rowSds(cov) %>% round(digits = 5)
        meth <- getMeth(bs, regions = regions[,c("chr", "start", "end")], type = "raw", what = "perRegion")
        regions$methMean <- DelayedMatrixStats::rowMeans2(meth, na.rm = TRUE) %>% round(digits = 5)
        regions$methSD <- DelayedMatrixStats::rowSds(meth, na.rm = TRUE) %>% round(digits = 5)
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

plotRegionStats <- function(regions, maxQuantile = 1, bins = 75, histCol = "#132B43", lineCol = "red", nBreaks = 4, 
                            save = TRUE, file = "Region_Plots.pdf", width = 11, height = 8.5, verbose = TRUE){
        if(verbose){
                message("[plotRegionStats] Plotting histograms of region statistics")
        }
        variables <- c("width", "n", "covMin", "covMean", "methMean", "methSD")
        medians <- data.frame(variable = factor(variables, levels = variables), 
                              value = sapply(regions[,variables], FUN = median)) 
        if(maxQuantile < 1){
                if(verbose){
                        message("[plotRegionStats] Limiting x-axis to values in bottom ", maxQuantile * 100, 
                                "% for width, n, covMin, and covMean")
                }
                regions$width[regions$width >= quantile(regions$width, probs = maxQuantile)] <- NA
                regions$n[regions$n >= quantile(regions$n, probs = maxQuantile)] <- NA
                regions$covMin[regions$covMin >= quantile(regions$covMin, probs = maxQuantile)] <- NA
                regions$covMean[regions$covMean >= quantile(regions$covMean, probs = maxQuantile)] <- NA
        }
        regions <- reshape2::melt(regions[,c("RegionID", "width", "n", "covMin", "covMean", "methMean", "methSD")], 
                                  id.vars = "RegionID")
        gg <- ggplot(data = regions)
        gg <- gg +
                geom_histogram(aes(x = value), bins = bins, fill = histCol, color = histCol, na.rm = TRUE) +
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

plotSDstats <- function(regions, maxQuantile = 1, bins = 75, nBreaks = 4, legend.position = c(1.09,0.9), save = TRUE, 
                        file = "SD_Plots.pdf", width = 8.5, height = 8.5, verbose = TRUE){
        if(verbose){
                message("[plotSDstats] Plotting methylation SD vs region statistics")
        }
        if(maxQuantile < 1){
                if(verbose){
                        message("[plotSDstats] Limiting x-axis to values in bottom ", maxQuantile * 100, 
                                "% for n, covMin, and covMean")
                }
                regions$n[regions$n >= quantile(regions$n, probs = maxQuantile)] <- NA
                regions$covMin[regions$covMin >= quantile(regions$covMin, probs = maxQuantile)] <- NA
                regions$covMean[regions$covMean >= quantile(regions$covMean, probs = maxQuantile)] <- NA
        }
        regions <- reshape2::melt(regions[,c("RegionID", "n", "covMin", "covMean", "methMean", "methSD")], 
                                  id.vars = c("RegionID", "methSD"))
        gg <- ggplot(data = regions)
        gg <- gg +
                geom_bin2d(aes(x = value, y = methSD, color = ..count..), bins = bins, na.rm = TRUE) +
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

.regionTotals <- function(regions, covMin, methSD){
        regions <- regions[regions$covMin >= covMin & regions$methSD >= methSD,]
        totals <- c("covMin" = covMin, "methSD" = methSD, "totalRegions_K" = nrow(regions)/10^3, 
                    "totalWidth_Mb" = sum(regions$width)/10^6, "totalN_M" = sum(regions$n)/10^6)
        return(totals)
}

getRegionTotals <- function(regions, covMin = seq(0,20,2), methSD = seq(0,0.1,0.01),
                            save = TRUE, file = "Region_Totals.txt", verbose = TRUE){
        if(verbose){
                message("[getRegionTotals] Calculating region totals at specified covMin and methSD cutoffs")
        }
        covMin <- rep(covMin, each = length(methSD))
        methSD <- rep(methSD, times = length(covMin))
        regionTotals <- mapply(FUN = .regionTotals, covMin = covMin, methSD = methSD, 
                               MoreArgs = list(regions = regions)) %>% t() %>% as.data.frame()
        if(save){
                if(verbose){
                        message("[getRegionTotals] Saving file as ", file)
                        write.table(regionTotals, file = file, sep = "\t", row.names = FALSE)
                }
        }
        return(regionTotals)
}

plotRegionTotals <- function(regionTotals, nBreaks = 4, legend.position = c(1.08,0.897), save = TRUE, 
                             file = "Region_Totals.pdf", width = 11, height = 11, verbose = TRUE){
        if(verbose){
                message("[plotRegionTotals] Plotting region totals")
        }
        regionTotals <- reshape2::melt(regionTotals, id.vars = c("covMin", "methSD"))
        regionTotals$variable <- as.character(regionTotals$variable) %>% 
                str_replace_all(pattern = c("totalRegions_K" = "Total Regions (Thousands)", "totalWidth_Mb" = "Total Width (Mb)", 
                                            "totalN_M" = "Total CpGs (Millions)")) %>%
                factor(levels = c("Total Regions (Thousands)", "Total Width (Mb)", "Total CpGs (Millions)"))
        gg <- ggplot(data = regionTotals)
        gg <- gg +
                geom_line(aes(x = methSD, y = value, group = covMin, color = covMin)) +
                geom_text(data = subset(regionTotals, methSD == min(methSD)), 
                          aes(x = methSD, y = value, group = covMin, color = covMin, label = covMin), 
                          size = 4.5, check_overlap = TRUE, nudge_x = -0.001, hjust = 1) +
                facet_wrap(vars(variable), nrow = 3, ncol = 1, scales = "free_y", strip.position = "left") +
                xlab("SD Cutoff") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks), expand = expand_scale(mult = c(0.05, 0.03))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_color_gradient("Minimum\nCoverage\nCutoff", breaks = breaks_pretty(n = nBreaks - 1)) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 14, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"), axis.title.x = element_text(size = 18),
                      axis.title.y = element_blank(), legend.background = element_blank(), legend.position = legend.position, 
                      legend.title = element_text(size = 18), legend.text = element_text(size = 14), 
                      panel.border = element_rect(color = "black", size = 1.25), 
                      panel.grid = element_blank(), panel.spacing.x = unit(0.3, "lines"),
                      panel.spacing.y = unit(0.8, "lines"), plot.margin = unit(c(1,7,0.7,0.2), "lines"),
                      strip.background = element_blank(), strip.placement = "outside", 
                      strip.switch.pad.wrap = unit(0, "lines"), strip.text = element_text(size = 18))
        if(save){
                if(verbose){
                        message("[plotRegionTotals] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

filterRegions <- function(regions, covMin = 10, methSD = 0.05, save = TRUE, file = "Filtered_Regions.txt", verbose = TRUE){
        if(verbose){
                message("[getRegions] Filtering regions for at least ", covMin, 
                        " reads in all samples and methylation SD of at least ", methSD * 100, "%")
        }
        regions <- regions[regions$covMin >= covMin & regions$methSD >= methSD,]
        if(verbose){
                message("[getRegions] Creating new Region IDs")
        }
        regions$RegionID <- paste("Region", 1:nrow(regions), sep = "_")
        if(save){
                if(verbose){
                        message("[getRegions] Saving file as ", file)
                }
                write.table(regions, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
        }
        return(regions)
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
                        message("[getDendro] Clustering by euclidean distance")
                }
                dist <- dist(x)
        } else {
                if(distance == "pearson"){
                        if(verbose){
                                message("[getDendro] Clustering with pearson correlation as the distance")
                        }
                        dist <- (1 - WGCNA::cor(x, use = "pairwise.complete.obs")) %>% as.dist()
                } else {
                        if(distance == "bicor"){
                                if(verbose){
                                        message("[getDendro] Clustering with bicor correlation as the distance")
                                }
                                dist <- (1 - bicor(x, maxPOutliers = maxPOutliers, use = "pairwise.complete.obs")) %>% 
                                        as.dist()
                        } else {
                                stop("[getDendro] Error: Distance must be either euclidean, pearson, or bicor")
                        }
                }
        }
        dendro <- hclust(dist, method = "average")
        dendro$labels <- gsub("ME", replacement = "", x = dendro$labels, fixed = TRUE)
        return(dendro)
}

plotDendro <- function(dendro, label = TRUE, labelSize = 2.5, expandX = c(0.03,0.03), expandY = c(0.3,0.08), 
                       nBreaks = 4, save = TRUE, file = "Dendrogram.pdf", width = 11, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotDendro] Plotting dendrogram")
        }
        dendroPlot <- dendro_data(dendro)
        fix <- dendroPlot$segments$yend == 0
        dendroPlot$segments$yend[fix] <- dendroPlot$segments$y[fix] - max(dendroPlot$segments$y) * 0.05
        dendroPlot$labels$y <- dendroPlot$segments$yend[fix] - max(dendroPlot$segments$y) * 0.01
        dendroPlot$labels$label <- gsub("ME", replacement = "", x = dendroPlot$labels$label, fixed = TRUE)
        gg <- ggplot()
        gg <- gg +
                geom_segment(data = dendroPlot$segments, aes(x = x, y = y, xend = xend, yend = yend), lwd = 0.3, lineend = "square") +
                scale_x_continuous(expand = expand_scale(mult = expandX)) +
                scale_y_continuous(expand = expand_scale(mult = expandY), breaks = breaks_pretty(n = nBreaks)) +
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
                         RsquaredCut = 0.8, blockSize = 40000, verbose = TRUE){
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
                                         networkType = "signed", corFnc = "cor", blockSize = blockSize, verbose = verboseNum)
        } else {
                if(corType == "bicor"){
                        sft <- pickSoftThreshold(meth, RsquaredCut = RsquaredCut, powerVector = powerVector, 
                                                 networkType = "signed", corFnc = "bicor", 
                                                 corOptions = list(maxPOutliers = maxPOutliers), blockSize = blockSize, 
                                                 verbose = verboseNum)
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
                geom_text(aes(x = powerEstimate, y = min(0, min(fitIndices$value)), label = powerEstimate), color = lineCol, 
                          nudge_x = -1) +
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

getModules <- function(meth, maxBlockSize = 40000, power = NULL, corType = c("pearson", "bicor"), maxPOutliers = 0.1, deepSplit = 4, 
                       minModuleSize = 10, mergeCutHeight = 0.1, nThreads = 4, save = TRUE, file = "Modules.rds", verbose = TRUE){
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
        modules <- blockwiseModules(meth, checkMissingData = FALSE, maxBlockSize = maxBlockSize, corType = corType, 
                                    maxPOutliers = maxPOutliers, power = power, networkType = "signed", TOMtype = "signed", 
                                    deepSplit = deepSplit, minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, 
                                    nThreads = nThreads, verbose = verboseNum)
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

getModuleBED <- function(regions, modules, grey = TRUE, save = TRUE, file = "Modules.bed", verbose = TRUE){
        if(verbose){
                message("[getModuleBED] Creating bed file of regions annotated with identified modules")
        }
        regions$module <- modules$colors[match(regions$RegionID, names(modules$colors))]
        if(!grey){
                if(verbose){
                        message("[getModuleBED] Excluding regions in grey (unassigned) module")
                }
                regions <- regions[!regions$module == "grey",]
        }
        regions$rgb <- col2rgb(regions$module) %>% apply(2, paste, collapse = ",")
        bed <- cbind(regions[c("chr", "start", "end", "RegionID")], score = 0, strand = ".", thickStart = 0, thickEnd = 0, 
                     rgb = regions$rgb)
        if(save){
                if(verbose){
                        message("[getModuleBED] Saving file as ", file)
                }
                name <- gsub(".bed", replacement = "", file)
                write(paste("track name='", name, "' description='", name, "' itemRgb='On'", sep = ""), file = file)
                write.table(bed, file = file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        }
        return(bed)
}

plotHeatmap <- function(x, rowDendro, colDendro, colors = blueWhiteRed(100, gamma = 0.3), limit = max(abs(x)), 
                        axis.text.size = 8, legend.title = "Bicor", legend.title.size = 16, legend.text.size = 12,
                        legend.position = c(0.3,0.905), rowDendroMargins = c(-1.55,1,-0.1,-1.1), 
                        colDendroMargins = c(1,-0.5,-1,0.8), rowColorMargins = c(-1.85,-1.5,0.55,1.7),
                        colColorMargins = c(-1.6,-0.85,1.8,0.55), save = TRUE, file = "Heatmap.pdf", width = 11, 
                        height = 9.5, verbose = TRUE){
        if(verbose){
                message("[plotHeatmap] Plotting heatmap with dendrograms")
        }
        limits <- c(-limit, limit)
        x <- as.data.frame(x)
        rownames(x) <- gsub("ME", replacement = "", x = rownames(x), fixed = TRUE)
        colnames(x) <- gsub("ME", replacement = "", x = colnames(x), fixed = TRUE)
        rowModules <- sum(!is.na(col2hcl(rownames(x)))) == length(rownames(x))
        colModules <- sum(!is.na(col2hcl(colnames(x)))) == length(colnames(x))
        x$rowID <- factor(rownames(x), levels = rowDendro$labels[rev(rowDendro$order)])
        x <- reshape2::melt(x, id.vars = "rowID")
        x$variable <- factor(x$variable, levels = colDendro$labels[colDendro$order])
        hmMarginL <- ifelse(rowModules, yes = 2, no = -1)
        hmMarginB <- ifelse(colModules, yes = 2, no = -1)
        heatmap <- ggplot(data = x) +
                geom_tile(aes(x = variable, y = rowID, color = value, fill = value)) +
                scale_fill_gradientn(legend.title, colors = colors, limits = limits, aesthetics = c("color", "fill")) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_text(size = axis.text.size, color = "black", angle = 90, vjust = 0.5), 
                      axis.text.y = element_text(size = axis.text.size, color = "black"), 
                      axis.ticks = element_line(size = 0.8, color = "black"), 
                      axis.title = element_blank(), legend.position = "none", panel.background = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25), panel.grid = element_blank(), 
                      plot.background = element_blank(), plot.margin = unit(c(0,1,hmMarginB,hmMarginL), "lines"))
        legend <- get_legend(heatmap + theme(legend.position = legend.position, legend.background = element_blank(),
                                             legend.title = element_text(size = legend.title.size), 
                                             legend.text = element_text(size = legend.text.size)))
        rowDendroPlot <- ggplot(data = dendro_data(rowDendro)$segments) +
                geom_segment(aes(x = -x, y = y, xend = -xend, yend = yend), lwd = 0.5, lineend = "square") +
                coord_flip() +
                theme_dendro() +
                theme(plot.margin = unit(rowDendroMargins, "lines"))
        colDendroPlot <- ggplot(data = dendro_data(colDendro)$segments) +
                geom_segment(aes(x = x, y = y, xend = xend, yend = yend), lwd = 0.5, lineend = "square") +
                theme_dendro() +
                theme(plot.margin = unit(colDendroMargins, "lines"))
        rowColors <- NULL
        colColors <- NULL
        if(rowModules){
                if(verbose){
                        message("[plotHeatmap] Using colors in row names for y-axis labels")
                }
                rowColors <- ggplot(data = data.frame(x = 0, y = 1:length(levels(x$rowID)), color = levels(x$rowID))) +
                        geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                        scale_fill_identity(aesthetics = c("color", "fill")) +
                        theme_void() +
                        theme(legend.position = "none", plot.margin = unit(rowColorMargins, "lines"))
                heatmap <- heatmap + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        }
        if(colModules){
                if(verbose){
                        message("[plotHeatmap] Using colors in column names for x-axis labels")
                }
                colColors <- ggplot(data = data.frame(x = 1:length(levels(x$variable)), y = 0, color = levels(x$variable))) +
                        geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                        scale_fill_identity(aesthetics = c("color", "fill")) +
                        theme_void() +
                        theme(legend.position = "none", plot.margin = unit(colColorMargins, "lines"))
                heatmap <- heatmap + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        }
        gg <- plot_grid(NULL, colDendroPlot, NULL, NULL, rowColors, heatmap, rowDendroPlot, legend, NULL, colColors, 
                        NULL, NULL, nrow = 3, ncol = 4, rel_widths = c(0.045, 1, 0.15, 0.15), 
                        rel_heights = c(0.15, 1, 0.045))
        if(save){
                if(verbose){
                        message("[plotHeatmap] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

getCor <- function(x, y = NULL, transpose = FALSE, corType = c("bicor", "pearson"), maxPOutliers = 0.1, verbose = TRUE){
        if(transpose){
                if(verbose){
                        message("[getCor] Transposing data")
                }
                x <- t(x)
        }
        corType <- match.arg(corType)
        if(verbose){
                message("[getCor] Calculating correlations using ", corType, " correlation")
        }
        if(corType == "bicor"){
                cor <- bicor(x, y = y, use = "pairwise.complete.obs", maxPOutliers = maxPOutliers, 
                             pearsonFallback = "none")
        } else {
                if(corType == "pearson"){
                        cor <- WGCNA::cor(x, y = y, use = "pairwise.complete.obs")
                } else {
                        stop("[getCor] corType must be either bicor or pearson")
                }
        }
        return(cor)
}

getMEtraitCor <- function(MEs, colData, corType = c("bicor", "pearson"), maxPOutliers = 0.1, robustY = FALSE, 
                          save = TRUE, file = "ME_Trait_Correlation_Stats.txt", verbose = TRUE){
        corType <- match.arg(corType)
        if(verbose){
                message("[getMEtraitCor] Testing associations between module eigennodes and sample traits using ",
                        corType, " correlation")
        }
        colData <- colData[rownames(MEs),]
        if(corType == "bicor"){
                cor <- bicorAndPvalue(x = MEs, y = colData, maxPOutliers = maxPOutliers, robustY = robustY, 
                                      pearsonFallback = "none")
        } else {
                if(corType == "pearson"){
                        cor <- corAndPvalue(x = MEs, y = colData)
                } else {
                        stop("[getMEtraitCor] corType must be either bicor or pearson")
                }
        }
        stats <- list.rbind(cor) %>% as.data.frame()
        stats$module <- rownames(cor$p) %>% gsub(pattern = "ME", replacement = "", x = .) %>% rep(length(cor)) %>%
                factor(levels = unique(.))
        stats$stat <- names(cor) %>% rep(each = nrow(cor$p)) %>% factor(levels = unique(.))
        stats <- reshape2::melt(stats, id.vars = c("module", "stat")) %>%
                reshape2::dcast(formula = module + variable ~ stat, value.var = "value")
        stats$q <- p.adjust(stats$p, method = "fdr")
        if(corType == "bicor"){
                stats <- stats[,c("module", "variable", "nObs", "bicor", "Z", "t", "p", "q")]
        } else {
                stats <- stats[,c("module", "variable", "nObs", "cor", "Z", "t", "p", "q")]
        }
        if(save){
                if(verbose){
                        message("[getMEtraitCor] Saving file as ", file)
                }
                write.table(stats, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
        }
        return(stats)
}

plotMEtraitCor <- function(MEtraitCor, sigOnly = FALSE, star.size = 8, star.nudge_y = -0.38,
                           colors = blueWhiteRed(100, gamma = 0.9), limit = max(abs(MEtraitCor$bicor)), 
                           axis.text.size = 12, legend.position = c(1.08, 0.915), legend.text.size = 12, 
                           legend.title.size = 16, colColorMargins = c(-0.7,4.21,1.2,11.07), save = TRUE, 
                           file = "ME_Trait_Correlation_Heatmap.pdf", width = 11, height = 9.5, verbose = TRUE){
        if(verbose){
                message("[plotMEtraitCor] Plotting ME trait correlation heatmap")
        }
        MEtraitCor$Significant <- (MEtraitCor$q < 0.05) %>% factor(levels = c("TRUE", "FALSE"))
        if(sigOnly){
                sigVars <- MEtraitCor$variable[MEtraitCor$Significant == "TRUE"] %>% unique() %>% as.character()
                sigMods <- MEtraitCor$module[MEtraitCor$Significant == "TRUE"] %>% unique() %>% as.character()
                MEtraitCor <- subset(MEtraitCor, variable %in% sigVars & module %in% sigMods)
                MEtraitCor$module <- factor(MEtraitCor$module, 
                                            levels = levels(MEtraitCor$module)[levels(MEtraitCor$module) %in% sigMods])
        }
        heatmap <- ggplot(data = MEtraitCor) +
                geom_tile(aes(x = module, y = variable, color = bicor, fill = bicor)) +
                geom_text(aes(x = module, y = variable, alpha = Significant), label = "*", color = "black", 
                          size = star.size, nudge_y = star.nudge_y) +
                scale_fill_gradientn("Bicor", colors = colors, limits = c(-limit, limit), 
                                     aesthetics = c("color", "fill")) +
                scale_x_discrete(expand = expand_scale(mult = 0.01)) +
                scale_y_discrete(expand = expand_scale(mult = 0.01)) +
                scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_blank(), 
                      axis.text.y = element_text(size = axis.text.size, color = "black"), 
                      axis.ticks.x = element_blank(), axis.ticks.y = element_line(size = 0.8, color = "black"), 
                      axis.title = element_blank(), legend.background = element_blank(), 
                      legend.position = legend.position, 
                      legend.text = element_text(size = legend.text.size), 
                      legend.title = element_text(size = legend.title.size), 
                      panel.background = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
                      panel.grid = element_blank(), plot.background = element_blank(), 
                      plot.margin = unit(c(1,6,1,1), "lines"))
        colColors <- ggplot(data = data.frame(x = 1:length(levels(MEtraitCor$module)), y = 0, 
                                              color = levels(MEtraitCor$module))) +
                geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                scale_fill_identity(aesthetics = c("color", "fill")) +
                theme_void() +
                theme(legend.position = "none", plot.margin = unit(colColorMargins, "lines"))
        gg <- plot_grid(heatmap, colColors, nrow = 2, rel_heights = c(1, 0.045))
        if(save){
                if(verbose){
                        message("[plotMEtraitCor] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

# Set Global Options ####
options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
enableWGCNAThreads(nThreads = 4)

# Read Bismark CpG Reports ####
colData <- read.xlsx("sample_info.xlsx", rowNames = TRUE)
bs <- getCpGs(colData)

# Examine CpG Totals at Different Cutoffs ####
CpGtotals <- getCpGtotals(bs)
plotCpGtotals(CpGtotals)

# Filter BSobject ####
bs <- filterCpGs(bs, cov = 2, perSample = 0.75)

# Call Regions ####
regions <- getRegions(bs)
plotRegionStats(regions, maxQuantile = 0.99, file = "Unfiltered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99, file = "Unfiltered_SD_Plots.pdf")

# Examine Region Totals at Different Cutoffs ####
regionTotals <- getRegionTotals(regions)
plotRegionTotals(regionTotals)

# Filter Regions ####
regions <- filterRegions(regions, covMin = 10, methSD = 0.05)
plotRegionStats(regions, maxQuantile = 0.99, file = "Filtered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99, file = "Filtered_SD_Plots.pdf")

# Adjust Methylation Data for PCs ####
meth <- getRegionMeth(regions, bs = bs)
mod <- model.matrix(~1, data = pData(bs))
methAdj <- adjustRegionMeth(meth, mod = mod)
getDendro(methAdj, distance = "euclidean") %>% plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))

# Select Soft Power Threshold ####
sft <- getSoftPower(methAdj, corType = "pearson")
plotSoftPower(sft)

# Get Comethylation Modules ####
modules <- getModules(methAdj, power = sft$powerEstimate, corType = "pearson")
plotRegionDendro(modules)
BED <- getModuleBED(regions, modules = modules)

# Examine Correlations between Modules and Samples ####
MEs <- modules$MEs
moduleDendro <- getDendro(MEs, distance = "bicor")
plotDendro(moduleDendro, file = "Module_ME_Dendrogram.pdf", labelSize = 4, nBreaks = 5)
moduleCor <- getCor(MEs, corType = "bicor")
plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro, file = "Module_Correlation_Heatmap.pdf")

sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
plotDendro(sampleDendro, file = "Sample_ME_Dendrogram.pdf", labelSize = 3, nBreaks = 5)
sampleCor <- getCor(MEs, transpose = TRUE, corType = "bicor")
plotHeatmap(sampleCor, rowDendro = sampleDendro, colDendro = sampleDendro, file = "Sample_Correlation_Heatmap.pdf")

plotHeatmap(MEs, rowDendro = sampleDendro, colDendro = moduleDendro, file = "Sample_ME_Heatmap.pdf",
            legend.title = "Module\nEigennode", legend.position = c(0.37,0.89))

# Test Correlations between Module Eigennodes and Sample Traits ####
MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor", file = "ME_Trait_Correlation_Stats.txt")
traitDendro <- reshape2::acast(data = MEtraitCor, formula = variable ~ module, value.var = "bicor") %>%
        dist() %>% hclust(method = "average")
MEtraitCor$module <- factor(MEtraitCor$module, levels = moduleDendro$labels[moduleDendro$order])
MEtraitCor$variable <- factor(MEtraitCor$variable, levels = rev(traitDendro$labels[traitDendro$order]))
plotMEtraitCor(MEtraitCor, file = "ME_Trait_Correlation_Heatmap.pdf")
plotMEtraitCor(MEtraitCor, sigOnly = TRUE, star.size = 11, star.nudge_y = -0.27, legend.position = c(1.14, 0.745),
               colColorMargins = c(-1,5.1,0.5,10.47), file = "Sig_ME_Trait_Correlation_Heatmap.pdf", 
               width = 7, height = 3.5)
