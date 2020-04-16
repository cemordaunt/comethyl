# Comethyl Analysis Pipeline ------------------------------------------
# Charles Mordaunt

# Load Packages ####
.libPaths("/share/lasallelab/programs/comethyl/R_3.6")
AnnotationHub::setAnnotationHubOption("CACHE", value = "/share/lasallelab/programs/comethyl/R_3.6")
sapply(c("scales", "openxlsx", "rlist", "tidyverse", "ggdendro", "cowplot", "annotatr", "rtracklayer", "bsseq", "dmrseq", 
         "WGCNA", "sva", "rGREAT", "R.devices", "biomaRt"), require, character.only = TRUE)

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
        perCpGs <- (nCpGs / length(bs))
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

getRegions <- function(bs, annotation = NULL, genome = c("hg38", "hg19", "mm10", "mm9", "rn6", "rn5", "rn4", "dm6", "dm3", "galGal5"), 
                       upstream = 5000, downstream = 1000, custom = NULL, maxGap = 150, n = 3, save = TRUE, 
                       file = "Unfiltered_Regions.txt", verbose = TRUE){
        if(!is.null(annotation) & !is.null(custom)){
                stop("[getRegions] annotation and custom cannot both have values")
        }
        if(!is.null(annotation)){
                genome <- match.arg(genome)
                if(verbose){
                        message("[getRegions] Using ", annotation, " annotation for the ", genome, 
                                " genome as regions")
                }
                if(annotation %in% c("genes", "promoters", "transcripts")){
                        txdb <- annotatr:::get_txdb_name(genome)
                        if(requireNamespace(txdb, quietly = TRUE)){
                                library(txdb, character.only = TRUE)
                        }
                        txdb <- get(txdb)
                        orgdb <- annotatr:::get_orgdb_name(genome)
                        if (requireNamespace(sprintf("org.%s.eg.db", orgdb), quietly = TRUE)) {
                                library(sprintf("org.%s.eg.db", orgdb), character.only = TRUE)
                        }
                        regions <- suppressWarnings(switch(annotation, 
                                                           genes = genes(txdb),
                                                           promoters = promoters(txdb, upstream = upstream, 
                                                                                 downstream = downstream, 
                                                                                 use.names = FALSE),
                                                           transcripts = transcripts(txdb)))
                } else {
                        if(annotation %in% builtin_annotations()){
                                if(!grepl(genome, x = annotation, fixed = TRUE)){
                                        stop("[getRegions] Annotation must match genome")
                                }
                                regions <- build_annotations(genome, annotations = annotation)
                        } else {
                                stop("[getRegions] Annotation not supported")
                        }
                }
                regions <- keepStandardChromosomes(regions, pruning.mode = "coarse") %>% trim()
                regions$n <- countOverlaps(regions, subject = bs)
                regions <- as.data.frame(regions)
                colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
        } else {
                if(!is.null(custom)){
                        if(verbose){
                                message("[getRegions] Using custom regions")
                        }
                        regions <- custom
                        regions$n <- countOverlaps(regions, subject = bs)
                        regions <- as.data.frame(regions)
                        colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
                } else {
                        if(verbose){
                                message("[getRegions] Calling regions when at least ", n, " CpGs are no more than ", 
                                        maxGap, " bases apart")
                        }
                        regions <- bsseq:::regionFinder3(x = as.integer(rep(1, length(bs))), chr = as.character(seqnames(bs)), 
                                                         positions = start(bs), maxGap = maxGap, verbose = FALSE)[["up"]]
                }
        }
        regions <- regions[regions$n >= n,] %>% .[with(., order(chr, start, end)),]
        regions$RegionID <- paste("Region", 1:nrow(regions), sep = "_")
        regions$chr <- as.character(regions$chr)
        if(verbose){
                message("[getRegions] Calculating region statistics")
        }
        regions$width <- regions$end - regions$start
        cov <- getCoverage(bs, regions = regions[,c("chr", "start", "end")], what = "perRegionTotal")
        regions$covMin <- DelayedArray::rowMins(cov)
        regions$covMean <- DelayedMatrixStats::rowMeans2(cov)
        regions$covSD <- DelayedMatrixStats::rowSds(cov)
        meth <- getMeth(bs, regions = regions[,c("chr", "start", "end")], type = "raw", what = "perRegion")
        regions$methMean <- DelayedMatrixStats::rowMeans2(meth, na.rm = TRUE)
        regions$methSD <- DelayedMatrixStats::rowSds(meth, na.rm = TRUE)
        colnames <- c("RegionID", "chr", "start", "end", "width", "n", "covMin", "covMean", "covSD", "methMean", "methSD")
        regions <- regions[,c(colnames, colnames(regions)[!colnames(regions) %in% colnames & 
                                                                  !colnames(regions) %in% c("idxStart", "idxEnd", "cluster")])]
        if(save){
                if(verbose){
                        message("[getRegions] Saving file as ", file)
                }
                write.table(regions, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
        }
        return(regions)
}

plotRegionStats <- function(regions, maxQuantile = 1, bins = 30, histCol = "#132B43", lineCol = "red", nBreaks = 4, 
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

plotSDstats <- function(regions, maxQuantile = 1, bins = 30, nBreaks = 4, legend.position = c(1.09,0.9), save = TRUE, 
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
        length_covMin <- length(covMin)
        length_methSD <- length(methSD)
        covMin <- rep(covMin, each = length_methSD)
        methSD <- rep(methSD, times = length_covMin)
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
        dendro$labels <- str_remove_all(dendro$labels, pattern = "ME")
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
        dendroPlot$labels$label <- str_remove_all(dendroPlot$labels$label, pattern = "ME")
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
                         RsquaredCut = 0.8, blockSize = 40000, gcInterval = blockSize - 1, save = TRUE,
                         file = "Soft_Power.rds", verbose = TRUE){
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
        if(save){
                if(verbose){
                        message("[getSoftPower] Saving file as ", file)
                }
                saveRDS(sft, file = file)
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

getModules <- function(meth, power, regions, maxBlockSize = 40000, corType = c("pearson", "bicor"), 
                       maxPOutliers = 0.1, deepSplit = 4, minModuleSize = 10, mergeCutHeight = 0.1, nThreads = 4, save = TRUE, 
                       file = "Modules.rds", verbose = TRUE){
        if(is.null(power)){
                stop("[getModules] You must select a soft power threshold")
        }
        if(is.null(regions)){
                stop("[getModules] You must include regions")
        }
        corType <- match.arg(corType)
        if(!corType %in% c("pearson", "bicor")){
                stop("[getModules] corType must be either pearson or bicor")
        }
        if(verbose){
                message("[getModules] Constructing network and detecting modules in blocks using ", corType, " correlation")
                verboseNum <- 10
        } else {
                verboseNum <- 0
        }
        modules <- blockwiseModules(meth, checkMissingData = FALSE, maxBlockSize = maxBlockSize, corType = corType, 
                                    maxPOutliers = maxPOutliers, power = power, networkType = "signed", TOMtype = "signed", 
                                    deepSplit = deepSplit, minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, 
                                    nThreads = nThreads, verbose = verboseNum)
        if(verbose){
                message("[getModules] Assigning modules and calculating module membership using ", corType, " correlation")
        }
        colnames(modules$MEs) <- str_remove_all(colnames(modules$MEs), pattern = "ME")
        if(corType == "pearson"){
                membership <- WGCNA::cor(x = meth, y = modules$MEs, use = "pairwise.complete.obs", nThreads = nThreads)
        } else {
                membership <- bicor(x = meth, y = modules$MEs, use = "pairwise.complete.obs", maxPOutliers = maxPOutliers,
                                    nThreads = nThreads)
        }
        regions$module <- modules$colors[match(regions$RegionID, names(modules$colors))]
        regions <- lapply(unique(regions$module), function(x){
                regions <- regions[regions$module == x,]
                regions$membership <- membership[regions$RegionID,x]
                regions$hubRegion <- regions$membership == max(regions$membership)
                return(regions)
        })
        regions <- list.rbind(regions) %>% .[order(as.integer(str_remove_all(.$RegionID, pattern = "Region_"))),]
        modules$regions <- regions
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

getModuleBED <- function(regions, grey = FALSE, save = TRUE, file = "Modules.bed", verbose = TRUE){
        if(!"module" %in% colnames(regions)){
                stop("[getModuleBED] Regions must have module annotation")
        }
        if(verbose){
                message("[getModuleBED] Creating bed file of regions annotated with identified modules")
        }
        if(!grey){
                if(verbose){
                        message("[getModuleBED] Excluding regions in grey (unassigned) module")
                }
                regions <- regions[!regions$module == "grey",]
        }
        regions$RegionID <- paste(regions$RegionID, regions$module, sep = "_")
        regions$rgb <- col2rgb(regions$module) %>% apply(2, paste, collapse = ",")
        BED <- cbind(regions[c("chr", "start", "end", "RegionID")], score = 0, strand = ".", thickStart = 0, thickEnd = 0, 
                     rgb = regions$rgb)
        if(save){
                if(verbose){
                        message("[getModuleBED] Saving file as ", file)
                }
                name <- basename(file) %>% str_remove_all(pattern = ".bed")
                write(paste("track name='", name, "' description='", name, "' itemRgb='On'", sep = ""), file = file)
                write.table(BED, file = file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        }
        return(BED)
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
        rownames(x) <- str_remove_all(rownames(x), pattern = "ME")
        colnames(x) <- str_remove_all(colnames(x), pattern = "ME")
        rowModules <- sum(rownames(x) %in% colors()) == length(rownames(x))
        colModules <- sum(colnames(x) %in% colors()) == length(colnames(x))
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

getCor <- function(x, y = NULL, transpose = FALSE, corType = c("bicor", "pearson"), maxPOutliers = 0.1, robustY = TRUE,
                   verbose = TRUE){
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
                             robustY = robustY, pearsonFallback = "none")
        } else {
                if(corType == "pearson"){
                        cor <- WGCNA::cor(x, y = y, use = "pairwise.complete.obs")
                } else {
                        stop("[getCor] corType must be either bicor or pearson")
                }
        }
        return(cor)
}

getMEtraitCor <- function(MEs, colData, corType = c("bicor", "pearson"), maxPOutliers = 0.1, robustY = FALSE, adjMethod = "fdr", 
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
        stats$module <- rownames(cor$p) %>% str_remove_all(pattern = "ME") %>% rep(length(cor)) %>%
                factor(levels = unique(.))
        stats$stat <- names(cor) %>% rep(each = nrow(cor$p)) %>% factor(levels = unique(.))
        stats <- reshape2::melt(stats, id.vars = c("module", "stat"), variable.name = "trait") %>%
                reshape2::dcast(formula = module + trait ~ stat, value.var = "value")
        stats$adj_p <- p.adjust(stats$p, method = adjMethod)
        if(corType == "bicor"){
                stats <- stats[,c("module", "trait", "nObs", "bicor", "Z", "t", "p", "adj_p")]
        } else {
                stats <- stats[,c("module", "trait", "nObs", "cor", "Z", "t", "p", "adj_p")]
        }
        if(save){
                if(verbose){
                        message("[getMEtraitCor] Saving file as ", file)
                }
                write.table(stats, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
        }
        return(stats)
}

plotMEtraitCor <- function(MEtraitCor, moduleOrder = 1:length(unique(MEtraitCor$module)), 
                           traitOrder = 1:length(unique(MEtraitCor$trait)), sigOnly = FALSE, adj_p = 0.05, star.size = 8, 
                           star.nudge_y = -0.38, colors = blueWhiteRed(100, gamma = 0.9), limit = NULL, 
                           axis.text.size = 12, legend.position = c(1.08, 0.915), legend.text.size = 12, 
                           legend.title.size = 16, colColorMargins = c(-0.7,4.21,1.2,11.07), save = TRUE, 
                           file = "ME_Trait_Correlation_Heatmap.pdf", width = 11, height = 9.5, verbose = TRUE){
        if(verbose){
                message("[plotMEtraitCor] Plotting ME trait correlation heatmap")
        }
        MEtraitCor$module <- factor(MEtraitCor$module, levels = levels(MEtraitCor$module)[moduleOrder])
        MEtraitCor$trait <- factor(MEtraitCor$trait, levels = levels(MEtraitCor$trait)[rev(traitOrder)])
        MEtraitCor$Significant <- (MEtraitCor$adj_p < adj_p & !is.na(MEtraitCor$adj_p)) %>% factor(levels = c("TRUE", "FALSE"))
        if(sigOnly){
                sigModules <- MEtraitCor$module[MEtraitCor$Significant == "TRUE"] %>% unique() %>% as.character()
                sigTraits <- MEtraitCor$trait[MEtraitCor$Significant == "TRUE"] %>% unique() %>% as.character()
                MEtraitCor <- subset(MEtraitCor, module %in% sigModules & trait %in% sigTraits)
                MEtraitCor$module <- factor(MEtraitCor$module, 
                                            levels = levels(MEtraitCor$module)[levels(MEtraitCor$module) %in% sigModules])
        }
        if("bicor" %in% colnames(MEtraitCor)){
                corType <- "bicor"
        } else {
                if("cor" %in% colnames(MEtraitCor)){
                        corType <- "cor"
                } else {
                        stop("[plotMEtraitCor] corType unknown, must be either bicor or cor")
                }
        }
        if(is.null(limit)){
                limit <- max(abs(MEtraitCor[[corType]]))
        }
        heatmap <- ggplot(data = MEtraitCor) +
                geom_tile(aes(x = module, y = trait, color = MEtraitCor[[corType]], fill = MEtraitCor[[corType]])) +
                geom_text(aes(x = module, y = trait, alpha = Significant), label = "*", color = "black", 
                          size = star.size, nudge_y = star.nudge_y) +
                scale_fill_gradientn(str_to_title(corType), colors = colors, limits = c(-limit, limit), 
                                     aesthetics = c("color", "fill")) +
                scale_x_discrete(expand = expand_scale(mult = 0.01)) +
                scale_y_discrete(expand = expand_scale(mult = 0.01)) +
                scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c("TRUE" = 1, "FALSE" = 0), guide = FALSE) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_blank(), axis.text.y = element_text(size = axis.text.size, color = "black"), 
                      axis.ticks.x = element_blank(), axis.ticks.y = element_line(size = 0.8, color = "black"), 
                      axis.title = element_blank(), legend.background = element_blank(), legend.position = legend.position, 
                      legend.text = element_text(size = legend.text.size), legend.title = element_text(size = legend.title.size), 
                      panel.background = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
                      panel.grid = element_blank(), plot.background = element_blank(), plot.margin = unit(c(1,6,1,1), "lines"))
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

plotMEtraitDot <- function(ME, trait, traitCode = NULL, colors = NULL, fun.data = "median_hilow", 
                           fun.args = list(conf.int = 0.5), binwidth = 0.01, stackratio = 1.4, dotsize = 0.85, ylim = NULL,
                           nBreaks = 4, axis.title.size = 20, axis.text.size = 16, xlab = "Trait", ylab = "Module Eigennode",
                           save = TRUE, file = "ME_Trait_Dotplot.pdf", width = 6, height = 6, verbose = TRUE){
        if(verbose){
                message("[plotMEtraitDot] Plotting module eigennode by categorical trait")
        }
        if(!is.null(traitCode)){
                trait <- names(traitCode)[match(trait, traitCode)] %>% factor(levels = names(traitCode))
        } else {
                trait <- as.factor(trait)
        }
        dotplot <- ggplot() +
                stat_summary(aes(x = trait, y = ME, group = trait), fun.data = fun.data, geom = "crossbar", 
                             color = "black", size = 0.5, fun.args = fun.args) +
                geom_dotplot(aes(x = trait, y = ME, fill = trait, color = trait), binwidth = binwidth, binaxis = "y", 
                             stackdir = "center", position = "dodge", stackratio = stackratio, dotsize = dotsize) +
                coord_cartesian(ylim = ylim) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                theme_bw(base_size = 25) +
                theme(legend.position = "none", panel.grid.major = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
                      panel.grid.minor = element_blank(), strip.background = element_blank(),
                      axis.text = element_text(color = "black", size = axis.text.size), 
                      axis.title = element_text(size = axis.title.size), plot.margin = unit(c(1,1,1,1), "lines")) +
                xlab(xlab) +
                ylab(ylab)
        if(!is.null(colors)){
                dotplot <- dotplot +
                        scale_color_manual(breaks = names(colors), values = colors, aesthetics = c("color", "fill"))
        }
        if(verbose){
                message("[plotMEtraitDot] Saving file as ", file)
        }
        ggsave(file, plot = dotplot, dpi = 600, width = width, height = height, units = "in")
}

plotMEtraitScatter <- function(ME, trait, color = "#132B43", xlim = NULL, ylim = NULL, nBreaks = 4, point.size = 2.5, 
                               axis.title.size = 20, axis.text.size = 16, xlab = "Trait", ylab = "Module Eigennode", 
                               save = TRUE, file = "ME_Trait_Scatterplot.pdf", width = 6, height = 6, verbose = TRUE){
        if(verbose){
                message("[plotMEtraitScatter] Plotting module eigennode by continuous trait")
        }
        scatterplot <- ggplot() +
                geom_smooth(aes(x = trait, y = ME), method = MASS::rlm, color = "#56B1F7", fill = "#336A98") +   
                geom_point(aes(x = trait, y = ME), color = color, size = point.size) +   
                coord_cartesian(xlim = xlim, ylim = ylim) +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                theme_bw(base_size = 25) +
                theme(legend.position = "none", panel.grid.major = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
                      panel.grid.minor = element_blank(), strip.background = element_blank(),
                      axis.text = element_text(color = "black", size = axis.text.size), 
                      axis.title = element_text(size = axis.title.size), plot.margin = unit(c(1,1,1,1), "lines")) +
                xlab(xlab) +
                ylab(ylab)
        if(verbose){
                message("[plotMEtraitScatter] Saving file as ", file)
        }
        ggsave(file, plot = scatterplot, dpi = 600, width = width, height = height, units = "in")
}

plotMethTrait <- function(module, regions, meth, trait, discrete = NULL, traitCode = NULL, traitColors = NULL, 
                          heatmapColors = blueWhiteRed(100, gamma = 0.3), limit = NULL, expandY = 0.05, axis.text.size = 11,
                          heatmap.legend.position = c(1.1,0.743), trait.legend.position = c(1.017,4.39),
                          heatmap.legend.title = "Relative\nMethylation (%)", trait.legend.title = "Trait",
                          legend.text.size = 11, legend.title.size = 14, heatmapMargins = c(1,8,0,1),
                          traitMargins = c(0,6,1,5.15), save = TRUE, file = "Module_Methylation_Trait_Heatmap.pdf",
                          width = 11, height = 4, verbose = TRUE){
        if(verbose){
                message("[plotMethTrait] Plotting ", module, " module region methylation by ", trait.legend.title)
        }
        RegionIDs <- rev(regions$RegionID[regions$module == module])
        if(is.null(discrete)){
                discrete <- length(unique(trait)) <= 5
        }
        if(discrete){
                if(!is.null(traitCode)){
                        trait <- names(traitCode)[match(trait, traitCode)] %>% factor(levels = names(traitCode))
                } else {
                        trait <- as.factor(trait)
                }
        }
        meth <- meth[RegionIDs,order(trait)]
        meth <- (meth - DelayedMatrixStats::rowMeans2(meth)) %>% reshape2::melt()
        if(is.null(limit)){
                limit <- max(abs(meth$value))
        }
        heatmap <- ggplot(data = meth) +
                geom_tile(aes(x = Var2, y = Var1, color = value, fill = value)) +
                scale_fill_gradientn(heatmap.legend.title, colors = heatmapColors, limits = c(-limit,limit), 
                                     aesthetics = c("color", "fill")) +
                scale_y_discrete(expand = expand_scale(mult = expandY)) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_blank(), axis.text.y = element_text(size = axis.text.size, color = "black"),
                      axis.ticks.x = element_blank(), axis.ticks.y = element_line(size = 0.8, color = "black"),
                      axis.title = element_blank(), legend.background = element_blank(), 
                      legend.position = heatmap.legend.position, legend.text = element_text(size = legend.text.size), 
                      legend.title = element_text(size = legend.title.size), panel.background = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), panel.grid = element_blank(), 
                      plot.background = element_blank(), plot.margin = unit(heatmapMargins, "lines"))
        colColors <- ggplot(data = data.frame(x = 1:length(trait), y = 0, color = sort(trait))) +
                geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                theme_void() +
                theme(legend.position = trait.legend.position, legend.text = element_text(size = legend.text.size), 
                      legend.title = element_text(size = legend.title.size), plot.margin = unit(traitMargins, "lines"))
        if(discrete){
                if(!is.null(traitColors)){
                        colColors <- colColors +
                                scale_color_manual(trait.legend.title, breaks = names(traitColors), values = traitColors, 
                                                   aesthetics = c("color", "fill"))
                } else {
                        colColors <- colColors +
                                scale_color_discrete(trait.legend.title, aesthetics = c("color", "fill"))
                }
        } else {
                colColors <- colColors +
                        scale_color_continuous(trait.legend.title, aesthetics = c("color", "fill"))
        }
        gg <- plot_grid(heatmap, colColors, nrow = 2, rel_heights = c(1,0.15))
        if(save){
                if(verbose){
                        message("[plotMethTrait] Saving file as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
}

annotateModule <- function(regions, module = NULL, grey = FALSE, genome = c("hg38", "hg19", "hg18", "mm10", "mm9", "danRer7"), 
                           includeCuratedRegDoms = FALSE, rule = c("basalPlusExt", "twoClosest", "oneClosest"), 
                           adv_upstream = 5, adv_downstream = 1, adv_span = 1000, adv_twoDistance = 1000, 
                           adv_oneDistance = 1000, version = c("4.0.4", "3.0.0", "2.0.2"), save = TRUE, 
                           file = "Annotated_Module_Regions.txt", verbose =  TRUE){
        if(!is.null(module)){
                if(verbose){
                        message("[annotateModule] Filtering for regions in ", paste(module, collapse = ", "), " module(s)")
                }
                regions <- regions[regions$module %in% module,]
        }
        if(!grey){
                if(verbose){
                        message("[annotateModule] Excluding regions in grey (unassigned) module")
                }
                regions <- regions[!regions$module == "grey",]
        }
        genome <- match.arg(genome)
        rule <- match.arg(rule)
        version <- match.arg(version)
        if((version == "4.0.4" & genome %in% c("hg18", "danRer7")) | (version == "3.0.0" & genome %in% c("hg38", "hg18")) |
           (version == "2.0.2" & genome %in% c("hg38", "mm10"))){
                stop("[annotateModule] The ", genome, " genome assembly is not supported for GREAT v", version)
        }
        if(verbose){
                message("[annotateModule] Using the ", genome, " genome assembly for annotations")
                message("[annotateModule] Adding genes to regions using GREAT v", version, " with the ", rule, " rule")
        }
        GR_regions <- with(regions, GRanges(chr, ranges = IRanges(start, end = end), RegionID = RegionID))
        job <- submitGreatJob(GR_regions, species = genome, includeCuratedRegDoms = includeCuratedRegDoms, rule = rule,
                              adv_upstream = adv_upstream, adv_downstream = adv_downstream, adv_span = adv_span,
                              adv_twoDistance = adv_twoDistance, adv_oneDistance = adv_oneDistance, request_interval = 0,
                              version = version)
        regions_genes <- suppressGraphics(plotRegionGeneAssociationGraphs(job, type = 1, request_interval = 0)) %>% 
                as.data.frame() %>%
                merge(x = regions[,c("RegionID", "chr", "start", "end")], 
                      y = .[,c("seqnames", "start", "end", "gene", "distTSS")], by.x = c("chr", "start", "end"), 
                      by.y = c("seqnames", "start", "end"), all = TRUE, sort = FALSE)
        regions_genes$RegionID <- factor(regions_genes$RegionID, levels = unique(regions_genes$RegionID))
        colnames(regions_genes) <- str_replace_all(colnames(regions_genes), pattern = c("gene" = "gene_symbol", 
                                                                                        "distTSS" = "distance_to_TSS"))
        if(verbose){
                message("[annotateModule] Adding gene info from BioMart")
        }
        dataset <- str_sub(genome, start = 1, end = 2) %>% switch(hg = "hsapiens_gene_ensembl", mm = "mmusculus_gene_ensembl",
                                                                  da = "drerio_gene_ensembl")
        ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset)
        genes_annotated <- suppressMessages(getBM(attributes = c("external_gene_name", "description", "ensembl_gene_id", "entrezgene_id"), 
                                                  filters = c("external_gene_name"), values = regions_genes$gene, mart = ensembl))
        colnames(genes_annotated) <- colnames(genes_annotated) %>%
                str_replace_all(pattern = c("external_gene_name" = "gene_symbol", "description" = "gene_description",
                                            "ensembl_gene_id" = "gene_ensemblID", "entrezgene_id" = "gene_entrezID"))
        genes_annotated$gene_description <- str_split_fixed(genes_annotated$gene_description, 
                                                            pattern = fixed(" ["), n = 2)[,1] %>%
                str_remove_all(",")
        regions_annotated <- merge(x = regions_genes, y = genes_annotated, by = "gene_symbol", all.x = TRUE, all.y = FALSE, 
                                   sort = FALSE) %>% 
                unique() %>%
                aggregate(formula = cbind(gene_symbol, distance_to_TSS, gene_description, gene_ensemblID, gene_entrezID) ~ RegionID, 
                          data = ., FUN = function(x) paste(x, collapse = " | "), simplify = TRUE, drop = FALSE) %>%
                merge(x = regions, y = ., by = "RegionID", all.x = TRUE, all.y = FALSE, sort = FALSE)
        if(verbose){
                message("[annotateModule] Getting gene context from annotatr")
        }
        annotations <- paste(genome, c("basicgenes", "genes_intergenic", "enhancers_fantom"), sep = "_")
        regions_GeneReg <- suppressWarnings(suppressMessages(build_annotations(genome = genome, annotations = annotations))) %>%
                GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
                annotate_regions(regions = GR_regions, annotations = ., ignore.strand = TRUE, quiet = TRUE) %>% 
                as.data.frame()
        colnames(regions_GeneReg)[colnames(regions_GeneReg) == "annot.type"] <- "gene_context"
        pattern <- rep("", times = 5)
        names(pattern) <- c(genome, "genes", "s", "fantom", "_")
        regions_GeneReg$gene_context <- str_replace_all(regions_GeneReg$gene_context, pattern = pattern)
        regions_annotated <- aggregate(formula = gene_context ~ RegionID, data = regions_GeneReg, 
                                       FUN = function(x) paste(unique(x), collapse = ", "), simplify = TRUE) %>%
                merge(x = regions_annotated, y = ., by = "RegionID", all.x = TRUE, all.y = FALSE, sort = FALSE) 
        if(verbose){
                message("[annotateModule] Getting CpG context from annotatr")
        }
        regions_CpGs <- suppressMessages(build_annotations(genome = genome, annotations = paste(genome, "cpgs", sep = "_"))) %>%
                GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
                annotate_regions(regions = GR_regions, annotations = ., ignore.strand = TRUE, quiet = TRUE) %>% 
                as.data.frame()
        colnames(regions_CpGs)[colnames(regions_CpGs) == "annot.type"] <- "CpG_context"
        pattern <- c("island", "shore", "shelf", "open sea")
        names(pattern) <- paste(genome, "cpg", c("islands", "shores", "shelves", "inter"), sep = "_")
        regions_CpGs$CpG_context <- str_replace_all(regions_CpGs$CpG_context, pattern = pattern)
        regions_annotated <- aggregate(formula = CpG_context ~ RegionID, data = regions_CpGs, 
                                       FUN = function(x) paste(unique(x), collapse = ", "), simplify = TRUE) %>%
                merge(x = regions_annotated, y = ., by = "RegionID", all.x = TRUE, all.y = FALSE, sort = FALSE) %>%
                .[order(.$module, as.integer(str_remove_all(.$RegionID, pattern = "Region_"))),]
        if(save){
                if(verbose){
                        message("[annotateModule] Saving file as ", file)
                }
                write.table(regions_annotated, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
        }
        return(regions_annotated)
}

getGeneList <- function(regions_annotated, module = NULL,type = c("symbol", "description", "ensemblID", "entrezID"), 
                        verbose = TRUE){
        if(!is.null(module)){
                if(verbose){
                        message("[getGeneList] Filtering for regions in ", paste(module, collapse = ", "), " module(s)")
                }
                regions_annotated <- regions_annotated[regions_annotated$module %in% module,]
        }
        type <- match.arg(type)
        if(verbose){
                message("[getGeneList] Getting gene ", type, "s")
        }
        geneList <- str_split(regions_annotated[,paste("gene", type, sep = "_")], pattern = fixed(" | ")) %>%
                as_vector() %>% unique() %>% sort()
        return(geneList)
}

enrichModule <- function(regions, module = NULL, genome = c("hg38", "hg19", "hg18", "mm10", "mm9", "danRer7"), 
                         includeCuratedRegDoms = FALSE, rule = c("basalPlusExt", "twoClosest", "oneClosest"), 
                         adv_upstream = 5, adv_downstream = 1, adv_span = 1000, adv_twoDistance = 1000, 
                         adv_oneDistance = 1000, version = c("4.0.4", "3.0.0", "2.0.2"), 
                         ontologies = c("GO Molecular Function", "GO Biological Process", "GO Cellular Component",
                                        "Mouse Phenotype", "Human Phenotype"), min_background_region_hits = 5,
                         adjMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"), 
                         min_region_hits = 2, pvalue_threshold = 0.01, save = TRUE, file = "Module_Enrichment_Results.txt", 
                         verbose =  TRUE){
        if(!is.null(module)){
                if(verbose){
                        message("[enrichModule] Analyzing functional enrichment for regions in ", 
                                paste(module, collapse = ", "), " module(s)")
                }
                moduleRegions <- regions[regions$module %in% module,]
        } else {
                if(verbose){
                        message("[enrichModule] Analyzing functional enrichment for all regions assigned to a module")
                }
                moduleRegions <- regions[!regions$module == "grey",]
        }
        genome <- match.arg(genome)
        rule <- match.arg(rule)
        version <- match.arg(version)
        if((version == "4.0.4" & genome %in% c("hg18", "danRer7")) | (version == "3.0.0" & genome %in% c("hg38", "hg18")) |
           (version == "2.0.2" & genome %in% c("hg38", "mm10"))){
                stop("[enrichModule] The ", genome, " genome assembly is not supported for GREAT v", version)
        }
        adjMethod <- match.arg(adjMethod)
        if(verbose){
                message("[enrichModule] Using GREAT v", version, " with the ", genome, " genome assembly and the ", rule, 
                        " rule")
        }
        if(genome == "danRer7"){
                ontologies <- ontologies[!ontologies %in% c("Mouse Phenotype", "Human Phenotype")]
        }
        gr <- with(moduleRegions, GRanges(chr, ranges = IRanges(start, end = end)))
        bg <- with(regions, GRanges(chr, ranges = IRanges(start, end = end)))
        job <- submitGreatJob(gr, bg = bg, species = genome, includeCuratedRegDoms = includeCuratedRegDoms, rule = rule,
                              adv_upstream = adv_upstream, adv_downstream = adv_downstream, adv_span = adv_span,
                              adv_twoDistance = adv_twoDistance, adv_oneDistance = adv_oneDistance, request_interval = 0,
                              version = version)
        if(verbose){
                message("[enrichModule] Getting results and adjusting p-values using the ", adjMethod, 
                        " method")
        }
        enrichTables <- getEnrichmentTables(job, ontology = ontologies)
        results <- list.rbind(enrichTables)
        rownames(results) <- 1:nrow(results)
        pattern <- c("hyper_" = "", "id" = "ID", "name" = "term", "total_regions" = "background_region_hits", 
                     "expected" = "expected_region_hits", "foreground_region_hits" = "region_hits", 
                     "region_set_coverage" = "region_coverage", "foreground_gene_hits" = "gene_hits", 
                     "total_genes_annotated" = "term_genes", "raw_pvalue" = "p")
        colnames(results) <- colnames(results) %>% str_to_lower() %>% str_replace_all(pattern = pattern)
        results$ontology <- rep(names(enrichTables), sapply(enrichTables, nrow))
        results <- subset(results, background_region_hits >= min_background_region_hits)
        results$log_p <- -log10(results$p)
        results$adj_p <- p.adjust(results$p, method = adjMethod)
        results$log_adj_p <- -log10(results$adj_p)
        results <- subset(results, region_hits >= min_region_hits & p < pvalue_threshold)
        if(verbose){
                message("[enrichModule] Extracting genes for enriched terms")
        }
        results$genes <- mapply(function(x,y){
                suppressGraphics(plotRegionGeneAssociationGraphs(job, type = 1, ontology = x, termID = y, request_interval = 0,
                                                                 max_tries = 10000, verbose = FALSE)) %>% .$gene %>% 
                        unique() %>% sort() %>% paste(collapse = ", ")
        }, x = results$ontology, y = results$ID, USE.NAMES = FALSE)
        results <- results[order(results$p), c("ID", "term", "ontology", "background_region_hits", "expected_region_hits", 
                                               "region_hits", "fold_enrichment", "region_coverage", "term_region_coverage", 
                                               "gene_hits", "genes", "background_gene_hits", "term_genes", "p", "log_p", 
                                               "adj_p", "log_adj_p")]
        if(save){
                if(verbose){
                        message("[enrichModule] Saving file as ", file)
                }
                write.table(results, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
        }
        return(results)
}

listOntologies <- function(genome = c("hg38", "hg19", "hg18", "mm10", "mm9", "danRer7"), 
                           version = c("4.0.4", "3.0.0", "2.0.2"), verbose =  TRUE){
        genome <- match.arg(genome)
        version <- match.arg(version)
        if((version == "4.0.4" & genome %in% c("hg18", "danRer7")) | (version == "3.0.0" & genome %in% c("hg38", "hg18")) |
           (version == "2.0.2" & genome %in% c("hg38", "mm10"))){
                stop("[listOntologies] The ", genome, " genome assembly is not supported for GREAT v", version)
        }
        if(verbose){
                message("[listOntologies] Getting available ontologies for GREAT v", version, " with the ", genome, 
                        " genome assembly")
        }
        ontologies <- GRanges("chr1", ranges = IRanges(1, end = 2)) %>%
                submitGreatJob(species = genome, request_interval = 0, version = version) %>% 
                availableOntologies()
        return(ontologies)
}

# Set Global Options ####
options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
enableWGCNAThreads(nThreads = 4)

# Read Bismark CpG Reports ####
colData <- read.xlsx("sample_info.xlsx", rowNames = TRUE)
bs <- getCpGs(colData, file = "Unfiltered_BSseq.rds")

# Examine CpG Totals at Different Cutoffs ####
CpGtotals <- getCpGtotals(bs, file = "CpG_Totals.txt")
plotCpGtotals(CpGtotals, file = "CpG_Totals.pdf")

# Filter BSobject ####
bs <- filterCpGs(bs, cov = 2, perSample = 0.75, file = "Filtered_BSseq.rds")

# Call Regions ####
regions <- getRegions(bs, file = "Unfiltered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99, file = "Unfiltered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99, file = "Unfiltered_SD_Plots.pdf")

# Examine Region Totals at Different Cutoffs ####
regionTotals <- getRegionTotals(regions, file = "Region_Totals.txt")
plotRegionTotals(regionTotals, file = "Region_Totals.pdf")

# Filter Regions ####
regions <- filterRegions(regions, covMin = 10, methSD = 0.05, file = "Filtered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99, file = "Filtered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99, file = "Filtered_SD_Plots.pdf")

# Adjust Methylation Data for PCs ####
meth <- getRegionMeth(regions, bs = bs, file = "Region_Methylation.rds")
mod <- model.matrix(~1, data = pData(bs))
methAdj <- adjustRegionMeth(meth, mod = mod, file = "Adjusted_Region_Methylation.rds")
getDendro(methAdj, distance = "euclidean") %>% plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))

# Select Soft Power Threshold ####
sft <- getSoftPower(methAdj, corType = "pearson", file = "Soft_Power.rds")
plotSoftPower(sft, file = "Soft_Power_Plots.pdf")

# Get Comethylation Modules ####
modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions, corType = "pearson", file = "Modules.rds")
plotRegionDendro(modules, file = "Region_Dendrograms.pdf")
BED <- getModuleBED(modules$regions, file = "Modules.bed")

# Examine Correlations between Modules and Samples ####
MEs <- modules$MEs
moduleDendro <- getDendro(MEs, distance = "bicor")
plotDendro(moduleDendro, labelSize = 4, nBreaks = 5, file = "Module_ME_Dendrogram.pdf")
moduleCor <- getCor(MEs, corType = "bicor")
plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro, file = "Module_Correlation_Heatmap.pdf")

sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
plotDendro(sampleDendro, labelSize = 3, nBreaks = 5, file = "Sample_ME_Dendrogram.pdf")
sampleCor <- getCor(MEs, transpose = TRUE, corType = "bicor")
plotHeatmap(sampleCor, rowDendro = sampleDendro, colDendro = sampleDendro, file = "Sample_Correlation_Heatmap.pdf")

plotHeatmap(MEs, rowDendro = sampleDendro, colDendro = moduleDendro, legend.title = "Module\nEigennode", 
            legend.position = c(0.37,0.89), file = "Sample_ME_Heatmap.pdf")

# Test Correlations between Module Eigennodes and Sample Traits ####
MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor", file = "ME_Trait_Correlation_Stats.txt")
traitDendro <- getCor(MEs, y = colData, corType = "bicor", robustY = FALSE) %>% getDendro(transpose = TRUE)
plotDendro(traitDendro, labelSize = 3.5, expandY = c(0.65,0.08), file = "Trait_Dendrogram.pdf")
plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order, traitOrder = traitDendro$order,
               file = "ME_Trait_Correlation_Heatmap.pdf")
plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order, traitOrder = traitDendro$order, sigOnly = TRUE, star.size = 11, 
               star.nudge_y = -0.27, legend.position = c(1.14, 0.745), colColorMargins = c(-1,5.1,0.5,10.47), 
               file = "Sig_ME_Trait_Correlation_Heatmap.pdf", width = 7, height = 3.5)

# Explore Significant ME-Trait Correlations ####
# Plot Module Eigennodes vs Traits
plotMEtraitDot(MEs$bisque4, trait = colData$Diagnosis_ASD, traitCode = c("TD" = 0, "ASD" = 1), 
               colors = c("TD" = "#3366CC", "ASD" = "#FF3366"), ylim = c(-0.2,0.2), xlab = "Diagnosis", 
               ylab = "Bisque 4 Module Eigennode", file = "bisque4_ME_Diagnosis_Dotplot.pdf")
plotMEtraitScatter(MEs$paleturquoise, trait = colData$Gran, ylim = c(-0.15,0.15), xlab = "Granulocytes", 
                   ylab = "Pale Turquoise Module Eigennode", file = "paleturquoise_ME_Granulocytes_Scatterplot.pdf")
plotMEtraitScatter(MEs$paleturquoise, trait = colData$Bcell, ylim = c(-0.15,0.15), xlab = "B-cells", 
                   ylab = "Pale Turquoise Module Eigennode", file = "paleturquoise_ME_Bcells_Scatterplot.pdf")

# Plot Region Methylation vs Traits
regions <- modules$regions
plotMethTrait("bisque4", regions = regions, meth = meth, trait = colData$Diagnosis_ASD, traitCode = c("TD" = 0, "ASD" = 1),
              traitColors = c("TD" = "#3366CC", "ASD" = "#FF3366"), trait.legend.title = "Diagnosis",
              file = "bisque4_Module_Methylation_Diagnosis_Heatmap.pdf")
plotMethTrait("paleturquoise", regions = regions, meth = meth, trait = colData$Gran, expandY = 0.04,
              trait.legend.title = "Granulocytes", trait.legend.position = c(1.034,3.35), 
              file = "paleturquoise_Module_Methylation_Granulocytes_Heatmap.pdf")
plotMethTrait("paleturquoise", regions = regions, meth = meth, trait = colData$Bcell, expandY = 0.04,
              trait.legend.title = "B-cells", trait.legend.position = c(1.004,3.35), 
              file = "paleturquoise_Module_Methylation_Bcells_Heatmap.pdf")

# Annotate Modules ####
regionsAnno <- annotateModule(regions, module = c("bisque4", "paleturquoise"), genome = "hg38",
                              file = "Annotated_bisque4_paleturquoise_Module_Regions.txt")
geneList_bisque4 <- getGeneList(regionsAnno, module = "bisque4")
geneList_paleturquoise <- getGeneList(regionsAnno, module = "paleturquoise")

# Analyze Functional Enrichment ####
ontologies <- listOntologies("hg38", version = "4.0.4")
enrich_bisque4 <- enrichModule(regions, module = "bisque4", genome = "hg38", file = "bisque4_Module_Enrichment.txt")
enrich_paleturquoise <- enrichModule(regions, module = "paleturquoise", genome = "hg38", 
                                     file = "paleturquoise_Module_Enrichment.txt")

