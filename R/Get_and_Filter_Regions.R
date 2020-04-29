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
        regions$chr <- factor(regions$chr, levels = seqlevels(bs))
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
                scale_y_continuous(expand = expansion(mult = c(0.008, 0.05))) +
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
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks), expand = expansion(mult = c(0.0062, 0.05))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks), expand = expansion(mult = c(0.006, 0.05))) +
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
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks), expand = expansion(mult = c(0.05, 0.03))) +
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

.regionTotals <- function(regions, covMin, methSD){
        regions <- regions[regions$covMin >= covMin & regions$methSD >= methSD,]
        totals <- c("covMin" = covMin, "methSD" = methSD, "totalRegions_K" = nrow(regions)/10^3,
                    "totalWidth_Mb" = sum(regions$width)/10^6, "totalN_M" = sum(regions$n)/10^6)
        return(totals)
}

