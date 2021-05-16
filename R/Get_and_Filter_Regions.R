#' Generate Regions from CpGs
#'
#' \code{getRegions()} generates a set of regions and some statistics based on
#' the CpGs in a \code{\link[bsseq:`BSseq-class`]{BSseq}} object and then saves
#' it as a tab-delimited text file. Regions can be defined based on CpG
#' locations (for CpG clusters), built-in genomic annotations from
#' \pkg{annotatr}, or a custom genomic annotation.
#'
#' These regions still need to be filtered for minimum coverage and methylation
#' standard deviation.
#'
#' @param bs A \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#' @param annotation A \code{character(1)} giving the built-in genomic
#'         annotation to use for defining regions. Shortcuts are available for
#'         \code{genes}, \code{promoters}, and \code{transcripts}. Get the
#'         entire list of possible annotations with
#'         \code{\link[annotatr]{builtin_annotations()}}, which also includes
#'         CpG islands, enhancers, and chromatin states.
#' @param genome A \code{character(1)} with the genome build to use for
#'         built-in annotations. Available builds include \code{hg38},
#'         \code{hg19}, \code{mm10}, \code{mm9}, \code{rn6}, \code{rn5},
#'         \code{rn4}, \code{dm6}, \code{dm3}, and \code{galGal5}.
#' @param upstream A \code{numeric(1)} giving the number of bases upstream of a
#'         transcription start site to specify a promoter. Used for the
#'         \code{promoters} built-in annotation.
#' @param downstream A \code{numeric(1)} giving the number of bases downstream
#'         of a transcription start site to specify a promoter. Used for the
#'         \code{promoters} built-in annotation.
#' @param custom A \code{\link[GenomicRanges:`GRanges-class`]{GRanges}} object
#'         with a custom genomic annotation for defining regions. Construct this
#'         using \code{\link[GenomicRanges]{GRanges()}}.
#' @param maxGap A \code{numeric(1)} specifying the maximum number of bases
#'         between CpGs to be included in the same CpG cluster.
#' @param n A \code{numeric(1)} giving the minimum number of CpGs for a region
#'         to be returned. This applies to CpG clusters, built-in, and custom,
#'         annotations.
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{data.frame}.
#' @param file A \code{character(1)} giving the file name (.txt) for the saved
#'         \code{data.frame}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{data.frame} with the region genomic locations along with some
#'         statistics, including number of CpGs, coverage minimum, mean, and
#'         standard deviation, and methylation mean and standard deviation.
#'
#' @seealso \itemize{
#'         \item \code{\link{plotRegionStats()}}, \code{\link{plotSDstats()}},
#'                 \code{\link{getRegionTotals()}}, and
#'                 \code{\link{plotRegionTotals()}} for help visualizing region
#'                 characteristics and setting cutoffs for filtering.
#'        \item \code{\link{filterRegions()}} for filtering regions by minimum
#'                coverage and methylation standard deviation.
#' }
#'
#' @examples \dontrun{
#'
#' # Call Regions
#' regions <- getRegions(bs, file = "Unfiltered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Unfiltered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Unfiltered_SD_Plots.pdf")
#'
#' # Examine Region Totals at Different Cutoffs
#' regionTotals <- getRegionTotals(regions, file = "Region_Totals.txt")
#' plotRegionTotals(regionTotals, file = "Region_Totals.pdf")
#'
#' # Filter Regions
#' regions <- filterRegions(regions, covMin = 10, methSD = 0.05,
#'                          file = "Filtered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Filtered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Filtered_SD_Plots.pdf")
#' }
#'
#' @export
#'
#' @import bsseq
#' @importFrom magrittr %>%

getRegions <- function(bs, annotation = NULL,
                       genome = c("hg38", "hg19", "mm10", "mm9", "rn6", "rn5",
                                  "rn4", "dm6", "dm3", "galGal5"),
                       upstream = 5000, downstream = 1000, custom = NULL,
                       maxGap = 150, n = 3, save = TRUE,
                       file = "Unfiltered_Regions.txt", verbose = TRUE){
        if(!is.null(annotation) & !is.null(custom)){
                stop("[getRegions] annotation and custom cannot both have values")
        }
        if(!is.null(annotation)){
                genome <- match.arg(genome)
                if(verbose){
                        message("[getRegions] Using ", annotation,
                                " annotation for the ", genome,
                                " genome as regions")
                }
                if(annotation %in% c("genes", "promoters", "transcripts")){
                        txdb <- annotatr:::get_txdb_name(genome)
                        if(requireNamespace(txdb, quietly = TRUE)){
                                library(txdb, character.only = TRUE)
                        }
                        txdb <- get(txdb)
                        orgdb <- annotatr:::get_orgdb_name(genome)
                        if (requireNamespace(sprintf("org.%s.eg.db", orgdb),
                                             quietly = TRUE)) {
                                library(sprintf("org.%s.eg.db", orgdb),
                                        character.only = TRUE)
                        }
                        regions <- suppressWarnings(
                                switch(annotation,
                                       genes = GenomicFeatures::genes(txdb),
                                       promoters = GenomicFeatures::promoters(txdb,
                                                                              upstream = upstream,
                                                                              downstream = downstream,
                                                                              use.names = FALSE),
                                       transcripts = GenomicFeatures::transcripts(txdb))
                        )
                } else {
                        if(annotation %in% annotatr::builtin_annotations()){
                                if(!grepl(genome, x = annotation, fixed = TRUE)){
                                        stop("[getRegions] Annotation must match genome")
                                }
                                regions <- annotatr::build_annotations(genome,
                                                                       annotations = annotation)
                        } else {
                                stop("[getRegions] Annotation not supported")
                        }
                }
                regions <- GenomeInfoDb::keepStandardChromosomes(regions,
                                                                 pruning.mode = "coarse") %>%
                        trim()
                regions$n <- GenomicRanges::countOverlaps(regions, subject = bs)
                regions <- as.data.frame(regions)
                colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
        } else {
                if(!is.null(custom)){
                        if(verbose){
                                message("[getRegions] Using custom regions")
                        }
                        regions <- custom
                        regions$n <- GenomicRanges::countOverlaps(regions,
                                                                  subject = bs)
                        regions <- as.data.frame(regions)
                        colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
                } else {
                        if(verbose){
                                message("[getRegions] Calling regions when at least ",
                                        n, " CpGs are no more than ", maxGap,
                                        " bases apart")
                        }
                        regions <- bsseq:::regionFinder3(x = as.integer(rep(1, length(bs))),
                                                         chr = as.character(GenomeInfoDb::seqnames(bs)),
                                                         positions = BiocGenerics::start(bs),
                                                         maxGap = maxGap,
                                                         verbose = FALSE)[["up"]]
                }
        }
        regions$chr <- factor(regions$chr, levels = GenomeInfoDb::seqlevels(bs))
        regions <- regions[regions$n >= n,] %>%
                .[with(., order(chr, start, end)),]
        regions$RegionID <- paste("Region", 1:nrow(regions), sep = "_")
        regions$chr <- as.character(regions$chr)
        if(verbose){
                message("[getRegions] Calculating region statistics")
        }
        regions$width <- regions$end - regions$start
        cov <- getCoverage(bs, regions = regions[,c("chr", "start", "end")],
                           what = "perRegionTotal")
        regions$covMin <- DelayedArray::rowMins(cov)
        regions$covMean <- DelayedMatrixStats::rowMeans2(cov)
        regions$covSD <- DelayedMatrixStats::rowSds(cov)
        meth <- getMeth(bs, regions = regions[,c("chr", "start", "end")],
                        type = "raw", what = "perRegion")
        regions$methMean <- DelayedMatrixStats::rowMeans2(meth, na.rm = TRUE)
        regions$methSD <- DelayedMatrixStats::rowSds(meth, na.rm = TRUE)
        colnames <- c("RegionID", "chr", "start", "end", "width", "n", "covMin",
                      "covMean", "covSD", "methMean", "methSD")
        regions <- regions[,c(colnames,
                              colnames(regions)[!colnames(regions) %in% colnames &
                                                        !colnames(regions) %in%
                                                        c("idxStart", "idxEnd", "cluster")])]
        if(save){
                if(verbose){
                        message("[getRegions] Saving file as ", file)
                }
                write.table(regions, file = file, quote = FALSE, sep = "\t",
                            row.names = FALSE)
        }
        return(regions)
}

#' Plot Histograms of Region Statistics
#'
#' \code{plotRegionStats()} takes a set of regions from \code{\link{getRegions()}},
#' generates histograms of region characteristics, and saves it as a pdf.
#' Region-level statistics include width, number of CpGs, minimum coverage, mean
#' coverage, mean methylation, and methylation standard deviation.
#'
#' It's recommended to examine region characteristics before and after filtering.
#' The vertical line on each histogram indicates the median value for that feature.
#' A \code{ggplot} object is produced and can be edited outside of this function
#' if desired.
#'
#' @param regions A \code{data.frame} output from \code{\link{getRegions()}}
#'         giving the set of regions and statistics for each region.
#' @param maxQuantile A \code{numeric(1)} giving the maximum quantile of each
#'         feature to plot.
#' @param bins A \code{numeric(1)} specifying the number of bins in each histogram.
#' @param histCol A \code{character(1)} giving the color of the histogram.
#' @param lineCol A \code{character(1)} giving the color of the vertical median
#'         line.
#' @param nBreaks A \code{numeric(1)} specifying the number of breaks for the
#'         x-axis.
#' @param save A \code{logical(1)} indicating whether to save the plot.
#' @param file A \code{character(1)} giving the file name (.pdf) for the plot.
#' @param width A \code{numeric(1)} specifying the width in inches of the saved
#'         plot.
#' @param height A \code{numeric(1)} specifying the height in inches of the saved
#'         plot.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \itemize{
#'         \item \code{\link{getRegions()}} to generate the set of regions.
#'         \item \code{\link{plotSDstats()}}, \code{\link{getRegionTotals()}},
#'                 and \code{\link{plotRegionTotals()}} for more help visualizing
#'                 region characteristics and setting cutoffs for filtering.
#'        \item \code{\link{filterRegions()}} for filtering regions by minimum
#'                 coverage and methylation standard deviation.
#' }
#'
#' @examples \dontrun{
#'
#' # Call Regions
#' regions <- getRegions(bs, file = "Unfiltered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Unfiltered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Unfiltered_SD_Plots.pdf")
#'
#' # Examine Region Totals at Different Cutoffs
#' regionTotals <- getRegionTotals(regions, file = "Region_Totals.txt")
#' plotRegionTotals(regionTotals, file = "Region_Totals.pdf")
#'
#' # Filter Regions
#' regions <- filterRegions(regions, covMin = 10, methSD = 0.05,
#'                          file = "Filtered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Filtered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Filtered_SD_Plots.pdf")
#' }
#'
#' @export
#'
#' @import ggplot2
#' @importFrom scales breaks_pretty

plotRegionStats <- function(regions, maxQuantile = 1, bins = 30, histCol = "#132B43",
                            lineCol = "red", nBreaks = 4, save = TRUE,
                            file = "Region_Plots.pdf", width = 11, height = 8.5,
                            verbose = TRUE){
        if(verbose){
                message("[plotRegionStats] Plotting histograms of region statistics")
        }
        variables <- c("width", "n", "covMin", "covMean", "methMean", "methSD")
        medians <- data.frame(variable = factor(variables, levels = variables),
                              value = sapply(regions[,variables], FUN = median))
        if(maxQuantile < 1){
                if(verbose){
                        message("[plotRegionStats] Limiting x-axis to values in bottom ",
                                maxQuantile * 100, "% for width, n, covMin, and covMean")
                }
                regions$width[regions$width >= quantile(regions$width,
                                                        probs = maxQuantile)] <- NA
                regions$n[regions$n >= quantile(regions$n,
                                                probs = maxQuantile)] <- NA
                regions$covMin[regions$covMin >= quantile(regions$covMin,
                                                          probs = maxQuantile)] <- NA
                regions$covMean[regions$covMean >= quantile(regions$covMean,
                                                            probs = maxQuantile)] <- NA
        }
        regions <- reshape2::melt(regions[,c("RegionID", "width", "n", "covMin",
                                             "covMean", "methMean", "methSD")],
                                  id.vars = "RegionID")
        gg <- ggplot(data = regions)
        gg <- gg +
                geom_histogram(aes(x = value), bins = bins, fill = histCol,
                               color = histCol, na.rm = TRUE) +
                geom_vline(data = medians, aes(xintercept = value), color = lineCol) +
                facet_wrap(vars(variable), nrow = 2, ncol = 3, scales = "free",
                           strip.position = "bottom") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(expand = expansion(mult = c(0.008, 0.05))) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_text(size = 12, color = "black"),
                      axis.text.y = element_blank(),
                      axis.ticks.x = element_line(size = 1.25, color = "black"),
                      axis.ticks.y = element_blank(), axis.title = element_blank(),
                      legend.position = "none",
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(), panel.spacing.x = unit(0.6, "lines"),
                      panel.spacing.y = unit(0.8, "lines"),
                      plot.margin = unit(c(1,1,0.3,0.6), "lines"),
                      strip.background = element_blank(), strip.placement = "outside",
                      strip.switch.pad.wrap = unit(0, "lines"),
                      strip.text.x = element_text(size = 16))
        if(save){
                if(verbose){
                        message("[plotRegionStats] Saving plots as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width,
                       height = height, units = "in")
        }
        return(gg)
}

#' Plot Heatmaps of Region Standard Deviation vs Features
#'
#' \code{plotSDstats()} takes a set of regions from \code{\link{getRegions()}},
#' generates heatmaps of methylation standard deviation against region features,
#' and saves it as a pdf. Compared features include number of CpGs, minimum
#' coverage, mean coverage, and mean methylation.
#'
#' It's recommended examine these plots before and after filtering to ensure
#' removal of regions with high variability due to insufficient data. Plots are
#' heatmaps of 2D bin counts, with the color indicating the number of regions in
#' that bin on the log10 scale. A \code{ggplot} object is produced and can be
#' edited outside of this function if desired.
#'
#' @param regions A \code{data.frame} output from \code{\link{getRegions()}}
#'         giving the set of regions and statistics for each region.
#' @param maxQuantile A \code{numeric(1)} giving the maximum quantile of each
#'         feature to plot.
#' @param bins A \code{numeric(1)} specifying the number of bins for both axes
#'         in each heatmap.
#' @param nBreaks A \code{numeric(1)} specifying the number of breaks for both
#'         axes.
#' @param legend.position A \code{numeric(2)} specifying the position of the
#'         legend, as x-axis, y-axis. May also be a \code{character(1)}
#'         indicating "none", "left", "right", "bottom", or "top".
#' @param save A \code{logical(1)} indicating whether to save the plot.
#' @param file A \code{character(1)} giving the file name (.pdf) for the plot.
#' @param width A \code{numeric(1)} specifying the width in inches of the saved
#'         plot.
#' @param height A \code{numeric(1)} specifying the height in inches of the saved
#'         plot.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \itemize{
#'         \item \code{\link{getRegions()}} to generate the set of regions.
#'         \item \code{\link{plotRegionStats()}}, \code{\link{getRegionTotals()}},
#'                 and \code{\link{plotRegionTotals()}} for more help visualizing
#'                 region characteristics and setting cutoffs for filtering.
#'        \item \code{\link{filterRegions()}} for filtering regions by minimum
#'                 coverage and methylation standard deviation.
#' }
#'
#' @examples \dontrun{
#'
#' # Call Regions
#' regions <- getRegions(bs, file = "Unfiltered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Unfiltered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Unfiltered_SD_Plots.pdf")
#'
#' # Examine Region Totals at Different Cutoffs
#' regionTotals <- getRegionTotals(regions, file = "Region_Totals.txt")
#' plotRegionTotals(regionTotals, file = "Region_Totals.pdf")
#'
#' # Filter Regions
#' regions <- filterRegions(regions, covMin = 10, methSD = 0.05,
#'                          file = "Filtered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Filtered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Filtered_SD_Plots.pdf")
#' }
#'
#' @export
#'
#' @import ggplot2
#' @importFrom scales breaks_pretty

plotSDstats <- function(regions, maxQuantile = 1, bins = 30, nBreaks = 4,
                        legend.position = c(1.09,0.9), save = TRUE,
                        file = "SD_Plots.pdf", width = 8.5, height = 8.5,
                        verbose = TRUE){
        if(verbose){
                message("[plotSDstats] Plotting methylation SD vs region statistics")
        }
        if(maxQuantile < 1){
                if(verbose){
                        message("[plotSDstats] Limiting x-axis to values in bottom ",
                                maxQuantile * 100, "% for n, covMin, and covMean")
                }
                regions$n[regions$n >= quantile(regions$n,
                                                probs = maxQuantile)] <- NA
                regions$covMin[regions$covMin >= quantile(regions$covMin,
                                                          probs = maxQuantile)] <- NA
                regions$covMean[regions$covMean >= quantile(regions$covMean,
                                                            probs = maxQuantile)] <- NA
        }
        regions <- reshape2::melt(regions[,c("RegionID", "n", "covMin", "covMean",
                                             "methMean", "methSD")],
                                  id.vars = c("RegionID", "methSD"))
        gg <- ggplot(data = regions)
        gg <- gg +
                geom_bin2d(aes(x = value, y = methSD, color = ..count..),
                           bins = bins, na.rm = TRUE) +
                facet_wrap(vars(variable), nrow = 2, ncol = 2, scales = "free_x",
                           strip.position = "bottom") +
                scale_fill_continuous(name = "Count", trans = "log10") +
                scale_color_continuous(guide = FALSE, trans = "log10") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks),
                                   expand = expansion(mult = c(0.0062, 0.05))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks),
                                   expand = expansion(mult = c(0.006, 0.05))) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 12, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size = 16, color = "black"),
                      legend.background = element_blank(),
                      legend.position = legend.position,
                      legend.text = element_text(size = 12),
                      legend.title = element_text(size = 16),
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(),
                      panel.spacing.x = unit(1, "lines"),
                      panel.spacing.y = unit(0.8, "lines"),
                      plot.margin = unit(c(1,6,0.3,1), "lines"),
                      strip.background = element_blank(),
                      strip.placement = "outside",
                      strip.switch.pad.wrap = unit(0, "lines"),
                      strip.text.x = element_text(size = 16))
        if(save){
                if(verbose){
                        message("[plotSDstats] Saving plots as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width,
                       height = height, units = "in")
        }
        return(gg)
}

#' Get Region Totals at Different Cutoffs
#'
#' \code{getRegionTotals()} calculates the total number of regions, as well as
#' the total width and number of CpGs remaining in a region set after filtering
#' at different \code{covMin} and \code{methSD} cutoffs and then saves it as a
#' tab-delimited text file.
#'
#' The purpose of this function is to help balance cutoffs for minimum coverage
#' and methylation standard deviation and identify a robust set of variably
#' methylated regions. It's recommended to input multiple \code{covMin} and
#' \code{methSD} cutoffs for comparison. Computational resources are also a
#' consideration for network construction, with region sets of 250K or less
#' generally performing well.
#'
#' @param regions A \code{data.frame} output from \code{\link{getRegions()}}
#'         giving the set of regions and statistics for each region.
#' @param covMin A \code{numeric} specifying the minimum number of reads at CpGs
#'         in a region in any sample for that region to be included in the total.
#' @param methSD A \code{numeric} specifying the minimum methylation standard
#'         deviation for that region to be included in the total.
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{data.frame}
#' @param file A \code{character(1)} giving the file name (.txt) for the saved
#'         \code{data.frame}
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{data.frame} giving the total number of regions, width, and
#'         number of CpGs at all combinations of \code{covMin} and \code{methSD}.
#'
#' @seealso \itemize{
#'         \item \code{\link{getRegions()}} to generate the set of regions.
#'         \item \code{\link{plotRegionStats()}}, \code{\link{plotSDstats()}},
#'                 and \code{\link{plotRegionTotals()}} for more help visualizing
#'                 region characteristics and setting cutoffs for filtering.
#'        \item \code{\link{filterRegions()}} for filtering regions by minimum
#'                 coverage and methylation standard deviation.
#' }
#'
#' @examples \dontrun{
#'
#' # Call Regions
#' regions <- getRegions(bs, file = "Unfiltered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Unfiltered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Unfiltered_SD_Plots.pdf")
#'
#' # Examine Region Totals at Different Cutoffs
#' regionTotals <- getRegionTotals(regions, file = "Region_Totals.txt")
#' plotRegionTotals(regionTotals, file = "Region_Totals.pdf")
#'
#' # Filter Regions
#' regions <- filterRegions(regions, covMin = 10, methSD = 0.05,
#'                          file = "Filtered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Filtered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Filtered_SD_Plots.pdf")
#' }
#'
#' @export
#'
#' @importFrom magrittr %>%

getRegionTotals <- function(regions, covMin = seq(0,20,2), methSD = seq(0,0.1,0.01),
                            save = TRUE, file = "Region_Totals.txt",
                            verbose = TRUE){
        if(verbose){
                message("[getRegionTotals] Calculating region totals at specified covMin and methSD cutoffs")
        }
        length_covMin <- length(covMin)
        length_methSD <- length(methSD)
        covMin <- rep(covMin, each = length_methSD)
        methSD <- rep(methSD, times = length_covMin)
        regionTotals <- mapply(FUN = .regionTotals, covMin = covMin,
                               methSD = methSD,
                               MoreArgs = list(regions = regions)) %>%
                t() %>% as.data.frame()
        if(save){
                if(verbose){
                        message("[getRegionTotals] Saving file as ", file)
                        write.table(regionTotals, file = file, sep = "\t",
                                    row.names = FALSE)
                }
        }
        return(regionTotals)
}

#' Visualize Region Totals at Different Cutoffs
#'
#' \code{plotRegionTotals()} plots the total number of regions, width, and total
#' number of CpGs remaining after filtering by different combinations of
#' \code{covMin} and \code{methSD} in a line plot and then saves it as a .pdf.
#'
#' \code{plotRegionTotals()} is designed to be used in combination with
#' \code{\link{getRegionTotals()}}. A \code{ggplot} object is produced and can
#' be edited outside of this function if desired.
#'
#' @param regionTotals A \code{data.frame}, output from
#'         \code{\link{getRegionTotals()}}
#' @param nBreaks A \code{numeric(1)} specifying the number of breaks used for
#'         both axes and the legend.
#' @param legend.position A \code{numeric(2)} specifying the position of the
#'         legend, as x-axis, y-axis. May also be a \code{character(1)}
#'         indicating "none", "left", "right", "bottom", or "top".
#' @param save A \code{logical(1)} indicating whether to save the plot.
#' @param file A \code{character(1)} giving the file name (.pdf) for the saved
#'         plot.
#' @param width A \code{numeric(1)} specifying the width in inches of the saved
#'         plot.
#' @param height A \code{numeric(1)} specifying the height in inches of the
#'         saved plot.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \itemize{
#'         \item \code{\link{getRegions()}} to generate the set of regions.
#'         \item \code{\link{plotRegionStats()}}, \code{\link{plotSDstats()}},
#'                 and \code{\link{getRegionTotals()}} for more help visualizing
#'                 region characteristics and setting cutoffs for filtering.
#'        \item \code{\link{filterRegions()}} for filtering regions by minimum
#'                 coverage and methylation standard deviation.
#' }
#'
#' @examples \dontrun{
#'
#' # Call Regions
#' regions <- getRegions(bs, file = "Unfiltered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Unfiltered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Unfiltered_SD_Plots.pdf")
#'
#' # Examine Region Totals at Different Cutoffs
#' regionTotals <- getRegionTotals(regions, file = "Region_Totals.txt")
#' plotRegionTotals(regionTotals, file = "Region_Totals.pdf")
#'
#' # Filter Regions
#' regions <- filterRegions(regions, covMin = 10, methSD = 0.05,
#'                          file = "Filtered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Filtered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Filtered_SD_Plots.pdf")
#' }
#'
#' @export
#'
#' @import ggplot2
#' @importFrom scales breaks_pretty
#' @importFrom magrittr %>%

plotRegionTotals <- function(regionTotals, nBreaks = 4,
                             legend.position = c(1.08,0.897), save = TRUE,
                             file = "Region_Totals.pdf", width = 11, height = 11,
                             verbose = TRUE){
        if(verbose){
                message("[plotRegionTotals] Plotting region totals")
        }
        regionTotals <- reshape2::melt(regionTotals,
                                       id.vars = c("covMin", "methSD"))
        regionTotals$variable <- as.character(regionTotals$variable) %>%
                stringr::str_replace_all(pattern = c("totalRegions_K" = "Total Regions (Thousands)",
                                                     "totalWidth_Mb" = "Total Width (Mb)",
                                                     "totalN_M" = "Total CpGs (Millions)")) %>%
                factor(levels = c("Total Regions (Thousands)", "Total Width (Mb)",
                                  "Total CpGs (Millions)"))
        gg <- ggplot(data = regionTotals)
        gg <- gg +
                geom_line(aes(x = methSD, y = value, group = covMin,
                              color = covMin)) +
                geom_text(data = subset(regionTotals, methSD == min(methSD)),
                          aes(x = methSD, y = value, group = covMin,
                              color = covMin, label = covMin),
                          size = 4.5, check_overlap = TRUE, nudge_x = -0.001,
                          hjust = 1) +
                facet_wrap(vars(variable), nrow = 3, ncol = 1, scales = "free_y",
                           strip.position = "left") +
                xlab("SD Cutoff") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks),
                                   expand = expansion(mult = c(0.05, 0.03))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_color_gradient("Minimum\nCoverage\nCutoff",
                                     breaks = breaks_pretty(n = nBreaks - 1)) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 14, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"),
                      axis.title.x = element_text(size = 18),
                      axis.title.y = element_blank(),
                      legend.background = element_blank(),
                      legend.position = legend.position,
                      legend.title = element_text(size = 18),
                      legend.text = element_text(size = 14),
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(),
                      panel.spacing.x = unit(0.3, "lines"),
                      panel.spacing.y = unit(0.8, "lines"),
                      plot.margin = unit(c(1,7,0.7,0.2), "lines"),
                      strip.background = element_blank(),
                      strip.placement = "outside",
                      strip.switch.pad.wrap = unit(0, "lines"),
                      strip.text = element_text(size = 18))
        if(save){
                if(verbose){
                        message("[plotRegionTotals] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width,
                       height = height, units = "in")
        }
        return(gg)
}

#' Filter Regions
#'
#' \code{filterRegions()} filters a region set using specified \code{covMin} and
#' \code{methSD} thresholds and then saves it as a tab-delimited text file.
#'
#' The purpose of this function is to use cutoffs for minimum coverage
#' and methylation standard deviation guided by other functions and obtain a
#' robust set of variably methylated regions. Computational resources are also a
#' consideration for network construction, with region sets of 250K or less
#' generally performing well.
#'
#' @param regions A \code{data.frame} output from \code{\link{getRegions()}}
#'         giving the set of regions and statistics for each region.
#' @param covMin A \code{numeric(1)} specifying the minimum number of reads at CpGs
#'         in a region in any sample
#' @param methSD A \code{numeric(1)} specifying the minimum methylation standard
#'         deviation in a region
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{data.frame}
#' @param file A \code{character(1)} giving the file name (.txt) for the saved
#'         \code{data.frame}
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A filtered version of the \code{regions} \code{data.frame}
#'
#' @seealso \itemize{
#'         \item \code{\link{getRegions()}} to generate the set of regions.
#'         \item \code{\link{plotRegionStats()}}, \code{\link{plotSDstats()}},
#'                 \code{\link{getRegionTotals()}}, and
#'                 \code{\link{plotRegionTotals()}} for help visualizing region
#'                 characteristics and setting cutoffs for filtering.
#'         \item \code{\link{getRegionMeth()}} to get methylation values for
#'                 these regions in all samples.
#' }
#'
#' @examples \dontrun{
#'
#' # Call Regions
#' regions <- getRegions(bs, file = "Unfiltered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Unfiltered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Unfiltered_SD_Plots.pdf")
#'
#' # Examine Region Totals at Different Cutoffs
#' regionTotals <- getRegionTotals(regions, file = "Region_Totals.txt")
#' plotRegionTotals(regionTotals, file = "Region_Totals.pdf")
#'
#' # Filter Regions
#' regions <- filterRegions(regions, covMin = 10, methSD = 0.05,
#'                          file = "Filtered_Regions.txt")
#' plotRegionStats(regions, maxQuantile = 0.99,
#'                 file = "Filtered_Region_Plots.pdf")
#' plotSDstats(regions, maxQuantile = 0.99,
#'             file = "Filtered_SD_Plots.pdf")
#' }
#'
#' @export

filterRegions <- function(regions, covMin = 10, methSD = 0.05, save = TRUE,
                          file = "Filtered_Regions.txt", verbose = TRUE){
        if(verbose){
                message("[getRegions] Filtering regions for at least ", covMin,
                        " reads in all samples and methylation SD of at least ",
                        methSD * 100, "%")
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
                write.table(regions, file = file, quote = FALSE, sep = "\t",
                            row.names = FALSE)
        }
        return(regions)
}

.regionTotals <- function(regions, covMin, methSD){
        regions <- regions[regions$covMin >= covMin & regions$methSD >= methSD,]
        totals <- c("covMin" = covMin, "methSD" = methSD,
                    "totalRegions_K" = nrow(regions)/10^3,
                    "totalWidth_Mb" = sum(regions$width)/10^6,
                    "totalN_M" = sum(regions$n)/10^6)
        return(totals)
}

