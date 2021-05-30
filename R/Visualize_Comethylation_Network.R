#' Plot Region Dendrograms
#'
#' \code{plotRegionDendro()} extracts plotting data from a \code{modules} list,
#' plots a region dendrogram with module assignments, and then saves it as a .pdf.
#'
#' \code{plotRegionDendro()} is designed to be used in combination with
#' \code{\link{getModules()}}. This function does not produce a \code{ggplot}
#' object, but instead uses the \code{WGCNA} function
#' \code{\link[WGCNA]{plotDendroAndColors()}} to plot the dendrogram.
#'
#' @param modules A \code{list} of module assignments and statistics produced by
#'         \code{\link{getModules()}}.
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
#' @return None, produces a plot as a side effect.
#'
#' @seealso \itemize{
#'         \item \code{\link{getModules()}} to build a comethylation network and
#'                 identify modules of comethylated regions.
#'         \item \code{\link{getModuleBED()}} to visualize genomic locations and
#'                 module assignments.
#'         \item \code{\link{getDendro()}} and \code{\link{plotDendro()}} to
#'                 generate and visualize dendrograms for samples, modules, and
#'                 traits.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Visualize Comethylation Modules
#' plotRegionDendro(modules, file = "Region_Dendrograms.pdf")
#' BED <- getModuleBED(modules$regions, file = "Modules.bed")
#' }
#'
#' @export
#'
#' @import WGCNA

plotRegionDendro <- function(modules, save = TRUE, file = "Region_Dendrograms.pdf",
                             width = 11, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotRegionDendro] Plotting region dendrograms and modules for each block")
        }
        blockColors <- lapply(modules$blockGenes, function(x) modules$colors[x])
        blockNames <- paste("Block ", 1:length(modules$dendrograms), " (",
                            sapply(modules$blockGenes, length), " regions)",
                            sep = "")
        if(save){
                if(verbose){
                        message("[plotRegionDendro] Saving plot as ", file)
                }
                pdf(file = file, width = width, height = height)
        }
        invisible(mapply(FUN = plotDendroAndColors, dendro = modules$dendrograms,
                         colors = blockColors, main = blockNames,
                         MoreArgs = list(groupLabels = "Modules", dendroLabels = FALSE,
                                         marAll = c(1.5,5,3,1.5), saveMar = FALSE,
                                         cex.lab = 1.2, cex.colorLabels = 1.2,
                                         autoColorHeight = FALSE, lwd = 0.8,
                                         colorHeight = 0.15, cex.axis = 1,
                                         frame.plot = TRUE)))
        invisible(dev.off())
}

#' Get a Module BED file
#'
#' \code{getModuleBED()} takes a \code{data.frame} of regions with module annotations,
#' converts it to the BED file format suitable for viewing it on the UCSC Genome
#' Browser, and then saves it.
#'
#' \code{getModuleBED()} is designed to be used in combination with
#' \code{\link{getModules()}}.The BED file produced includes a header line to
#' enable single-step viewing on the UCSC Genome Browser. Each region is labeled
#' by its \code{RegionID} and assigned module, and is colored by the module color.
#' "Grey" (unassigned) regions are excluded by default, but can be optionally
#' included.
#'
#' @param regions A \code{data.frame} of regions with module assignments, typically
#'         obtained from \code{\link{getModules()}}.
#' @param grey A \code{logical(1)} specifying whether to include "grey" (unassigned)
#'         regions in the BED file.
#' @param save A \code{logical(1)} indicating whether to save the BED file.
#' @param file A \code{character(1)} giving the file name (.BED).
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A BED file.
#'
#' @seealso \itemize{
#'         \item \code{\link{getModules()}} to build a comethylation network and
#'                 identify modules of comethylated regions.
#'         \item \code{\link{plotRegionDendro()}} to visualize region similarity
#'                 and module assignments.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Visualize Comethylation Modules
#' plotRegionDendro(modules, file = "Region_Dendrograms.pdf")
#' BED <- getModuleBED(modules$regions, file = "Modules.bed")
#' }
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_remove_all

getModuleBED <- function(regions, grey = FALSE, save = TRUE, file = "Modules.bed",
                         verbose = TRUE){
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
        BED <- cbind(regions[c("chr", "start", "end", "RegionID")], score = 0,
                     strand = ".", thickStart = 0, thickEnd = 0, rgb = regions$rgb)
        if(save){
                if(verbose){
                        message("[getModuleBED] Saving file as ", file)
                }
                name <- basename(file) %>% str_remove_all(pattern = ".bed")
                write(paste("track name='", name, "' description='", name,
                            "' itemRgb='On'", sep = ""), file = file)
                write.table(BED, file = file, append = TRUE, quote = FALSE,
                            sep = "\t", row.names = FALSE, col.names = FALSE)
        }
        return(BED)
}


#' Calculate Correlations
#'
#' \code{getCor()} calculates correlation coefficients using either \code{pearson}
#' or \code{bicor} methods. Calculations can be done between columns of a single
#' matrix or between two vectors or matrices.
#'
#' The first input argument can be optionally transposed. The correlation
#' calculations are performed by \code{WGCNA} functions \code{\link[WGCNA]{cor()}}
#' and \code{\link[WGCNA]{bicor()}}.
#'
#' @param x A \code{numeric vector} or \code{matrix}. \code{x} must be a \code{matrix}
#'         if \code{y} is null.
#' @param y A \code{numeric vector} or \code{matrix}. If null, correlations will
#'         be calculated for columns of \code{x}.
#' @param transpose A \code{logical(1)} specifying whether to transpose the
#'         \code{matrix}.
#' @param corType A \code{character(1)} indicating which correlation statistic
#'         to use in the calculation. Potential values include \code{pearson} and
#'         \code{bicor}.
#' @param maxPOutliers A \code{numeric(1)} specifying the maximum percentile that
#'         can be considered outliers on each side of the median for the
#'         \code{bicor} statistic.
#' @param robustY A \code{logical(1)} indicating whether to use robust calculation
#'         for \code{y} for the \code{bicor} statistic. \code{FALSE} is recommended
#'         if \code{y} is a binary variable.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{numeric matrix}.
#'
#' @seealso \itemize{
#'         \item \code{\link{getModules()}} to build a comethylation network and
#'                 identify modules of comethylated regions.
#'         \item \code{\link{getDendro()}} and \code{\link{plotDendro()}} to
#'                 generate and visualize dendrograms.
#'         \item \code{\link{plotHeatmap()}} to visualize correlations between
#'                 samples and modules.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Examine Correlations between Modules
#' MEs <- modules$MEs
#' moduleDendro <- getDendro(MEs, distance = "bicor")
#' plotDendro(moduleDendro, labelSize = 4, nBreaks = 5,
#'            file = "Module_ME_Dendrogram.pdf")
#' moduleCor <- getCor(MEs, corType = "bicor")
#' plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro,
#'             file = "Module_Correlation_Heatmap.pdf")
#'
#' # Examine Correlations between Samples
#' sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
#' plotDendro(sampleDendro, labelSize = 3, nBreaks = 5,
#'            file = "Sample_ME_Dendrogram.pdf")
#' sampleCor <- getCor(MEs, transpose = TRUE, corType = "bicor")
#' plotHeatmap(sampleCor, rowDendro = sampleDendro, colDendro = sampleDendro,
#'             file = "Sample_Correlation_Heatmap.pdf")
#'
#' # Visualize Module Eigennode Values
#' plotHeatmap(MEs, rowDendro = sampleDendro, colDendro = moduleDendro,
#'             legend.title = "Module\nEigennode",
#'             legend.position = c(0.37,0.89),
#'             file = "Sample_ME_Heatmap.pdf")
#' }
#'
#' @export
#'
#' @import WGCNA

getCor <- function(x, y = NULL, transpose = FALSE, corType = c("bicor", "pearson"),
                   maxPOutliers = 0.1, robustY = TRUE, verbose = TRUE){
        if(transpose){
                if(verbose){
                        message("[getCor] Transposing data")
                }
                x <- t(x)
        }
        corType <- match.arg(corType)
        if(verbose){
                message("[getCor] Calculating correlations using ", corType,
                        " correlation")
        }
        if(corType == "bicor"){
                cor <- bicor(x, y = y, use = "pairwise.complete.obs",
                             maxPOutliers = maxPOutliers, robustY = robustY,
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

#' Plot a Heatmap with Dendrograms
#'
#' \code{plotHeatmap()} takes a \code{numeric matrix} and two dendrograms, plots them,
#' and then saves it all as a .pdf. The row names and column names of \code{x}
#' must include the same set of values as the labels of \code{rowDendro} and
#' \code{colDendro}, respectively.
#'
#' \code{plotHeatmap()} is designed to be used in combination with
#' \code{\link{getCor()}} and \code{\link{getDendro()}}. The function will
#' check to see if module color names are in the row and column names and then
#' plot a color bar with the module colors. A \code{ggplot} object is produced
#' and can be edited outside of this function if desired.
#'
#' @param x A \code{numeric matrix}. The row names and column names of \code{x}
#'         must include the same set of values as the labels of \code{rowDendro}
#'         and \code{colDendro}, respectively.
#' @param rowDendro An \code{\link[stats]{hclust}} object generated by
#'         \code{\link{getDendro()}}.
#' @param colDendro An \code{\link[stats]{hclust}} object generated by
#'         \code{\link{getDendro()}}.
#' @param colors A \code{character} giving a vector of colors to use for the
#'         gradient on the heatmap. The default uses \code{\link[WGCNA]{blueWhiteRed()}}
#'         to generate these colors.
#' @param limit A \code{numeric(1)} giving the maximum value (symmetric) for the
#'         heatmap color scale.
#' @param axis.text.size A \code{numeric(1)} specifying the size of the text for
#'         both axes.
#' @param legend.title A \code{character(1)} with the title of the legend.
#' @param legend.title.size A \code{numeric(1)} giving the size of the text for
#'         the legend title.
#' @param legend.text.size A \code{numeric(1)} specifying the size of the text for
#'         the legend axis.
#' @param legend.position A \code{numeric(2)} with the position of the legend,
#'         as x-axis, y-axis. May also be a \code{character(1)} indicating "none",
#'         "left", "right", "bottom", or "top".
#' @param rowDendroMargins A \code{numeric(4)} giving the width of the margins
#'         for the row (vertical) dendrogram as top, right, bottom, and left.
#' @param colDendroMargins A \code{numeric(4)} giving the width of the margins
#'         for the column (horizontal) dendrogram as top, right, bottom, and left.
#' @param rowColorMargins A \code{numeric(4)} giving the width of the margins
#'         for the row (vertical) color bar as top, right, bottom, and left.
#' @param colColorMargins A \code{numeric(4)} giving the width of the margins
#'         for the column (horizontal) color bar as top, right, bottom, and left.
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
#'         \item \code{\link{getModules()}} to build a comethylation network and
#'                 identify modules of comethylated regions.
#'         \item \code{\link{getDendro()}} and \code{\link{plotDendro()}} to
#'                 generate and visualize dendrograms.
#'         \item \code{\link{getCor()}} to calculate correlation coefficients.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Examine Correlations between Modules
#' MEs <- modules$MEs
#' moduleDendro <- getDendro(MEs, distance = "bicor")
#' plotDendro(moduleDendro, labelSize = 4, nBreaks = 5,
#'            file = "Module_ME_Dendrogram.pdf")
#' moduleCor <- getCor(MEs, corType = "bicor")
#' plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro,
#'             file = "Module_Correlation_Heatmap.pdf")
#'
#' # Examine Correlations between Samples
#' sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
#' plotDendro(sampleDendro, labelSize = 3, nBreaks = 5,
#'            file = "Sample_ME_Dendrogram.pdf")
#' sampleCor <- getCor(MEs, transpose = TRUE, corType = "bicor")
#' plotHeatmap(sampleCor, rowDendro = sampleDendro, colDendro = sampleDendro,
#'             file = "Sample_Correlation_Heatmap.pdf")
#'
#' # Visualize Module Eigennode Values
#' plotHeatmap(MEs, rowDendro = sampleDendro, colDendro = moduleDendro,
#'             legend.title = "Module\nEigennode",
#'             legend.position = c(0.37,0.89),
#'             file = "Sample_ME_Heatmap.pdf")
#' }
#'
#' @export
#'
#' @import ggplot2
#' @import ggdendro
#' @import cowplot
#' @import WGCNA

plotHeatmap <- function(x, rowDendro, colDendro,
                        colors = blueWhiteRed(100, gamma = 0.3),
                        limit = max(abs(x)), axis.text.size = 8,
                        legend.title = "Bicor",  legend.title.size = 16,
                        legend.text.size = 12, legend.position = c(0.3,0.905),
                        rowDendroMargins = c(-1.55,1,-0.1,-1.1),
                        colDendroMargins = c(1,-0.5,-1,0.8),
                        rowColorMargins = c(-1.85,-1.5,0.55,1.7),
                        colColorMargins = c(-1.6,-0.85,1.8,0.55), save = TRUE,
                        file = "Heatmap.pdf", width = 11, height = 9.5,
                        verbose = TRUE){
        if(verbose){
                message("[plotHeatmap] Plotting heatmap with dendrograms")
        }
        limits <- c(-limit, limit)
        x <- as.data.frame(x)
        rownames(x) <- str_remove_all(rownames(x), pattern = "ME")
        colnames(x) <- str_remove_all(colnames(x), pattern = "ME")
        rowModules <- sum(rownames(x) %in% colors()) == length(rownames(x))
        colModules <- sum(colnames(x) %in% colors()) == length(colnames(x))
        x$rowID <- factor(rownames(x),
                          levels = rowDendro$labels[rev(rowDendro$order)])
        x <- reshape2::melt(x, id.vars = "rowID")
        x$variable <- factor(x$variable,
                             levels = colDendro$labels[colDendro$order])
        hmMarginL <- ifelse(rowModules, yes = 2, no = -1)
        hmMarginB <- ifelse(colModules, yes = 2, no = -1)
        heatmap <- ggplot(data = x) +
                geom_tile(aes(x = variable, y = rowID, color = value, fill = value)) +
                scale_fill_gradientn(legend.title, colors = colors, limits = limits,
                                     aesthetics = c("color", "fill")) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_text(size = axis.text.size,
                                                 color = "black", angle = 90,
                                                 vjust = 0.5),
                      axis.text.y = element_text(size = axis.text.size,
                                                 color = "black"),
                      axis.ticks = element_line(size = 0.8, color = "black"),
                      axis.title = element_blank(), legend.position = "none",
                      panel.background = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(), plot.background = element_blank(),
                      plot.margin = unit(c(0,1,hmMarginB,hmMarginL), "lines"))
        legend <- get_legend(heatmap + theme(legend.position = legend.position,
                                             legend.background = element_blank(),
                                             legend.title = element_text(size = legend.title.size),
                                             legend.text = element_text(size = legend.text.size)))
        rowDendroPlot <- ggplot(data = dendro_data(rowDendro)$segments) +
                geom_segment(aes(x = -x, y = y, xend = -xend, yend = yend),
                             lwd = 0.5, lineend = "square") +
                coord_flip() +
                theme_dendro() +
                theme(plot.margin = unit(rowDendroMargins, "lines"))
        colDendroPlot <- ggplot(data = dendro_data(colDendro)$segments) +
                geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                             lwd = 0.5, lineend = "square") +
                theme_dendro() +
                theme(plot.margin = unit(colDendroMargins, "lines"))
        rowColors <- NULL
        colColors <- NULL
        if(rowModules){
                if(verbose){
                        message("[plotHeatmap] Using colors in row names for y-axis labels")
                }
                rowColors <- ggplot(data = data.frame(x = 0, y = 1:length(levels(x$rowID)),
                                                      color = levels(x$rowID))) +
                        geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                        scale_fill_identity(aesthetics = c("color", "fill")) +
                        theme_void() +
                        theme(legend.position = "none",
                              plot.margin = unit(rowColorMargins, "lines"))
                heatmap <- heatmap + theme(axis.text.y = element_blank(),
                                           axis.ticks.y = element_blank())
        }
        if(colModules){
                if(verbose){
                        message("[plotHeatmap] Using colors in column names for x-axis labels")
                }
                colColors <- ggplot(data = data.frame(x = 1:length(levels(x$variable)),
                                                      y = 0, color = levels(x$variable))) +
                        geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                        scale_fill_identity(aesthetics = c("color", "fill")) +
                        theme_void() +
                        theme(legend.position = "none",
                              plot.margin = unit(colColorMargins, "lines"))
                heatmap <- heatmap + theme(axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank())
        }
        gg <- plot_grid(NULL, colDendroPlot, NULL, NULL, rowColors, heatmap,
                        rowDendroPlot, legend, NULL, colColors, NULL, NULL,
                        nrow = 3, ncol = 4, rel_widths = c(0.045, 1, 0.15, 0.15),
                        rel_heights = c(0.15, 1, 0.045))
        if(save){
                if(verbose){
                        message("[plotHeatmap] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width,
                       height = height, units = "in")
        }
        return(gg)
}

