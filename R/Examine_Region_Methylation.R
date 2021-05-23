#' Get Region Methylation Data
#'
#' \code{getRegionMeth()} extracts methylation values at specified regions for
#' all samples and then saves it as a .rds file.
#'
#' Methylation is summarized at the region level, and is estimated as the
#' methylated reads divided by the total reads. Methylation values are obtained
#' from a \code{\link[bsseq:`BSseq-class`]{BSseq}} object and can be either raw
#' or smoothed methylation.
#'
#' @param regions A \code{data.frame} of regions, typically after filtering with
#'         \code{\link{filterRegions()}}. Must have the columns \code{chr},
#'         \code{start}, and \code{end}.
#' @param bs A \code{\link[bsseq:`BSseq-class`]{BSseq}} object, typically after
#'         filtering with \code{\link{filterCpGs()}}.
#' @param type A \code{character(1)} specifying the type of methylation values
#'         to extract. Accepted values are \code{raw} and \code{smooth}
#' @param save A \code{logical(1)} indicating whether to save the \code{matrix}.
#' @param file A \code{character(1)} giving the file name (.rds) for the saved
#'         \code{matrix}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{numeric matrix}, where each row is a region and each column
#'         is a sample.
#'
#' @seealso \itemize{
#'         \item \code{\link{adjustRegionMeth()}} to adjust methylation for the
#'                 top principal components.
#'         \item \code{\link{getDendro()}} and \code{\link{plotDendro()}} to
#'                 generate and visualize dendrograms.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Methylation Data
#' meth <- getRegionMeth(regions, bs = bs, file = "Region_Methylation.rds")
#'
#' # Adjust Methylation Data for PCs
#' mod <- model.matrix(~1, data = pData(bs))
#' methAdj <- adjustRegionMeth(meth, mod = mod,
#'                             file = "Adjusted_Region_Methylation.rds")
#'
#' # Assess Sample Similarity
#' getDendro(methAdj, distance = "euclidean") %>%
#'         plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))
#' }
#'
#' @export
#'
#' @import bsseq

getRegionMeth <- function(regions, bs, type = c("raw", "smooth"), save = TRUE,
                          file = "Region_Methylation.rds", verbose = TRUE){
        if(verbose){
                message("[getRegionMeth] Calculating region methylation from BSseq object")
        }
        type <- match.arg(type)
        meth <- getMeth(bs, regions = regions[,c("chr", "start", "end")],
                        type = type, what = "perRegion")
        rownames(meth) <- regions$RegionID
        if(save){
                if(verbose){
                        message("[getRegionMeth] Saving file as ", file)
                }
                saveRDS(meth, file = file)
        }
        return(meth)
}

#' Adjust Methylation Data for Principal Components
#'
#' \code{adjustRegionMeth()} adjusts region methylation data for the top
#' principal components, transposes it, and then saves it as a .rds file.
#'
#' \code{adjustRegionMeth()} uses \code{\link[sva]{sva_network()}} to regress out
#' the top principal components. More information on the function and approach is
#' given in the documentation and publications related to the \code{sva} package.
#'
#' @param meth A \code{numeric matrix}, where each row is a region and each
#'         column is a sample. This is typically obtained from
#'         \code{\link{getRegionMeth()}}.
#' @param mod A \code{matrix} giving the model matrix being used to fit the data.
#'         See below for an example.
#' @param save A \code{logical(1)} indicating whether to save the \code{matrix}.
#' @param file A \code{character(1)} giving the file name (.rds) for the saved
#'         \code{matrix}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{numeric matrix}, where each row is a sample and each column
#'         is a region.
#'
#' @seealso \itemize{
#'         \item \code{\link{getRegionMeth()}} to extract region methylation
#'                 values.
#'         \item \code{\link{getDendro()}} and \code{\link{plotDendro()}} to
#'                 generate and visualize dendrograms.
#'         \item \code{\link{getSoftPower()}} and \code{\link{plotSoftPower()}}
#'                 to estimate the best soft-thresholding power and visualize
#'                 scale-free topology fit and connectivity.
#'         \item \code{\link{getModules()}} to build a comethylation network and
#'                 identify modules of comethylated regions.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Methylation Data
#' meth <- getRegionMeth(regions, bs = bs, file = "Region_Methylation.rds")
#'
#' # Adjust Methylation Data for PCs
#' mod <- model.matrix(~1, data = pData(bs))
#' methAdj <- adjustRegionMeth(meth, mod = mod,
#'                             file = "Adjusted_Region_Methylation.rds")
#'
#' # Assess Sample Similarity
#' getDendro(methAdj, distance = "euclidean") %>%
#'         plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))
#'
#' # Select Soft Power Threshold
#' sft <- getSoftPower(methAdj, corType = "pearson", file = "Soft_Power.rds")
#' plotSoftPower(sft, file = "Soft_Power_Plots.pdf")
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#' }
#'
#' @export
#'
#' @import sva
#' @importFrom magrittr %>%

adjustRegionMeth <- function(meth, mod = matrix(1, nrow = ncol(meth), ncol = 1),
                             save = TRUE,
                             file = "Adjusted_Region_Methylation.rds",
                             verbose = TRUE){
        if(verbose){
                message("[adjustRegionMeth] Determining number of principal components to adjust for")
        }
        n.pc <- num.sv(meth, mod = mod, seed = 5)
        if(verbose){
                message("[adjustRegionMeth] Adjusting region methylation for the top ",
                        n.pc, " principal components")
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

#' Perform Hierarchical Clustering
#'
#' \code{getDendro()} computes the distance between the rows of a matrix and
#' performs hierarchical clustering. Possible distance measures include
#' \code{euclidean}, \code{pearson}, and \code{bicor}. The function also optionally
#' transposes the matrix.
#'
#' Euclidean distance is calculated by \code{\link[stats]{dist(method = "euclidean")}},
#' while Pearson correlation and biweight midcorrelation (bicor) are computed by
#' the \code{WGCNA} functions \code{\link[WGCNA]{cor}} and \code{\link[WGCNA]{bicor}},
#' respectively. The \code{cor} and \code{bicor} are then subtracted from 1 to
#' calculate the dissimilarity. Hierarchical clustering is done by
#' \code{\link[stats]{hclust(method = average)}}.
#'
#' @param x A \code{numeric matrix}.
#' @param transpose A \code{logical(1)} specifying whether to transpose the
#'         \code{matrix}.
#' @param distance A \code{character(1)} indicating which distance measure to use.
#'         Possible values include \code{euclidean}, \code{pearson}, and \code{bicor}.
#' @param maxPOutliers A \code{numeric(1)} specifying the maximum percentile that
#'         can be considered outliers on each side of the median for the
#'         \code{bicor} statistic.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return An \code{\link[stats]{hclust}} object that describes the clustering
#'         tree.
#'
#' @seealso \itemize{
#'         \item \code{\link{plotDendro()}} to visualize dendrograms from
#'                 \code{getDendro()}.
#' }
#'
#' @examples \dontrun{
#'
#' # Assess Sample Similarity
#' getDendro(methAdj, distance = "euclidean") %>%
#'         plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))
#'
#' # Examine Correlations between Modules
#' moduleDendro <- getDendro(MEs, distance = "bicor")
#' plotDendro(moduleDendro, labelSize = 4, nBreaks = 5,
#'            file = "Module_ME_Dendrogram.pdf")
#'
#' # Characterize Correlations between Samples
#' sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
#' plotDendro(sampleDendro, labelSize = 3, nBreaks = 5,
#'            file = "Sample_ME_Dendrogram.pdf")
#'
#' # Examine Correlations between Traits
#' traitDendro <- getCor(MEs, y = colData, corType = "bicor",
#'                       robustY = FALSE) %>%
#'         getDendro(transpose = TRUE)
#' plotDendro(traitDendro, labelSize = 3.5, expandY = c(0.65,0.08),
#'            file = "Trait_Dendrogram.pdf")
#' }
#'
#' @export
#'
#' @import WGCNA
#' @importFrom magrittr %>%
#' @importFrom stringr str_remove_all

getDendro <- function(x, transpose = FALSE,
                      distance = c("euclidean", "pearson", "bicor"),
                      maxPOutliers = 0.1, verbose = TRUE){
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
                        dist <- (1 - WGCNA::cor(x, use = "pairwise.complete.obs")) %>%
                                as.dist()
                } else {
                        if(distance == "bicor"){
                                if(verbose){
                                        message("[getDendro] Clustering with bicor correlation as the distance")
                                }
                                dist <- (1 - bicor(x, maxPOutliers = maxPOutliers,
                                                   use = "pairwise.complete.obs")) %>%
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

#' Plot a Dendrogram
#'
#' \code{plotDendro()} extracts plotting data from an \code{\link[stats]{hclust}}
#' object, plots a dendrogram, and then saves it as a .pdf.
#'
#' \code{plotDendro()} is designed to be used in combination with
#' \code{\link{getDendro()}}. A \code{ggplot} object is produced and can be
#' edited outside of this function if desired.
#'
#' @param dendro An \code{\link[stats]{hclust}} object generated by
#'         \code{\link{getDendro()}}.
#' @param label A \code{logical(1)} indicating whether to add labels to the
#'         dendrogram.
#' @param labelSize A \code{numeric(1)} specifying the size of the label text.
#' @param expandX A \code{numeric(2)} giving the x-axis multiplicative range
#'         expansion factors as upper limit, lower limit.
#' @param expandY A \code{numeric(2)} giving the y-axis multiplicative range
#'         expansion factors as upper limit, lower limit.
#' @param nBreaks A \code{numeric(1)} specifying the number of breaks used for
#'         the y-axis.
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
#'         \item \code{\link{getDendro()}} to generate dendrograms for
#'                 \code{plotDendro()}.
#' }
#'
#' @examples \dontrun{
#'
#' # Assess Sample Similarity
#' getDendro(methAdj, distance = "euclidean") %>%
#'         plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))
#'
#' # Examine Correlations between Modules
#' moduleDendro <- getDendro(MEs, distance = "bicor")
#' plotDendro(moduleDendro, labelSize = 4, nBreaks = 5,
#'            file = "Module_ME_Dendrogram.pdf")
#'
#' # Characterize Correlations between Samples
#' sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
#' plotDendro(sampleDendro, labelSize = 3, nBreaks = 5,
#'            file = "Sample_ME_Dendrogram.pdf")
#'
#' # Examine Correlations between Traits
#' traitDendro <- getCor(MEs, y = colData, corType = "bicor",
#'                       robustY = FALSE) %>%
#'         getDendro(transpose = TRUE)
#' plotDendro(traitDendro, labelSize = 3.5, expandY = c(0.65,0.08),
#'            file = "Trait_Dendrogram.pdf")
#' }
#'
#' @export
#'
#' @import ggplot2
#' @import ggdendro
#' @importFrom scales breaks_pretty
#' @importFrom stringr str_remove_all

plotDendro <- function(dendro, label = TRUE, labelSize = 2.5,
                       expandX = c(0.03,0.03), expandY = c(0.3,0.08), nBreaks = 4,
                       save = TRUE, file = "Dendrogram.pdf", width = 11,
                       height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotDendro] Plotting dendrogram")
        }
        dendroPlot <- dendro_data(dendro)
        fix <- dendroPlot$segments$yend == 0
        dendroPlot$segments$yend[fix] <- dendroPlot$segments$y[fix] -
                max(dendroPlot$segments$y) * 0.05
        dendroPlot$labels$y <- dendroPlot$segments$yend[fix] -
                max(dendroPlot$segments$y) * 0.01
        dendroPlot$labels$label <- str_remove_all(dendroPlot$labels$label,
                                                  pattern = "ME")
        gg <- ggplot()
        gg <- gg +
                geom_segment(data = dendroPlot$segments,
                             aes(x = x, y = y, xend = xend, yend = yend),
                             lwd = 0.3, lineend = "square") +
                scale_x_continuous(expand = expansion(mult = expandX)) +
                scale_y_continuous(expand = expansion(mult = expandY),
                                   breaks = breaks_pretty(n = nBreaks)) +
                ylab("Height") +
                theme_dendro() +
                theme(plot.margin = unit(c(1,1,0,1), "lines"),
                      panel.background = element_rect(color = "black",
                                                      fill = "white", size = 1.1),
                      axis.ticks.y = element_line(),
                      axis.text.y = element_text(size = 12, color = "black"),
                      axis.title.y = element_text(size = 16, angle = 90,
                                                  vjust = 2))
        if(label){
                gg <- gg +
                        geom_text(data = dendroPlot$labels,
                                  aes(x = x, y = y, label = label), angle = 90,
                                  hjust = 1, size = labelSize)
        }
        if(save){
                if(verbose){
                        message("[plotDendro] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width,
                       height = height, units = "in")
        }
        return(gg)
}

