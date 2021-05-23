#' Estimate Soft Power Threshold
#'
#' \code{getSoftPower()} analyzes scale-free topology to estimate the best
#' soft-thresholding power from a vector of powers, calculate fit indices, and
#' then saves this as a .rds file. Possible correlation statistics include
#' \code{pearson} and \code{bicor}.
#'
#' Soft power is estimated by \code{\link[WGCNA]{pickSoftThreshold()}}, with
#' \code{corFnc} set to either \code{cor} or \code{bicor}. Calculations are
#' performed for a signed network in blocks of regions of size \code{blockSize}
#' (default = 40000). The best soft power threshold is chosen as the lowest power
#' where fit (R-squared) is greater than \code{RsquaredCut} (default = 0.8). More
#' information is given in the \code{WGCNA} package documentation for
#' \code{\link[WGCNA]{pickSoftThreshold()}}.
#'
#' @param meth A \code{numeric matrix}, where each row is a sample and each
#'         column is a region. This is typically obtained from
#'         \code{\link{adjustRegionMeth()}}.
#' @param powerVector A \code{numeric} specifying the soft power thresholds to
#'         examine for scale-free topology.
#' @param corType A \code{character(1)} indicating which correlation statistic
#'         to use in the adjacency calculation.
#' @param maxPOutliers A \code{numeric(1)} specifying the maximum percentile that
#'         can be considered outliers on each side of the median for the
#'         \code{bicor} statistic.
#' @param RsquaredCut A \code{numeric(1)} giving the minimum R-squared value for
#'         scale-free topology. Used to choose the best soft-thresholding power.
#' @param blockSize A \code{numeric(1)} specifying the number of regions in each
#'         block for the connectivity calculation. Decrease this if memory is
#'         insufficient.
#' @param gcInterval A \code{numeric(1)} indicating the interval for garbage
#'         collection.
#' @param save A \code{logical(1)} indicating whether to save the \code{list}.
#' @param file A \code{character(1)} giving the file name (.rds) for the saved
#'         \code{list}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{list} with two elements: \code{powerEstimate}, which gives the
#'         estimated best soft-thresholding power, and \code{fitIndices}, which
#'         is a \code{data.frame} with statistics on scale-free topology,
#'         including fit and connectivity, along with network density,
#'         centralization, and heterogeneity.
#'
#' @seealso \itemize{
#'         \item \code{\link{getRegionMeth()}} and \code{\link{adjustRegionMeth()}}
#'                 to extract methylation data and then adjust it for the top
#'                 principal components.
#'         \item \code{\link{plotSoftPower()}} to visualize fit and connectivity
#'                 for soft power estimation.
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
#' @import WGCNA

getSoftPower <- function(meth, powerVector = 1:20, corType = c("pearson", "bicor"),
                         maxPOutliers = 0.1, RsquaredCut = 0.8, blockSize = 40000,
                         gcInterval = blockSize - 1, save = TRUE,
                         file = "Soft_Power.rds", verbose = TRUE){
        corType <- match.arg(corType)
        if(verbose){
                message("[getSoftPower] Analyzing scale-free topology with ",
                        corType, " correlation to estimate best soft-thresholding power")
                verboseNum <- 1
        } else {
                verboseNum <- 0
        }
        if(corType == "pearson"){
                sft <- pickSoftThreshold(meth, RsquaredCut = RsquaredCut,
                                         powerVector = powerVector,
                                         networkType = "signed",
                                         moreNetworkConcepts = TRUE,
                                         corFnc = "cor",
                                         blockSize = blockSize,
                                         gcInterval = gcInterval,
                                         verbose = verboseNum)
        } else {
                if(corType == "bicor"){
                        sft <- pickSoftThreshold(meth, RsquaredCut = RsquaredCut,
                                                 powerVector = powerVector,
                                                 networkType = "signed",
                                                 moreNetworkConcepts = TRUE,
                                                 corFnc = "bicor",
                                                 corOptions = list(maxPOutliers = maxPOutliers),
                                                 blockSize = blockSize,
                                                 gcInterval = gcInterval,
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
                        ", fit = ",
                        round(sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == sft$powerEstimate], 3),
                        " and mean connectivity = ",
                        round(sft$fitIndices$mean.k.[sft$fitIndices$Power == sft$powerEstimate], 1))
        }
        if(save){
                if(verbose){
                        message("[getSoftPower] Saving file as ", file)
                }
                saveRDS(sft, file = file)
        }
        return(sft)
}

#' Plot Soft Power Fit and Connectivity
#'
#' \code{plotSoftPower()} visualizes scale-free topology fit and mean connectivity
#' for multiple soft power thresholds as a scatterplot, and then saves it as a
#' .pdf.
#'
#' \code{plotSoftPower()} is designed to be used in combination with
#' \code{\link{getSoftPower()}}. A \code{ggplot} object is produced and can be
#' edited outside of this function if desired.
#'
#' @param sft A \code{list} produced by \code{\link{getSoftPower()}} with two
#'         elements: \code{powerEstimate} and \code{fitIndices}.
#' @param pointCol A \code{character(1)} specifying the color of the points.
#' @param lineCol A \code{character(1)} giving the color of line and label for
#'         the estimated soft power threshold for scale-free topology.
#' @param nBreaks A \code{numeric(1)} specifying the number of breaks used for
#'         both axes.
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
#'         \item \code{\link{getRegionMeth()}} and \code{\link{adjustRegionMeth()}}
#'                 to extract methylation data and then adjust it for the top
#'                 principal components.
#'         \item \code{\link{getSoftPower()}} to calculate the best soft-thresholding
#'                 power and fit indices for scale-free topology.
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
#' @import ggplot2
#' @importFrom scales breaks_pretty
#' @importFrom stringr str_replace_all

plotSoftPower <- function(sft, pointCol = "#132B43", lineCol = "red", nBreaks = 4,
                          save = TRUE, file = "Soft_Power_Plots.pdf", width = 8.5,
                          height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotSoftPower] Plotting scale-free topology fit and mean connectivity by soft power threshold")
        }
        fitIndices <- data.frame(power = sft$fitIndices$Power,
                                 fit = -sign(sft$fitIndices[,"slope"]) *
                                         sft$fitIndices[,"SFT.R.sq"],
                                 log10_meanConnectivity = log10(sft$fitIndices$mean.k.),
                                 powerEstimate = sft$powerEstimate) %>%
                reshape2::melt(id.vars = c("power", "powerEstimate"))
        powerEstimateY <- min(0, fitIndices$value)
        fitIndices$variable <- str_replace_all(fitIndices$variable,
                                               c(fit = "Fit",
                                                 log10_meanConnectivity = "log[10]*(Mean~Connectivity)"))
        fitIndices$variable <- factor(fitIndices$variable,
                                      levels = c("Fit", "log[10]*(Mean~Connectivity)"))
        gg <- ggplot(data = fitIndices)
        gg <- gg +
                geom_vline(aes(xintercept = powerEstimate), color = lineCol) +
                geom_text(aes(x = powerEstimate, y = powerEstimateY,
                              label = powerEstimate), color = lineCol, nudge_x = -1) +
                geom_point(aes(x = power, y = value), color = pointCol, size = 1.2) +
                facet_wrap(vars(variable), nrow = 1, ncol = 2, scales = "free_y",
                           strip.position = "left", labeller = label_parsed) +
                xlab("Soft Power Threshold") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                expand_limits(x = 0, y = c(0,1)) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 12, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"),
                      axis.title.x = element_text(size = 16),
                      axis.title.y = element_blank(), legend.position = "none",
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(),
                      panel.spacing.x = unit(0.3, "lines"),
                      panel.spacing.y = unit(0.8, "lines"),
                      plot.margin = unit(c(1,1,0.7,0.2), "lines"),
                      strip.background = element_blank(), strip.placement = "outside",
                      strip.switch.pad.wrap = unit(0, "lines"),
                      strip.text.x = element_text(size = 16))
        if(save){
                if(verbose){
                        message("[plotSoftPower] Saving plots as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width,
                       height = height, units = "in")
        }
        return(gg)
}

#' Identify Modules of Comethylated Regions
#'
#' \code{getModules()} builds a comethylation network, identifies comethylated
#' modules, outputs a \code{list} with region module assignments, eigennode
#' values, dendrograms, and module membership, and then saves this as a .rds file.
#'
#' Comethylation networks are built and modules are identified by
#' \code{\link[WGCNA]{blockwiseModules()}}, with \code{corType} set to either
#' \code{pearson} or \code{bicor}. Calculations are performed for a signed network
#' in blocks of regions of maximum size \code{maxBlockSize} (default = 40000).
#' If there are more than \code{maxBlocksize} regions, then regions are
#' pre-clustered into blocks using projective K-means clustering. Region
#' correlations are performed within each block and regions are clustered with
#' average linkage hierarchical clustering. Modules are then identified with a
#' dynamic hybrid tree cut and highly correlated modules are merged together.
#' More information is given in the \code{WGCNA} package documentation for
#' \code{\link[WGCNA]{blockwiseModules()}}.
#'
#' @param meth A \code{numeric matrix}, where each row is a sample and each
#'         column is a region. This is typically obtained from
#'         \code{\link{adjustRegionMeth()}}.
#' @param power A \code{numeric(1)} giving the soft-thresholding power. This is
#'         typically obtained from \code{\link{getSoftPower()}}.
#' @param regions A \code{data.frame} of regions, typically after filtering with
#'         \code{\link{filterRegions()}}. Must have the column \code{RegionID}
#'         and correspond to the regions in \code{meth}.
#' @param maxBlockSize A \code{numeric(1)} specifying the maximum number of
#'         regions in a block. If there are more than this number regions, then
#'         regions are pre-clustered into blocks using projective K-means
#'         clustering. Decrease this if memory is insufficient.
#' @param corType A \code{character(1)} indicating which correlation statistic
#'         to use in the adjacency calculation.
#' @param maxPOutliers A \code{numeric(1)} specifying the maximum percentile that
#'         can be considered outliers on each side of the median for the
#'         \code{bicor} statistic.
#' @param deepSplit A \code{numeric(1)} specifying the sensitivity for module
#'         detection. Possible values are integers 0 to 4, with 4 having the
#'         highest sensitivity.
#' @param minModuleSize A \code{numeric(1)} giving the minimum number of regions
#'         to qualify as a module.
#' @param mergeCutHeight A \code{numeric(1)} specifying the cut height for
#'         merging correlated modules. Value is the maximum dissimilarity
#'         (1 - correlation) and ranges from 0 to 1.
#' @param nThreads A \code{numeric(1)} indicating the number of threads for
#'         correlation calculations.
#' @param save A \code{logical(1)} indicating whether to save the \code{list}.
#' @param file A \code{character(1)} giving the file name (.rds) for the saved
#'         \code{list}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{list} with 11 elements. See \code{\link[WGCNA]{blockwiseModules()}}
#'         for a description of these. Additional \code{regions} element is a
#'         \code{data.frame} with the region locations, statistics, module
#'         assignment, module membership, and hub region status.
#'
#' @seealso \itemize{
#'         \item \code{\link{getRegionMeth()}} and \code{\link{adjustRegionMeth()}}
#'                 to extract methylation data and then adjust it for the top
#'                 principal components.
#'         \item \code{\link{getSoftPower()}} and \code{\link{plotSoftPower()}}
#'                 to estimate the best soft-thresholding power and visualize
#'                 scale-free topology fit and connectivity.
#'         \item \code{\link{plotRegionDendro()}} and \code{\link{getModuleBED()}}
#'                 to visualize region similarity, genomic locations, and
#'                 module assignments.
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
#' # Select Soft Power Threshold
#' sft <- getSoftPower(methAdj, corType = "pearson", file = "Soft_Power.rds")
#' plotSoftPower(sft, file = "Soft_Power_Plots.pdf")
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
#' @importFrom magrittr %>%
#' @importFrom stringr str_remove_all

getModules <- function(meth, power, regions, maxBlockSize = 40000,
                       corType = c("pearson", "bicor"), maxPOutliers = 0.1,
                       deepSplit = 4, minModuleSize = 10, mergeCutHeight = 0.1,
                       nThreads = 4, save = TRUE, file = "Modules.rds",
                       verbose = TRUE){
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
                message("[getModules] Constructing network and detecting modules in blocks using ",
                        corType, " correlation")
                verboseNum <- 3
        } else {
                verboseNum <- 0
        }
        modules <- blockwiseModules(meth, checkMissingData = FALSE,
                                    maxBlockSize = maxBlockSize, corType = corType,
                                    maxPOutliers = maxPOutliers, power = power,
                                    networkType = "signed", TOMtype = "signed",
                                    deepSplit = deepSplit,
                                    minModuleSize = minModuleSize,
                                    mergeCutHeight = mergeCutHeight,
                                    nThreads = nThreads, verbose = verboseNum)
        if(verbose){
                message("[getModules] Assigning modules and calculating module membership using ",
                        corType, " correlation")
        }
        colnames(modules$MEs) <- str_remove_all(colnames(modules$MEs),
                                                pattern = "ME")
        if(corType == "pearson"){
                membership <- WGCNA::cor(x = meth, y = modules$MEs,
                                         use = "pairwise.complete.obs",
                                         nThreads = nThreads)
        } else {
                membership <- bicor(x = meth, y = modules$MEs,
                                    use = "pairwise.complete.obs",
                                    maxPOutliers = maxPOutliers,
                                    nThreads = nThreads)
        }
        regions$module <- modules$colors[match(regions$RegionID, names(modules$colors))]
        regions <- lapply(unique(regions$module), function(x){
                regions <- regions[regions$module == x,]
                regions$membership <- membership[regions$RegionID, x]
                regions$hubRegion <- regions$membership == max(regions$membership)
                return(regions)
        })
        regions <- list.rbind(regions) %>%
                .[order(as.integer(str_remove_all(.$RegionID, pattern = "Region_"))),]
        modules$regions <- regions
        if(save){
                if(verbose){
                        message("[getModules] Saving modules as ", file)
                }
                saveRDS(modules, file = file)
        }
        return(modules)
}

