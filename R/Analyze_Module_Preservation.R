#' Calculate Module Preservation
#'
#' \code{getModulePreservation()} examines module replication between a
#' reference and a test data set by estimating various preservation statistics,
#' which are then saved as a .txt file. Correlations are performed using either
#' \code{pearson} or \code{bicor} methods.
#'
#' Identical sets of regions should be assessed and assigned modules within
#' discovery (reference) and replication (test) data sets, though the replication
#' regions may be a subset of the discovery regions due to low coverage. It's
#' also recommended to filter CpGs so identical loci are also assessed within
#' regions. Network parameters should be as similar as possible, although
#' modules should be identified independently between the discovery and
#' replication datasets. Preservation statistics are calculated by
#' [WGCNA::modulePreservation()], with \code{corFnc} set to either \code{cor} or
#' \code{bicor}. More information is given in the documentation for
#' [WGCNA::modulePreservation()].
#'
#' @param meth_disc A \code{numeric matrix}, where each row is a sample and each
#'         column is a region. This is typically obtained from
#'         [adjustRegionMeth()] and is related to the discovery (reference) data
#'         set.
#' @param regions_disc A \code{data.frame} of regions with module assignments,
#'         which is typically obtained from [getModules()] and is related to the
#'         discovery data set.
#' @param meth_rep A \code{numeric matrix}, where each row is a sample and each
#'         column is a region. This is typically obtained from
#'         [adjustRegionMeth()] and is related to the replication (test) data
#'         set.
#' @param regions_rep A \code{data.frame} of regions with module assignments,
#'         which is typically obtained from [getModules()] and is related to the
#'         replication data set.
#' @param corType A \code{character(1)} indicating which correlation statistic
#'         to use in the adjacency calculation.
#' @param maxPOutliers A \code{numeric(1)} specifying the maximum percentile
#'         that can be considered outliers on each side of the median for the
#'         \code{bicor} statistic.
#' @param nPermutations A \code{numeric(1)} indicating the number of
#'         permutations to perform in the permutation test.
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{data.frame}.
#' @param file A \code{character(1)} giving the file name (.txt) for the saved
#'         \code{data.frame}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{data.frame} giving preservation statistics for each module in
#'         the discovery data set.
#'
#' @seealso \itemize{
#'         \item [getModules()] to build a comethylation network and identify
#'                 modules of comethylated regions.
#'         \item [plotModulePreservation()] to visualize module preservation
#'                 statistics.
#' }
#'
#' @examples \dontrun{
#' # Calculate Module Preservation
#' regions_disc <- modules_disc$regions
#' regions_rep <- modules_rep$regions
#' preservation <- getModulePreservation(methAdj_disc,
#'                                       regions_disc = regions_disc,
#'                                       meth_rep = methAdj_rep,
#'                                       regions_rep = regions_rep,
#'                                       corType = "pearson",
#'                                       file = "Module_Preservation_Stats.txt")
#'
#' # Visualize Module Preservation
#' plotModulePreservation(preservation, file = "Module_Preservation_Plots.pdf")
#' }
#'
#' @export
#'
#' @import WGCNA
#'
#' @importFrom magrittr %>%

getModulePreservation <- function(meth_disc, regions_disc, meth_rep, regions_rep,
                                  corType = c("pearson", "bicor"),
                                  maxPOutliers = 0.1, nPermutations = 100,
                                  save = TRUE,
                                  file = "Module_Preservation_Stats.txt",
                                  verbose = TRUE){
        corType <- match.arg(corType)
        if(!identical(colnames(meth_disc), regions_disc$RegionID)){
                stop("[getModulePreservation] RegionIDs must be identical between meth_disc and regions_disc")
        }
        if(!identical(colnames(meth_rep), regions_rep$RegionID)){
                stop("[getModulePreservation] RegionIDs must be identical between meth_rep and regions_rep")
        }
        if(!"module" %in% colnames(regions_disc) |
           !"module" %in% colnames(regions_rep)){
                stop("[getModulePreservation] Modules must be assigned to regions")
        }
        if(verbose){
                message("[getModulePreservation] Calculating module preservation statistics with ",
                        corType, " correlation")
                verboseNum <- 2
        } else {
                verboseNum <- 0
        }
        colnames(meth_disc) <- with(regions_disc,
                                    paste(chr, start, end, sep = "_"))
        colnames(meth_rep) <- with(regions_rep,
                                   paste(chr, start, end, sep = "_"))
        if(sum(!colnames(meth_rep) %in% colnames(meth_disc)) > 0){
                stop("[getModulePreservation] All replication regions must be included in discovery regions")
        }
        multiData <- list(Discovery = list(data = meth_disc),
                          Replication = list(data = meth_rep))
        multiColor <- list(Discovery = regions_disc$module,
                           Replication = regions_rep$module)
        if(corType == "pearson"){
                preservation <- modulePreservation(multiData,
                                                   multiColor = multiColor,
                                                   networkType = "signed",
                                                   corFnc = "cor",
                                                   nPermutations = nPermutations,
                                                   randomSeed = 5,
                                                   goldName = "random",
                                                   verbose = verboseNum)
        } else {
                if(corType == "bicor"){
                        preservation <- modulePreservation(multiData,
                                                           multiColor = multiColor,
                                                           networkType = "signed",
                                                           corFnc = "bicor",
                                                           corOptions = paste("use = 'p', maxPOutliers =",
                                                                              maxPOutliers),
                                                           nPermutations = nPermutations,
                                                           randomSeed = 5,
                                                           goldName = "random",
                                                           verbose = verboseNum)
                } else {
                        stop("[getModulePreservation] corType must be either pearson or bicor")
                }
        }
        categories <- c("quality", "preservation", "accuracy",
                        "referenceSeparability", "testSeparability")
        tables <- c("observed", "Z", "log.p")
        pattern <- c("Z." = "", "Z" = "", "log.p." = "", "log.p" = "")
        exclude <- c("medianRank.qual", "meanClusterCoeff.qual",
                     "medianRank.pres", "medianRankConnectivity.pres",
                     "medianRankDensity.pres", "cor.clusterCoeff",
                     "meanClusterCoeff.pres", "coClustering")
        preservationStats <- lapply(categories, function(x){
                category <- lapply(tables, function(y){
                        as.data.frame(preservation[[x]][[y]]$ref.Discovery$inColumnsAlsoPresentIn.Replication) %>%
                                tibble::rownames_to_column(var = "module") %>%
                                tidyr::pivot_longer(!module & !moduleSize,
                                                    names_to = "statistic",
                                                    values_to = y) %>%
                                dplyr::mutate(statistic = str_replace_all(statistic,
                                                                          pattern = fixed(pattern)))
                })
                names(category) <- tables
                category <- dplyr::full_join(category$observed, y = category$Z,
                                             by = c("module", "moduleSize",
                                                    "statistic")) %>%
                        dplyr::full_join(y = category$log.p,
                                         by = c("module", "moduleSize",
                                                "statistic")) %>%
                        dplyr::mutate(category = x)
        }) %>%
                dplyr::bind_rows() %>%
                dplyr::filter(!statistic %in% exclude) %>%
                dplyr::mutate(log.p = as.numeric(log.p)) %>%
                dplyr::select(module, moduleSize, category, statistic:log.p) %>%
                dplyr::arrange(module)
        if(save){
                if(verbose){
                        message("[getModulePreservation] Saving file as ", file)
                }
                write.table(preservationStats, file = file, quote = FALSE,
                            sep = "\t", row.names = FALSE)
        }
        return(preservationStats)
}

#' Visualize Module Preservation
#'
#' \code{plotModulePreservation()} plots Z-scores for various module
#' preservation statistics by module size as a scatterplot, and then saves it as
#' a .pdf.
#'
#' \code{plotModulePreservation()} is designed to be used in combination with
#' \code{getModulePreservation()}. A blue line is plotted at Z = 2 to indicate
#' weak to moderate evidence for preservation, while a green line is plotted at
#' Z = 10 to indicate strong evidence. A \code{ggplot} object is produced and
#' can be edited outside of this function if desired. More information is given
#' in the documentation for [WGCNA::modulePreservation()].
#'
#' @param preservation A \code{data.frame} of module preservation statistics for
#'         each module, generated by \code{getModulePreservation()}.
#' @param line.size A \code{numeric(1)} giving the size of the horizontal lines.
#' @param point.size A \code{numeric(1)} indicating the size of the points.
#' @param nBreaks A \code{numeric(1)} specifying the number of breaks used for
#'         both axes.
#' @param strip.text.size A \code{numeric(1)} indicating the size of the title
#'         text for each plot.
#' @param axis.text.size A \code{numeric(1)} specifying the size of the text for
#'         both axes.
#' @param axis.title.size A \code{numeric(1)} indicating the size of the title
#'         text for both axes.
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
#'         \item [getModules()] to build a comethylation network and identify
#'                 modules of comethylated regions.
#'         \item [getModulePreservation()] to calculate module preservation
#'                 statistics.
#' }
#'
#' @examples \dontrun{
#' # Calculate Module Preservation
#' regions_disc <- modules_disc$regions
#' regions_rep <- modules_rep$regions
#' preservation <- getModulePreservation(methAdj_disc,
#'                                       regions_disc = regions_disc,
#'                                       meth_rep = methAdj_rep,
#'                                       regions_rep = regions_rep,
#'                                       corType = "pearson",
#'                                       file = "Module_Preservation_Stats.txt")
#'
#' # Visualize Module Preservation
#' plotModulePreservation(preservation, file = "Module_Preservation_Plots.pdf")
#' }
#'
#' @export
#'
#' @import ggplot2
#'
#' @importFrom magrittr %>%
#' @importFrom scales breaks_pretty

plotModulePreservation <- function(preservation, line.size = 0.9,
                                   point.size = 1.2, nBreaks = 3,
                                   strip.text.size = 8.5, axis.text.size = 10,
                                   axis.title.size = 12, save = TRUE,
                                   file = "Module_Preservation_Plots.pdf",
                                   width = 9, height = 9, verbose = TRUE){
        if(verbose){
                message("[plotModulePreservation] Plotting module preservation statistic Z-scores by module size")
        }
        preservation <- dplyr::select(preservation, module, moduleSize,
                                      statistic, Z) %>%
                dplyr::mutate(statistic = factor(statistic,
                                                 levels = unique(statistic))) %>%
                dplyr::filter(!module %in% c("grey", "random") & is.finite(Z))
        scatterplot <- ggplot(preservation) +
                geom_hline(yintercept = 2, lty = 2, color = "blue",
                           size = line.size) +
                geom_hline(yintercept = 10, lty = 2, color = "darkgreen",
                           size = line.size) +
                geom_point(aes(x = moduleSize, y = Z, color = module),
                           size = point.size) +
                facet_wrap(vars(statistic), nrow = 6, scales = "free_y") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks),
                                   expand = expansion(0.15)) +
                scale_color_identity() +
                theme_bw(base_size = 25) +
                theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25),
                      axis.ticks = element_line(size = 0.9),
                      panel.grid.minor = element_blank(),
                      strip.background = element_blank(),
                      strip.text = element_text(margin = unit(c(1,0,0.2,0), "lines"),
                                                size = strip.text.size),
                      axis.text = element_text(color = "black",
                                               size = axis.text.size),
                      axis.title = element_text(size = axis.title.size),
                      panel.spacing.x = unit(0.2, "lines"),
                      panel.spacing.y = unit(-0.5, "lines"),
                      plot.margin = unit(c(0,1,1,1), "lines")) +
                xlab("Module Size") +
                ylab("Z-score")
        if(save){
                if(verbose){
                        message("[plotModulePreservation] Saving plots as ",
                                file)
                }
                ggsave(filename = file, plot = scatterplot, dpi = 600,
                       width = width, height = height, units = "in")
        }
        return(scatterplot)
}

