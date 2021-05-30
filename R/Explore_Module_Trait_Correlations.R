#' Calculate Correlation Statistics Between Module Eigennodes and Traits
#'
#' \code{getMEtraitCor()} calculates correlation coefficients and p-values between
#' eigennode values for all modules and all sample traits and saves it as a .txt
#' file. Correlations are performed using either \code{pearson} or \code{bicor}
#' methods.
#'
#' \code{getMEtraitCor()} is designed to be used in combination with
#' \code{\link{getModules()}}. The correlation calculations are performed by
#' \code{WGCNA} functions \code{\link[WGCNA]{corAndPvalue()}} and
#' \code{\link[WGCNA]{bicorAndPvalue()}}. P-values can be adjusted for multiple
#' comparisons by several methods, including \code{fdr}, \code{holm},
#' \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH}, and \code{BY}.
#'
#' @param MEs A \code{data.frame} of module eigennode values, where rows are
#'         samples and columns are modules. The row names of \code{MEs} must
#'         include the same set of values as the row names of \code{colData}.
#' @param colData A \code{data.frame} whose row names specify samples and whose
#'         columns are sample traits with numeric values.
#' @param corType A \code{character(1)} indicating which correlation statistic
#'         to use in the calculation. Potential values include \code{pearson} and
#'         \code{bicor}.
#' @param maxPOutliers A \code{numeric(1)} specifying the maximum percentile that
#'         can be considered outliers on each side of the median for the
#'         \code{bicor} statistic.
#' @param robustY A \code{logical(1)} indicating whether to use robust calculation
#'         for \code{y} for the \code{bicor} statistic. \code{FALSE} is recommended
#'         if \code{y} is a binary variable.
#' @param adjMethod A \code{character(1)} specifying the method for adjusting
#'         p-values, Potential values include \code{fdr}, \code{holm},
#'         \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH},
#'         \code{BY} and \code{none}.
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{data.frame}.
#' @param file A \code{character(1)} giving the file name (.txt) for the saved
#'         \code{data.frame}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{data.frame} giving correlation statistics for each module-trait
#'         pair.
#'
#' @seealso \itemize{
#'         \item \code{\link{getModules()}} to build a comethylation network and
#'                 identify modules of comethylated regions.
#'         \item \code{\link{getCor()}} to calculate correlation coefficients.
#'         \item \code{\link{getDendro()}} and \code{\link{plotDendro()}} to
#'                 generate and visualize dendrograms.
#'         \item \code{\link{plotMEtraitCor()}} to visualize ME-trait correlations.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Test Correlations between Module Eigennodes and Sample Traits
#' MEs <- modules$MEs
#' MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor",
#'                             file = "ME_Trait_Correlation_Stats.txt")
#'
#' # Examine Correlations between Sample Traits
#' traitDendro <- getCor(MEs, y = colData, corType = "bicor",
#'                       robustY = FALSE) %>%
#'         getDendro(transpose = TRUE)
#' plotDendro(traitDendro, labelSize = 3.5, expandY = c(0.65,0.08),
#'            file = "Trait_Dendrogram.pdf")
#'
#' # Visualize Correlations between Module Eigennodes and Sample Traits
#' moduleDendro <- getDendro(MEs, distance = "bicor")
#' plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order,
#'                traitOrder = traitDendro$order,
#'                file = "ME_Trait_Correlation_Heatmap.pdf")
#' plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order,
#'                traitOrder = traitDendro$order, sigOnly = TRUE,
#'                star.size = 11, star.nudge_y = -0.27,
#'                legend.position = c(1.14, 0.745),
#'                colColorMargins = c(-1,5.1,0.5,10.47),
#'                file = "Sig_ME_Trait_Correlation_Heatmap.pdf",
#'                width = 7, height = 3.5)
#' }
#'
#' @export
#'
#' @import WGCNA
#' @importFrom magrittr %>%
#' @importFrom stringr str_remove_all

getMEtraitCor <- function(MEs, colData, corType = c("bicor", "pearson"),
                          maxPOutliers = 0.1, robustY = FALSE,
                          adjMethod = c("fdr", "holm", "hochberg", "hommel",
                                        "bonferroni", "BH", "BY", "none"),
                          save = TRUE, file = "ME_Trait_Correlation_Stats.txt",
                          verbose = TRUE){
        corType <- match.arg(corType)
        adjMethod <- match.arg(adjMethod)
        if(verbose){
                message("[getMEtraitCor] Testing associations between module eigennodes and sample traits using ",
                        corType, " correlation")
        }
        colData <- colData[rownames(MEs),]
        if(corType == "bicor"){
                cor <- bicorAndPvalue(x = MEs, y = colData,
                                      maxPOutliers = maxPOutliers,
                                      robustY = robustY,
                                      pearsonFallback = "none")
        } else {
                if(corType == "pearson"){
                        cor <- corAndPvalue(x = MEs, y = colData)
                } else {
                        stop("[getMEtraitCor] corType must be either bicor or pearson")
                }
        }
        stats <- list.rbind(cor) %>% as.data.frame()
        stats$module <- rownames(cor$p) %>% str_remove_all(pattern = "ME") %>%
                rep(length(cor)) %>% factor(levels = unique(.))
        stats$stat <- names(cor) %>% rep(each = nrow(cor$p)) %>%
                factor(levels = unique(.))
        stats <- reshape2::melt(stats, id.vars = c("module", "stat"),
                                variable.name = "trait") %>%
                reshape2::dcast(formula = module + trait ~ stat,
                                value.var = "value")
        if(verbose){
                message("[getMEtraitCor] Adjusting p-values using the ",
                        adjMethod, " method")
        }
        stats$adj_p <- p.adjust(stats$p, method = adjMethod)
        if(corType == "bicor"){
                stats <- stats[,c("module", "trait", "nObs", "bicor", "Z", "t",
                                  "p", "adj_p")]
        } else {
                stats <- stats[,c("module", "trait", "nObs", "cor", "Z", "t",
                                  "p", "adj_p")]
        }
        if(save){
                if(verbose){
                        message("[getMEtraitCor] Saving file as ", file)
                }
                write.table(stats, file = file, quote = FALSE, sep = "\t",
                            row.names = FALSE)
        }
        return(stats)
}

#' Plot a Heatmap of Correlations Between Module Eigennodes and Traits
#'
#' \code{plotMEtraitCor()} takes a \code{data.frame} of correlation statistics
#' for module eigennodes and traits from \code{\link{getMEtraitCor()}}, plots it
#' as a heatmap, and saves it as a .pdf.
#'
#' \code{plotMEtraitCor()} is designed to be used in combination with
#' \code{\link{getMEtraitCor()}}. Stars are used to annotate significant
#' correlations, as defined by the \code{adj_p} threshold. The heatmap can
#' optionally be filtered to include only modules and traits with significant
#' correlations. A \code{ggplot} object is produced and can be edited outside of
#' this function if desired.
#'
#' @param MEtraitCor A \code{data.frame} of correlation statistics for module-trait
#'         pairs, typically generated by \code{\link{getMEtraitCore()}}.
#' @param moduleOrder A \code{numeric} specifying the order of modules in the
#'         heatmap. This must be the same length as the number of unique modules
#'         in \code{MEtraitCor}.
#' @param traitOrder A \code{numeric} specifying the order of traits in the
#'         heatmap. This must be the same length as the number of unique traits
#'         in \code{MEtraitCor}.
#' @param sigOnly A \code{logical(1)} indicating whether to plot only modules and
#'         traits with significant correlations.
#' @param adj_p A \code{numeric(1)} giving the threshold for the adjusted p-value
#'         to determine a significant correlation.
#' @param star.size A \code{numeric(1)} specifying the size of the star for
#'         annotating significant correlations.
#' @param star.nudge_y A \code{numeric(1)} giving the amount to adjust the position
#'         of the star on the y-axis.
#' @param colors A \code{character} giving a vector of colors to use for the
#'         gradient on the heatmap. The default uses \code{\link[WGCNA]{blueWhiteRed()}}
#'         to generate these colors.
#' @param limit A \code{numeric(1)} giving the maximum value (symmetric) for the
#'         heatmap color scale.
#' @param axis.text.size A \code{numeric(1)} specifying the size of the text for
#'         the y-axis.
#' @param legend.position A \code{numeric(2)} with the position of the legend,
#'         as x-axis, y-axis. May also be a \code{character(1)} indicating "none",
#'         "left", "right", "bottom", or "top".
#' @param legend.text.size A \code{numeric(1)} specifying the size of the text for
#'         the legend axis.
#' @param legend.title.size A \code{numeric(1)} giving the size of the text for
#'         the legend title.
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
#'         \item \code{\link{getMEtraitCor()}} to calculate ME-trait correlations.
#'         \item \code{\link{getCor()}} to calculate correlation coefficients.
#'         \item \code{\link{getDendro()}} and \code{\link{plotDendro()}} to
#'                 generate and visualize dendrograms.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Test Correlations between Module Eigennodes and Sample Traits
#' MEs <- modules$MEs
#' MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor",
#'                             file = "ME_Trait_Correlation_Stats.txt")
#'
#' # Examine Correlations between Sample Traits
#' traitDendro <- getCor(MEs, y = colData, corType = "bicor",
#'                       robustY = FALSE) %>%
#'         getDendro(transpose = TRUE)
#' plotDendro(traitDendro, labelSize = 3.5, expandY = c(0.65,0.08),
#'            file = "Trait_Dendrogram.pdf")
#'
#' # Visualize Correlations between Module Eigennodes and Sample Traits
#' moduleDendro <- getDendro(MEs, distance = "bicor")
#' plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order,
#'                traitOrder = traitDendro$order,
#'                file = "ME_Trait_Correlation_Heatmap.pdf")
#' plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order,
#'                traitOrder = traitDendro$order, sigOnly = TRUE,
#'                star.size = 11, star.nudge_y = -0.27,
#'                legend.position = c(1.14, 0.745),
#'                colColorMargins = c(-1,5.1,0.5,10.47),
#'                file = "Sig_ME_Trait_Correlation_Heatmap.pdf",
#'                width = 7, height = 3.5)
#' }
#'
#' @export
#'
#' @import ggplot2
#' @import cowplot
#' @import WGCNA
#' @importFrom magrittr %>%

plotMEtraitCor <- function(MEtraitCor,
                           moduleOrder = 1:length(unique(MEtraitCor$module)),
                           traitOrder = 1:length(unique(MEtraitCor$trait)),
                           sigOnly = FALSE, adj_p = 0.05, star.size = 8,
                           star.nudge_y = -0.38,
                           colors = blueWhiteRed(100, gamma = 0.9), limit = NULL,
                           axis.text.size = 12, legend.position = c(1.08, 0.915),
                           legend.text.size = 12, legend.title.size = 16,
                           colColorMargins = c(-0.7,4.21,1.2,11.07), save = TRUE,
                           file = "ME_Trait_Correlation_Heatmap.pdf", width = 11,
                           height = 9.5, verbose = TRUE){
        if(verbose){
                message("[plotMEtraitCor] Plotting ME trait correlation heatmap")
        }
        MEtraitCor$module <- factor(MEtraitCor$module,
                                    levels = levels(MEtraitCor$module)[moduleOrder])
        MEtraitCor$trait <- factor(MEtraitCor$trait,
                                   levels = levels(MEtraitCor$trait)[rev(traitOrder)])
        MEtraitCor$Significant <- (MEtraitCor$adj_p < adj_p & !is.na(MEtraitCor$adj_p)) %>%
                factor(levels = c("TRUE", "FALSE"))
        if(sigOnly){
                sigModules <- MEtraitCor$module[MEtraitCor$Significant == "TRUE"] %>%
                        unique() %>% as.character()
                sigTraits <- MEtraitCor$trait[MEtraitCor$Significant == "TRUE"] %>%
                        unique() %>% as.character()
                MEtraitCor <- subset(MEtraitCor, module %in% sigModules &
                                             trait %in% sigTraits)
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
        corData <- MEtraitCor[[corType]]
        if(is.null(limit)){
                limit <- max(abs(corData))
        }
        heatmap <- ggplot(data = MEtraitCor) +
                geom_tile(aes(x = module, y = trait, color = corData,
                              fill = corData)) +
                geom_text(aes(x = module, y = trait, alpha = Significant),
                          label = "*", color = "black",
                          size = star.size, nudge_y = star.nudge_y) +
                scale_fill_gradientn(str_to_title(corType), colors = colors,
                                     limits = c(-limit, limit),
                                     aesthetics = c("color", "fill")) +
                scale_x_discrete(expand = expansion(mult = 0.01)) +
                scale_y_discrete(expand = expansion(mult = 0.01)) +
                scale_alpha_manual(breaks = c("TRUE", "FALSE"),
                                   values = c("TRUE" = 1, "FALSE" = 0),
                                   guide = FALSE) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_text(size = axis.text.size,
                                                 color = "black"),
                      axis.ticks.x = element_blank(),
                      axis.ticks.y = element_line(size = 0.8, color = "black"),
                      axis.title = element_blank(),
                      legend.background = element_blank(),
                      legend.position = legend.position,
                      legend.text = element_text(size = legend.text.size),
                      legend.title = element_text(size = legend.title.size),
                      panel.background = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(),
                      plot.background = element_blank(),
                      plot.margin = unit(c(1,6,1,1), "lines"))
        colColors <- ggplot(data = data.frame(x = 1:length(levels(MEtraitCor$module)),
                                              y = 0,
                                              color = levels(MEtraitCor$module))) +
                geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                scale_fill_identity(aesthetics = c("color", "fill")) +
                theme_void() +
                theme(legend.position = "none",
                      plot.margin = unit(colColorMargins, "lines"))
        gg <- plot_grid(heatmap, colColors, nrow = 2, rel_heights = c(1, 0.045))
        if(save){
                if(verbose){
                        message("[plotMEtraitCor] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width,
                       height = height, units = "in")
        }
        return(gg)
}

#' Visualize a Module Eigennode - Trait Correlation as a Dot Plot
#'
#' \code{plotMEtraitDot()} takes a \code{vector} of module eigennode values and
#' a \code{vector} of categorical sample trait values, generates a dot plot, and
#' then saves it as a .pdf. \code{ME} and \code{trait} must be in the same order.
#'
#' \code{NA} values in the trait are removed if present, along with corresponding
#' \code{ME} values. Data is summarized like a box plot (median, Q1, Q3) by
#' default, but can also be summarized with other methods. See
#' \code{\link[ggplot2]{stat_summary()}} for more details. A \code{ggplot} object
#' is produced and can be edited outside of this function if desired.
#'
#' @param ME A \code{numeric} of module eigennode values. \code{ME} must be in the
#'         same order as \code{trait}.
#' @param trait A \code{numeric} of categorical sample trait values.
#' @param traitCode A named \code{numeric vector} matching each trait level to a
#'         numeric value. Example: c("Control" = 0, "Treatment" = 1).
#' @param colors A named \code{character vector} matching each trait level to a
#'         color. Example: c("Control" = "blue", "Treatment" = "red").
#' @param fun.data A \code{character(1)} specifying the summary function to use.
#'         Potential values include \code{median_hilow}, \code{mean_cl_boot},
#'         \code{mean_cl_normal}, and \code{mean_sdl}.
#' @param fun.args A \code{list} giving additional arguments to the summary
#'         function.
#' @param binwidth A \code{numeric(1)} specifying the maximum bin width.
#' @param stackratio A \code{numeric(1)} indicating how far apart to stack the dots.
#' @param dotsize A \code{numeric(1)} giving the size of the dots relative to
#'         \code{binwidth}.
#' @param ylim A \code{numeric(2)} specifying the limits of the y-axis.
#' @param nBreaks A \code{numeric(1)} giving the number of breaks on the y-axis.
#' @param axis.title.size A \code{numeric(1)} indicating the size of the title
#'         text for both axes.
#' @param axis.text.size A \code{numeric(1)} specifying the size of the text for
#'         both axes.
#' @param xlab A \code{character(1)} giving the x-axis title.
#' @param ylab A \code{character(1)} giving the y-axis title.
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
#'         \item \code{\link{getMEtraitCor()}} and \code{\link{plotMEtraitCor()}}
#'                 to calculate and visualize all ME-trait correlations.
#'         \item \code{\link{plotMEtraitScatter()}} and \code{\link{plotMethTrait()}}
#'                 for other methods to visualize a single ME-trait correlation.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Test Correlations between Module Eigennodes and Sample Traits
#' MEs <- modules$MEs
#' MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor",
#'                             file = "ME_Trait_Correlation_Stats.txt")
#' plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order,
#'                traitOrder = traitDendro$order,
#'                file = "ME_Trait_Correlation_Heatmap.pdf")
#'
#' # Explore Individual ME-Trait Correlations
#' plotMEtraitDot(MEs$bisque4, trait = colData$Diagnosis_ASD,
#'                traitCode = c("TD" = 0, "ASD" = 1),
#'                colors = c("TD" = "#3366CC", "ASD" = "#FF3366"),
#'                ylim = c(-0.2,0.2), xlab = "Diagnosis",
#'                ylab = "Bisque 4 Module Eigennode",
#'                file = "bisque4_ME_Diagnosis_Dotplot.pdf")
#' plotMEtraitScatter(MEs$paleturquoise, trait = colData$Gran,
#'                    ylim = c(-0.15,0.15), xlab = "Granulocytes",
#'                    ylab = "Pale Turquoise Module Eigennode",
#'                    file = "paleturquoise_ME_Granulocytes_Scatterplot.pdf")
#' regions <- modules$regions
#' plotMethTrait("bisque4", regions = regions, meth = meth,
#'               trait = colData$Diagnosis_ASD,
#'               traitCode = c("TD" = 0, "ASD" = 1),
#'               traitColors = c("TD" = "#3366CC", "ASD" = "#FF3366"),
#'               trait.legend.title = "Diagnosis",
#'               file = "bisque4_Module_Methylation_Diagnosis_Heatmap.pdf")
#' }
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%

plotMEtraitDot <- function(ME, trait, traitCode = NULL, colors = NULL,
                           fun.data = c("median_hilow", "mean_cl_boot",
                                        "mean_cl_normal", "mean_sdl"),
                           fun.args = list(conf.int = 0.5), binwidth = 0.01,
                           stackratio = 1.4, dotsize = 0.85, ylim = NULL,
                           nBreaks = 4, axis.title.size = 20, axis.text.size = 16,
                           xlab = "Trait", ylab = "Module Eigennode", save = TRUE,
                           file = "ME_Trait_Dotplot.pdf", width = 6, height = 6,
                           verbose = TRUE){
        if(verbose){
                message("[plotMEtraitDot] Plotting module eigennode by categorical trait")
        }
        if(sum(is.na(trait)) > 1){
                message("[plotMEtraitDot] Removing NA values")
                ME <- ME[!is.na(trait)]
                trait <- trait[!is.na(trait)]
        }
        if(!is.null(traitCode)){
                trait <- names(traitCode)[match(trait, traitCode)] %>%
                        factor(levels = names(traitCode))
        } else {
                trait <- as.factor(trait)
        }
        fun.data <- match.arg(fun.data)
        dotplot <- ggplot() +
                stat_summary(aes(x = trait, y = ME, group = trait),
                             fun.data = fun.data, geom = "crossbar",
                             color = "black", size = 0.5, fun.args = fun.args) +
                geom_dotplot(aes(x = trait, y = ME, fill = trait, color = trait),
                             binwidth = binwidth, binaxis = "y",
                             stackdir = "center", position = "dodge",
                             stackratio = stackratio, dotsize = dotsize) +
                coord_cartesian(ylim = ylim) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                theme_bw(base_size = 25) +
                theme(legend.position = "none", panel.grid.major = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25),
                      axis.ticks = element_line(size = 1.25),
                      panel.grid.minor = element_blank(),
                      strip.background = element_blank(),
                      axis.text = element_text(color = "black", size = axis.text.size),
                      axis.title = element_text(size = axis.title.size),
                      plot.margin = unit(c(1,1,1,1), "lines")) +
                xlab(xlab) +
                ylab(ylab)
        if(!is.null(colors)){
                dotplot <- dotplot +
                        scale_color_manual(breaks = names(colors),
                                           values = colors,
                                           aesthetics = c("color", "fill"))
        }
        if(verbose){
                message("[plotMEtraitDot] Saving file as ", file)
        }
        ggsave(file, plot = dotplot, dpi = 600, width = width, height = height,
               units = "in")
}

#' Visualize a Module Eigennode - Trait Correlation as a Scatter Plot
#'
#' \code{plotMEtraitScatter()} takes a \code{vector} of module eigennode values
#' and a \code{vector} of continuous sample trait values, generates a scatter plot,
#' and then saves it as a .pdf. \code{ME} and \code{trait} must be in the same
#' order.
#'
#' The values in \code{ME} and \code{trait} are plotted as points along with a
#' smoothed line with a shaded 95% confidence interval. The smoothed line is fit
#' using robust regression as implemented by \code{\link[MASS]{rlm}}. A
#' \code{ggplot} object is produced and can be edited outside of this function
#' if desired.
#'
#' @param ME A \code{numeric} of module eigennode values. \code{ME} must be in the
#'         same order as \code{trait}.
#' @param trait A \code{numeric} of continuous sample trait values.
#' @param color A \code{character(1)} giving the color of the points.
#' @param xlim A \code{numeric(2)} specifying the limits of the x-axis.
#' @param ylim A \code{numeric(2)} specifying the limits of the y-axis.
#' @param nBreaks A \code{numeric(1)} giving the number of breaks for both axes.
#' @param point.size A \code{numeric(1)} indicating the size of the points.
#' @param axis.title.size A \code{numeric(1)} indicating the size of the title
#'         text for both axes.
#' @param axis.text.size A \code{numeric(1)} specifying the size of the text for
#'         both axes.
#' @param xlab A \code{character(1)} giving the x-axis title.
#' @param ylab A \code{character(1)} giving the y-axis title.
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
#'         \item \code{\link{getMEtraitCor()}} and \code{\link{plotMEtraitCor()}}
#'                 to calculate and visualize all ME-trait correlations.
#'         \item \code{\link{plotMEtraitDot()}} and \code{\link{plotMethTrait()}}
#'                 for other methods to visualize a single ME-trait correlation.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Test Correlations between Module Eigennodes and Sample Traits
#' MEs <- modules$MEs
#' MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor",
#'                             file = "ME_Trait_Correlation_Stats.txt")
#' plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order,
#'                traitOrder = traitDendro$order,
#'                file = "ME_Trait_Correlation_Heatmap.pdf")
#'
#' # Explore Individual ME-Trait Correlations
#' plotMEtraitDot(MEs$bisque4, trait = colData$Diagnosis_ASD,
#'                traitCode = c("TD" = 0, "ASD" = 1),
#'                colors = c("TD" = "#3366CC", "ASD" = "#FF3366"),
#'                ylim = c(-0.2,0.2), xlab = "Diagnosis",
#'                ylab = "Bisque 4 Module Eigennode",
#'                file = "bisque4_ME_Diagnosis_Dotplot.pdf")
#' plotMEtraitScatter(MEs$paleturquoise, trait = colData$Gran,
#'                    ylim = c(-0.15,0.15), xlab = "Granulocytes",
#'                    ylab = "Pale Turquoise Module Eigennode",
#'                    file = "paleturquoise_ME_Granulocytes_Scatterplot.pdf")
#' regions <- modules$regions
#' plotMethTrait("bisque4", regions = regions, meth = meth,
#'               trait = colData$Diagnosis_ASD,
#'               traitCode = c("TD" = 0, "ASD" = 1),
#'               traitColors = c("TD" = "#3366CC", "ASD" = "#FF3366"),
#'               trait.legend.title = "Diagnosis",
#'               file = "bisque4_Module_Methylation_Diagnosis_Heatmap.pdf")
#' }
#'
#' @export
#'
#' @import ggplot2

plotMEtraitScatter <- function(ME, trait, color = "#132B43", xlim = NULL,
                               ylim = NULL, nBreaks = 4, point.size = 2.5,
                               axis.title.size = 20, axis.text.size = 16,
                               xlab = "Trait", ylab = "Module Eigennode",
                               save = TRUE, file = "ME_Trait_Scatterplot.pdf",
                               width = 6, height = 6, verbose = TRUE){
        if(verbose){
                message("[plotMEtraitScatter] Plotting module eigennode by continuous trait")
        }
        scatterplot <- ggplot() +
                geom_smooth(aes(x = trait, y = ME), method = MASS::rlm,
                            formula = y ~ x, color = "#56B1F7", fill = "#336A98") +
                geom_point(aes(x = trait, y = ME), color = color,
                           size = point.size) +
                coord_cartesian(xlim = xlim, ylim = ylim) +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                theme_bw(base_size = 25) +
                theme(legend.position = "none", panel.grid.major = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25),
                      axis.ticks = element_line(size = 1.25),
                      panel.grid.minor = element_blank(),
                      strip.background = element_blank(),
                      axis.text = element_text(color = "black",
                                               size = axis.text.size),
                      axis.title = element_text(size = axis.title.size),
                      plot.margin = unit(c(1,1,1,1), "lines")) +
                xlab(xlab) +
                ylab(ylab)
        if(verbose){
                message("[plotMEtraitScatter] Saving file as ", file)
        }
        ggsave(file, plot = scatterplot, dpi = 600, width = width,
               height = height, units = "in")
}

#' Plot Module Methylation Values By a Sample Trait
#'
#' \code{plotMethTrait()} extracts the methylation values for regions in a given
#' module, plots it against a sample trait in a heatmap, and then saves it as a
#' .pdf. \code{trait} must be in the same order as the columns in \code{meth}.
#'
#' \code{NA} values in the trait are removed if present, along with corresponding
#' \code{ME} values. If \code{discrete} is not provided \code{plotMethTrait()}
#' will guess if the trait is discrete (<= 5 unique values) and plot the trait
#' color as a discrete scale rather than a continuous one. Samples are ordered by
#' trait value in ascending order. Methylation values are plotted relative to the
#' mean methylation in that region.
#'
#' @param module A \code{character(1)} giving the name of the module to plot.
#' @param regions A \code{data.frame} of regions with module assignments, typically
#'         obtained from \code{\link{getModules()}}.
#' @param meth A \code{numeric matrix}, where each row is a region and each
#'         column is a sample. This is typically obtained from
#'         \code{\link{getRegionMeth()}}.
#' @param trait A \code{numeric} of sample trait values.
#' @param discrete A \code{logical(1)} identifying \code{trait} as a discrete
#'         variable or not. If null, \code{plotMethTrait()} will guess if the
#'         trait is discrete (<= 5 unique values).
#' @param traitCode A named \code{numeric vector} matching each trait level to a
#'         numeric value. Example: c("Control" = 0, "Treatment" = 1).
#' @param traitColors A named \code{character vector} matching each trait level
#'         to a color. Example: c("Control" = "blue", "Treatment" = "red").
#' @param heatmapColors A \code{character} giving a vector of colors to use for
#'         the gradient on the heatmap. The default uses
#'         \code{\link[WGCNA]{blueWhiteRed()}} to generate these colors.
#' @param limit A \code{numeric(1)} giving the maximum value (symmetric) for the
#'         heatmap color scale.
#' @param expandY A \code{numeric} specifying the multiplicative range expansion
#'         factors. Can be given as the symmetric lower and upper limit expansions
#'         or separately as a \code{vector} of length 2.
#' @param axis.text.size A \code{numeric(1)} giving the size of the text on the
#'         y-axis.
#' @param heatmap.legend.position A \code{numeric(2)} with the position of the
#'         heatmap legend, as x-axis, y-axis. May also be a  \code{character(1)}
#'         indicating "none", "left", "right", "bottom", or "top".
#' @param trait.legend.position A \code{numeric(2)} with the position of the
#'         color bar legend, as x-axis, y-axis. May also be a \code{character(1)}
#'         indicating "none", "left", "right", "bottom", or "top".
#' @param heatmap.legend.title A \code{character(1)} giving the title for the
#'         heatmap legend.
#' @param trait.legend.title A \code{character(1)} giving the title for the
#'         color bar legend.
#' @param legend.text.size A \code{numeric(1)} indicating the size the text in
#'         both legends.
#' @param legend.title.size A \code{numeric(1)} specifying the size of the text
#'         for both legend titles.
#' @param heatmapMargins A \code{numeric(4)} giving the width of the margins for
#'         the heatmap as top, right, bottom, and left.
#' @param traitMargins A \code{numeric(4)} giving the width of the margins for
#'         the color bar as top, right, bottom, and left.
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
#'         \item \code{\link{getMEtraitCor()}} and \code{\link{plotMEtraitCor()}}
#'                 to calculate and visualize all ME-trait correlations.
#'         \item \code{\link{plotMEtraitDot()}} and \code{\link{plotMEtraitScatter()}}
#'                 for other methods to visualize a single ME-trait correlation.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Test Correlations between Module Eigennodes and Sample Traits
#' MEs <- modules$MEs
#' MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor",
#'                             file = "ME_Trait_Correlation_Stats.txt")
#' plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order,
#'                traitOrder = traitDendro$order,
#'                file = "ME_Trait_Correlation_Heatmap.pdf")
#'
#' # Explore Individual ME-Trait Correlations
#' plotMEtraitDot(MEs$bisque4, trait = colData$Diagnosis_ASD,
#'                traitCode = c("TD" = 0, "ASD" = 1),
#'                colors = c("TD" = "#3366CC", "ASD" = "#FF3366"),
#'                ylim = c(-0.2,0.2), xlab = "Diagnosis",
#'                ylab = "Bisque 4 Module Eigennode",
#'                file = "bisque4_ME_Diagnosis_Dotplot.pdf")
#' plotMEtraitScatter(MEs$paleturquoise, trait = colData$Gran,
#'                    ylim = c(-0.15,0.15), xlab = "Granulocytes",
#'                    ylab = "Pale Turquoise Module Eigennode",
#'                    file = "paleturquoise_ME_Granulocytes_Scatterplot.pdf")
#' regions <- modules$regions
#' plotMethTrait("bisque4", regions = regions, meth = meth,
#'               trait = colData$Diagnosis_ASD,
#'               traitCode = c("TD" = 0, "ASD" = 1),
#'               traitColors = c("TD" = "#3366CC", "ASD" = "#FF3366"),
#'               trait.legend.title = "Diagnosis",
#'               file = "bisque4_Module_Methylation_Diagnosis_Heatmap.pdf")
#' }
#'
#' @export
#'
#' @import ggplot2
#' @import cowplot
#' @import WGCNA
#'
#' @importFrom magrittr %>%

plotMethTrait <- function(module, regions, meth, trait, discrete = NULL,
                          traitCode = NULL, traitColors = NULL,
                          heatmapColors = blueWhiteRed(100, gamma = 0.3),
                          limit = NULL, expandY = 0.05, axis.text.size = 11,
                          heatmap.legend.position = c(1.1,0.743),
                          trait.legend.position = c(1.017,4.39),
                          heatmap.legend.title = "Relative\nMethylation (%)",
                          trait.legend.title = "Trait", legend.text.size = 11,
                          legend.title.size = 14, heatmapMargins = c(1,8,0,1),
                          traitMargins = c(0,6,1,5.15), save = TRUE,
                          file = "Module_Methylation_Trait_Heatmap.pdf",
                          width = 11, height = 4, verbose = TRUE){
        if(verbose){
                message("[plotMethTrait] Plotting ", module,
                        " module region methylation by ", trait.legend.title)
        }
        if(sum(is.na(trait)) > 1){
                message("[plotMethTrait] Removing NA values")
                meth <- meth[,!is.na(trait)]
                trait <- trait[!is.na(trait)]
        }
        RegionIDs <- rev(regions$RegionID[regions$module == module])
        if(is.null(discrete)){
                discrete <- length(unique(trait)) <= 5
        }
        if(discrete){
                if(!is.null(traitCode)){
                        trait <- names(traitCode)[match(trait, traitCode)] %>%
                                factor(levels = names(traitCode))
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
                scale_fill_gradientn(heatmap.legend.title, colors = heatmapColors,
                                     limits = c(-limit,limit),
                                     aesthetics = c("color", "fill")) +
                scale_y_discrete(expand = expansion(mult = expandY)) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_text(size = axis.text.size,
                                                 color = "black"),
                      axis.ticks.x = element_blank(),
                      axis.ticks.y = element_line(size = 0.8, color = "black"),
                      axis.title = element_blank(),
                      legend.background = element_blank(),
                      legend.position = heatmap.legend.position,
                      legend.text = element_text(size = legend.text.size),
                      legend.title = element_text(size = legend.title.size),
                      panel.background = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(),
                      plot.background = element_blank(),
                      plot.margin = unit(heatmapMargins, "lines"))
        colColors <- ggplot(data = data.frame(x = 1:length(trait), y = 0,
                                              color = sort(trait))) +
                geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                theme_void() +
                theme(legend.position = trait.legend.position,
                      legend.text = element_text(size = legend.text.size),
                      legend.title = element_text(size = legend.title.size),
                      plot.margin = unit(traitMargins, "lines"))
        if(discrete){
                if(!is.null(traitColors)){
                        colColors <- colColors +
                                scale_color_manual(trait.legend.title,
                                                   breaks = names(traitColors),
                                                   values = traitColors,
                                                   aesthetics = c("color", "fill"))
                } else {
                        colColors <- colColors +
                                scale_color_discrete(trait.legend.title,
                                                     aesthetics = c("color", "fill"))
                }
        } else {
                colColors <- colColors +
                        scale_color_continuous(trait.legend.title,
                                               aesthetics = c("color", "fill"))
        }
        gg <- plot_grid(heatmap, colColors, nrow = 2, rel_heights = c(1,0.15))
        if(save){
                if(verbose){
                        message("[plotMethTrait] Saving file as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width,
                       height = height, units = "in")
        }
}

