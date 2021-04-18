#' Read Bismark CpG reports
#'
#' \code{getCpGs()} reads individual sample Bismark CpG reports into a single
#' \code{\link[bsseq:`BSseq-class`]{BSseq}} object and then saves it as a .rds
#' file.
#'
#' This \code{\link[bsseq:`BSseq-class`]{BSseq}} object still needs to be
#'         filtered for coverage at individual CpGs. More information on these
#'         arguments is given in the documentation for
#'         \code{\link[bsseq]{read.bismark()}}.
#'
#' @param colData A \code{data.frame} whose row names specify CpG reports to
#'         load into the \code{\link[bsseq:`BSseq-class`]{BSseq}} object and
#'         whose columns are sample traits with numeric values.
#' @param path A \code{character} giving the path(s) to the CpG reports.
#' @param pattern A regular expression used to filter for CpG reports.
#' @param sameLoci A \code{logical(1)} indicating whether CpG reports contain
#'         the same set of methylation loci. This is the case if the files are
#'         genome wide cytosine reports aligned to the same reference genome.
#'         The default \code{sameLoci = TRUE} speeds up \code{getCpGs()} by
#'         only having to parse each CpG report once.
#' @param chroms A \code{character} giving the chromosomes to include in the
#'         \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} instance
#'         providing the parallel back-end to use during evaluation.
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#' @param file A \code{character(1)} giving the file name (.rds) for the saved
#'         \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#'
#' @seealso \itemize{
#'         \item \code{\link{getCpGtotals()}} and \code{\link{plotCpGtotals()}}
#'                 for help with deciding coverage cutoffs.
#'         \item \code{\link{filterCpGs()}} to filter the
#'                 \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#'         \item \code{\link[bsseq]{read.bismark()}} for more details on the
#'                 arguments and the underlying functions.
#' }
#'
#' @examples \dontrun{
#'
#' # Read Bismark CpG Reports
#' colData <- read.xlsx("sample_info.xlsx", rowNames = TRUE)
#' bs <- getCpGs(colData, file = "Unfiltered_BSseq.rds")
#'
#' # Examine CpG Totals at Different Cutoffs
#' CpGtotals <- getCpGtotals(bs, file = "CpG_Totals.txt")
#' plotCpGtotals(CpGtotals, file = "CpG_Totals.pdf")
#'
#' # Filter BSseq Object
#' bs <- filterCpGs(bs, cov = 2, perSample = 0.75, file = "Filtered_BSseq.rds")
#' }
#'
#' @export
#'
#' @import bsseq
#' @importFrom magrittr %>%

getCpGs <- function(colData, path = getwd(), pattern = "*CpG_report.txt.gz",
                    sameLoci = TRUE,
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY",
                               "chrM"),
                    BPPARAM = BiocParallel::MulticoreParam(10), save = TRUE,
                    file = "Unfiltered_BSseq.rds", verbose = TRUE){
        if(verbose){
                message("[getCpGs] Loading CpG-level data")
        }
        files <- list.files(path, pattern = pattern) %>%
                .[pmatch(rownames(colData), table = .)]
        if(sameLoci){
                loci <- bsseq:::.readBismarkAsFWGRanges(files[1])
        } else {
                loci <- NULL
        }
        bs <- read.bismark(files, loci = loci, colData = colData,
                           BPPARAM = BPPARAM, verbose = verbose) %>%
                chrSelectBSseq(seqnames = chroms)
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

#' Get the Total CpGs at Different Coverage Cutoffs
#'
#' \code{getCpGtotals()} calculates the total number and percent of CpGs
#'         remaining in a \code{\link[bsseq:`BSseq-class`]{BSseq}} object after
#'         filtering at different \code{cov} and \code{perSample} cutoffs and
#'         then saves it as a tab-delimited text file.
#'
#' The purpose of this function is to help determine cutoffs to maximize the
#'         number of CpGs with sufficient data after filtering. Typically,
#'         the number of CpGs covered in 100% of samples decreases as the
#'         sample size increases, especially with low-coverage datasets.
#'
#' @param bs A \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#' @param cov A \code{numeric} specifying the minimum number of reads
#'         overlapping a CpG for it to be included in the total.
#' @param perSample A \code{numeric} specifying the minimum percent of samples
#'         with \code{cov} reads at a CpG for it to be included in the total.
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{data.frame}.
#' @param file A \code{character(1)} giving the file name (.txt) for the saved
#'         \code{data.frame}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{data.frame} giving the number of CpGs (in millions) and
#'         the percent of total CpGs at all combinations of \code{cov} and
#'         \code{perSample}. The number of samples corresponding to
#'         \code{perSample} is also given.
#'
#' @seealso \itemize{
#'         \item \code{\link{getCpGs()}} to generate the
#'                 \code{\link[bsseq:`BSseq-class`]{BSseq}} object from
#'                 individual Bismark CpG reports.
#'         \item \code{\link{plotCpGtotals()}} to visualize the CpG totals.
#'         \item \code{\link{filterCpGs()}} to filter the
#'                 \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#' }
#'
#' @examples \dontrun{
#'
#' # Read Bismark CpG Reports
#' colData <- read.xlsx("sample_info.xlsx", rowNames = TRUE)
#' bs <- getCpGs(colData, file = "Unfiltered_BSseq.rds")
#'
#' # Examine CpG Totals at Different Cutoffs
#' CpGtotals <- getCpGtotals(bs, file = "CpG_Totals.txt")
#' plotCpGtotals(CpGtotals, file = "CpG_Totals.pdf")
#'
#' # Filter BSseq Object
#' bs <- filterCpGs(bs, cov = 2, perSample = 0.75, file = "Filtered_BSseq.rds")
#' }
#'
#' @export
#'
#' @import bsseq
#' @importFrom magrittr %>%

getCpGtotals <- function(bs, cov = seq(0,10,1), perSample = seq(0.5,1,0.05),
                         save = TRUE, file = "CpG_Totals.txt", verbose = TRUE){
        if(verbose){
                message("[getCpGtotals] Calculating CpG totals at specified ",
                        "cov and perSample cutoffs")
        }
        nSample <- (perSample * ncol(bs)) %>% ceiling()
        bsCov <- getCoverage(bs)
        nCpGs <- sapply(cov, FUN = .nCpGsByCov, bsCov = bsCov,
                        nSample = nSample)
        nCpGs <- as.integer(nCpGs)
        perCpGs <- (nCpGs / length(bs))
        CpGtotals <- data.frame(cov = rep(cov, each = length(perSample)),
                                perSample = rep(perSample, times = length(cov)),
                                nSample = rep(nSample, times = length(cov)),
                                nCpGs_M = nCpGs / 10^6, perCpGs = perCpGs)
        if(save){
                if(verbose){
                        message("[getCpGtotals] Saving file as ", file)
                }
                write.table(CpGtotals, file = file, sep = "\t",
                            row.names = FALSE)
        }
        return(CpGtotals)
}

#' Visualize CpG Totals
#'
#' \code{plotCpGtotals()} plots the number of CpGs remaining after filtering by
#'         different combinations of \code{cov} and \code{perSample} in a line
#'         plot and then saves it as a pdf.
#'
#' \code{plotCpGtotals()} is designed to be used in combination with
#'         \code{\link{getCpGtotals()}}. A \code{ggplot} object is produced and
#'         can be further edited outside of this function if desired.
#'
#' @param CpGtotals A \code{data.frame}, output from \code{\link{getCpGs()}}.
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
#' @seealso \itemize{
#'         \item \code{\link{getCpGs()}} to generate the
#'                 \code{\link[bsseq:`BSseq-class`]{BSseq}} object from
#'                 individual Bismark CpG reports.
#'         \item \code{\link{getCpGtotals()}} to generate \code{CpGtotals}.
#'         \item \code{\link{filterCpGs()}} to filter the
#'                 \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#' }
#'
#' @examples \dontrun{
#'
#' # Read Bismark CpG Reports
#' colData <- read.xlsx("sample_info.xlsx", rowNames = TRUE)
#' bs <- getCpGs(colData, file = "Unfiltered_BSseq.rds")
#'
#' # Examine CpG Totals at Different Cutoffs
#' CpGtotals <- getCpGtotals(bs, file = "CpG_Totals.txt")
#' plotCpGtotals(CpGtotals, file = "CpG_Totals.pdf")
#'
#' # Filter BSseq Object
#' bs <- filterCpGs(bs, cov = 2, perSample = 0.75, file = "Filtered_BSseq.rds")
#' }
#'
#' @export
#'
#' @import ggplot2
#' @importFrom scales breaks_pretty

plotCpGtotals <- function(CpGtotals, nBreaks = 4,
                          legend.position = c(1.08,0.73), save = TRUE,
                          file = "CpG_Totals.pdf", width = 11, height = 4.25,
                          verbose = TRUE){
        if(verbose){
                message("[plotCpGtotals] Plotting CpG totals")
        }
        CpGtotals$perSample <- CpGtotals$perSample * 100
        gg <- ggplot(data = CpGtotals)
        gg <- gg +
                geom_line(aes(x = perSample, y = nCpGs_M, group = cov,
                              color = cov)) +
                geom_text(data = subset(CpGtotals, perSample == min(perSample)),
                          aes(x = perSample, y = nCpGs_M, group = cov,
                              color = cov, label = cov),
                          size = 4.5, check_overlap = TRUE, nudge_x = -0.5,
                          hjust = 1) +
                xlab("Samples (%) Cutoff") +
                ylab("Total CpGs (Millions)") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks + 1),
                                   expand = expansion(mult = c(0.05, 0.03))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_color_gradient("Coverage\nCutoff",
                                     breaks = breaks_pretty(n = nBreaks - 1)) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 14, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"),
                      axis.title = element_text(size = 18),
                      legend.background = element_blank(),
                      legend.position = legend.position,
                      legend.title = element_text(size = 18),
                      legend.text = element_text(size = 14),
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(),
                      plot.margin = unit(c(1,7,0.7,0.7), "lines"))
        if(save){
                if(verbose){
                        message("[plotCpGtotals] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width,
                       height = height, units = "in")
        }
        return(gg)
}

#' Filter BSseq Objects by Coverage
#'
#' \code{filterCpGs()} subsets a \code{\link[bsseq:`BSseq-class`]{BSseq}} object
#'         to include only those CpGs meeting \code{cov} and \code{perSample}
#'         cutoffs and then saves it as a .rds file.
#'
#' \code{filterCpGs()} is designed to be used after \code{cov} and
#'         \code{perSample} arguments have been optimized by
#'         \code{\link{getCpGtotals()}} and \code{\link{plotCpGtotals()}}.
#'
#' @param bs A \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#' @param cov A \code{numeric(1)} specifying the minimum number of reads
#'         overlapping a CpG for it to be included in the total.
#' @param perSample A \code{numeric(1)} specifying the minimum percent of
#'         samples with \code{cov} reads at a CpG for it to be included in the
#'         total.
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#' @param file A \code{character(1)} giving the file name (.rds) for the saved
#'         BSseq object.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{\link[bsseq:`BSseq-class`]{BSseq}} object.
#'
#' @seealso \itemize{
#'         \item \code{\link{getCpGs()}} to generate the
#'                 \code{\link[bsseq:`BSseq-class`]{BSseq}} object from
#'                 individual Bismark CpG reports.
#'         \item \code{\link{getCpGtotals()}} and \code{\link{plotCpGtotals()}}
#'                 for help with deciding coverage cutoffs.
#'         \item \code{\link{getRegions()}} to generate a set of regions based
#'                 on the CpGs.
#' }
#'
#' @examples \dontrun{
#'
#' # Read Bismark CpG Reports
#' colData <- read.xlsx("sample_info.xlsx", rowNames = TRUE)
#' bs <- getCpGs(colData, file = "Unfiltered_BSseq.rds")
#'
#' # Examine CpG Totals at Different Cutoffs
#' CpGtotals <- getCpGtotals(bs, file = "CpG_Totals.txt")
#' plotCpGtotals(CpGtotals, file = "CpG_Totals.pdf")
#'
#' # Filter BSseq Object
#' bs <- filterCpGs(bs, cov = 2, perSample = 0.75, file = "Filtered_BSseq.rds")
#' }
#'
#' @export
#'
#' @import bsseq
#' @importFrom magrittr %>%

filterCpGs <- function(bs, cov = 2, perSample = 0.75, save = TRUE,
                       file = "Filtered_BSseq.rds", verbose = TRUE){
        if(verbose){
                message("[getCpGs] Filtering CpG-level data for loci with at ",
                        "least ", cov, " reads in at least ", perSample * 100,
                        "% of samples")
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

.nCpGsBySample <- function(covSample, nSample){
        n <- length(covSample[covSample >= nSample])
        return(n)
}

.nCpGsByCov <- function(bsCov, cov, nSample){
        covSample <- (bsCov >= cov) %>% DelayedMatrixStats::rowSums2()
        nCpGs <- sapply(nSample, FUN = .nCpGsBySample, covSample = covSample)
        return(nCpGs)
}
