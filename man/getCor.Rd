% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualize_Comethylation_Network.R
\name{getCor}
\alias{getCor}
\title{Calculate Correlations}
\usage{
getCor(
  x,
  y = NULL,
  transpose = FALSE,
  corType = c("bicor", "pearson"),
  maxPOutliers = 0.1,
  robustY = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{x}{A \code{numeric vector} or \code{matrix}. \code{x} must be a
\code{matrix} if \code{y} is null.}

\item{y}{A \code{numeric vector} or \code{matrix}. If null, correlations will
be calculated for columns of \code{x}.}

\item{transpose}{A \code{logical(1)} specifying whether to transpose the
\code{matrix}.}

\item{corType}{A \code{character(1)} indicating which correlation statistic
to use in the calculation. Potential values include \code{pearson}
and \code{bicor}.}

\item{maxPOutliers}{A \code{numeric(1)} specifying the maximum percentile that
can be considered outliers on each side of the median for the
\code{bicor} statistic.}

\item{robustY}{A \code{logical(1)} indicating whether to use robust
calculation for \code{y} for the \code{bicor} statistic. \code{FALSE}
is recommended if \code{y} is a binary variable.}

\item{verbose}{A \code{logical(1)} indicating whether messages should be
printed.}
}
\value{
A \code{numeric matrix}.
}
\description{
\code{getCor()} calculates correlation coefficients using either
\code{pearson} or \code{bicor} methods. Calculations can be done between
columns of a single matrix or between two vectors or matrices.
}
\details{
The first input argument can be optionally transposed. The correlation
calculations are performed by \code{\link[WGCNA:cor]{WGCNA::cor()}} and \code{\link[WGCNA:bicor]{WGCNA::bicor()}}.
}
\examples{
\dontrun{

# Get Comethylation Modules
modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
                      corType = "pearson", file = "Modules.rds")

# Examine Correlations between Modules
MEs <- modules$MEs
moduleDendro <- getDendro(MEs, distance = "bicor")
plotDendro(moduleDendro, labelSize = 4, nBreaks = 5,
           file = "Module_ME_Dendrogram.pdf")
moduleCor <- getCor(MEs, corType = "bicor")
plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro,
            file = "Module_Correlation_Heatmap.pdf")
moduleCorStats <- getMEtraitCor(MEs, colData = MEs, corType = "bicor",
                                robustY = TRUE,
                                file = "Module_Correlation_Stats.txt")

# Examine Correlations between Samples
sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
plotDendro(sampleDendro, labelSize = 3, nBreaks = 5,
           file = "Sample_ME_Dendrogram.pdf")
sampleCor <- getCor(MEs, transpose = TRUE, corType = "bicor")
plotHeatmap(sampleCor, rowDendro = sampleDendro, colDendro = sampleDendro,
            file = "Sample_Correlation_Heatmap.pdf")

# Visualize Module Eigennode Values
plotHeatmap(MEs, rowDendro = sampleDendro, colDendro = moduleDendro,
            legend.title = "Module\nEigennode",
            legend.position = c(0.37,0.89),
            file = "Sample_ME_Heatmap.pdf")
}

}
\seealso{
\itemize{
\item \code{\link[=getModules]{getModules()}} to build a comethylation network and identify
modules of comethylated regions.
\item \code{\link[=getDendro]{getDendro()}} and \code{\link[=plotDendro]{plotDendro()}} to generate and visualize
dendrograms.
\item \code{\link[=plotHeatmap]{plotHeatmap()}} to visualize correlations between samples and
modules.
\item \code{\link[=getMEtraitCor]{getMEtraitCor()}} to calculate pairwise correlation
coefficients and p-values between module eigennode values.
}
}
