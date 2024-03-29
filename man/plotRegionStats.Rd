% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_and_Filter_Regions.R
\name{plotRegionStats}
\alias{plotRegionStats}
\title{Plot Histograms of Region Statistics}
\usage{
plotRegionStats(
  regions,
  maxQuantile = 1,
  bins = 30,
  histCol = "#132B43",
  lineCol = "red",
  nBreaks = 4,
  save = TRUE,
  file = "Region_Plots.pdf",
  width = 11,
  height = 8.5,
  verbose = TRUE
)
}
\arguments{
\item{regions}{A \code{data.frame} output from \code{\link[=getRegions]{getRegions()}} giving the set
of regions and statistics for each region.}

\item{maxQuantile}{A \code{numeric(1)} giving the maximum quantile of each
feature to plot.}

\item{bins}{A \code{numeric(1)} specifying the number of bins in each
histogram.}

\item{histCol}{A \code{character(1)} giving the color of the histogram.}

\item{lineCol}{A \code{character(1)} giving the color of the vertical median
line.}

\item{nBreaks}{A \code{numeric(1)} specifying the number of breaks for the
x-axis.}

\item{save}{A \code{logical(1)} indicating whether to save the plot.}

\item{file}{A \code{character(1)} giving the file name (.pdf) for the plot.}

\item{width}{A \code{numeric(1)} specifying the width in inches of the saved
plot.}

\item{height}{A \code{numeric(1)} specifying the height in inches of the saved
plot.}

\item{verbose}{A \code{logical(1)} indicating whether messages should be
printed.}
}
\value{
A \code{ggplot} object.
}
\description{
\code{plotRegionStats()} takes a set of regions from \code{\link[=getRegions]{getRegions()}},
generates histograms of region characteristics, and saves it as a pdf.
Region-level statistics include width, number of CpGs, minimum coverage, mean
coverage, mean methylation, and methylation standard deviation.
}
\details{
It's recommended to examine region characteristics before and after filtering.
The vertical line on each histogram indicates the median value for that
feature. A \code{ggplot} object is produced and can be edited outside of this
function if desired.
}
\examples{
\dontrun{

# Call Regions
regions <- getRegions(bs, file = "Unfiltered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99,
                file = "Unfiltered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99,
            file = "Unfiltered_SD_Plots.pdf")

# Examine Region Totals at Different Cutoffs
regionTotals <- getRegionTotals(regions, file = "Region_Totals.txt")
plotRegionTotals(regionTotals, file = "Region_Totals.pdf")

# Filter Regions
regions <- filterRegions(regions, covMin = 10, methSD = 0.05,
                         file = "Filtered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99,
                file = "Filtered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99,
            file = "Filtered_SD_Plots.pdf")
}

}
\seealso{
\itemize{
\item \code{\link[=getRegions]{getRegions()}} to generate the set of regions.
\item \code{\link[=plotSDstats]{plotSDstats()}}, \code{\link[=getRegionTotals]{getRegionTotals()}}, and \code{\link[=plotRegionTotals]{plotRegionTotals()}}
for more help visualizing region characteristics and setting
cutoffs for filtering.
\item \code{\link[=filterRegions]{filterRegions()}} for filtering regions by minimum coverage and
methylation standard deviation.
}
}
