% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualize_Comethylation_Network.R
\name{getModuleBED}
\alias{getModuleBED}
\title{Get a Module BED file}
\usage{
getModuleBED(
  regions,
  grey = FALSE,
  save = TRUE,
  file = "Modules.bed",
  verbose = TRUE
)
}
\arguments{
\item{regions}{A \code{data.frame} of regions with module assignments,
typically obtained from \code{\link[=getModules]{getModules()}}.}

\item{grey}{A \code{logical(1)} specifying whether to include "grey"
(unassigned) regions in the BED file.}

\item{save}{A \code{logical(1)} indicating whether to save the BED file.}

\item{file}{A \code{character(1)} giving the file name (.BED).}

\item{verbose}{A \code{logical(1)} indicating whether messages should be
printed.}
}
\value{
A BED file.
}
\description{
\code{getModuleBED()} takes a \code{data.frame} of regions with module
annotations, converts it to the BED file format suitable for viewing it on
the UCSC Genome Browser, and then saves it.
}
\details{
\code{getModuleBED()} is designed to be used in combination with \code{\link[=getModules]{getModules()}}.
The BED file produced includes a header line to enable single-step viewing on
the UCSC Genome Browser. Each region is labeled by its \code{RegionID} and
assigned module, and is colored by the module color. "Grey" (unassigned)
regions are excluded by default, but can be optionally included.
}
\examples{
\dontrun{

# Get Comethylation Modules
modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
                      corType = "pearson", file = "Modules.rds")

# Visualize Comethylation Modules
plotRegionDendro(modules, file = "Region_Dendrograms.pdf")
BED <- getModuleBED(modules$regions, file = "Modules.bed")
}

}
\seealso{
\itemize{
\item \code{\link[=getModules]{getModules()}} to build a comethylation network and identify
modules of comethylated regions.
\item \code{\link[=plotRegionDendro]{plotRegionDendro()}} to visualize region similarity and module
assignments.
}
}
