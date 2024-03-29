% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Analyze_Module_Preservation.R
\name{plotModulePreservation}
\alias{plotModulePreservation}
\title{Visualize Module Preservation}
\usage{
plotModulePreservation(
  preservation,
  line.size = 0.9,
  point.size = 1.2,
  nBreaks = 3,
  strip.text.size = 8.5,
  axis.text.size = 10,
  axis.title.size = 12,
  save = TRUE,
  file = "Module_Preservation_Plots.pdf",
  width = 9,
  height = 9,
  verbose = TRUE
)
}
\arguments{
\item{preservation}{A \code{data.frame} of module preservation statistics for
each module, generated by \code{getModulePreservation()}.}

\item{line.size}{A \code{numeric(1)} giving the size of the horizontal lines.}

\item{point.size}{A \code{numeric(1)} indicating the size of the points.}

\item{nBreaks}{A \code{numeric(1)} specifying the number of breaks used for
both axes.}

\item{strip.text.size}{A \code{numeric(1)} indicating the size of the title
text for each plot.}

\item{axis.text.size}{A \code{numeric(1)} specifying the size of the text for
both axes.}

\item{axis.title.size}{A \code{numeric(1)} indicating the size of the title
text for both axes.}

\item{save}{A \code{logical(1)} indicating whether to save the plot.}

\item{file}{A \code{character(1)} giving the file name (.pdf) for the saved
plot.}

\item{width}{A \code{numeric(1)} specifying the width in inches of the saved
plot.}

\item{height}{A \code{numeric(1)} specifying the height in inches of the
saved plot.}

\item{verbose}{A \code{logical(1)} indicating whether messages should be
printed.}
}
\value{
A \code{ggplot} object.
}
\description{
\code{plotModulePreservation()} plots Z-scores for various module
preservation statistics by module size as a scatterplot, and then saves it as
a .pdf.
}
\details{
\code{plotModulePreservation()} is designed to be used in combination with
\code{getModulePreservation()}. A blue line is plotted at Z = 2 to indicate
weak to moderate evidence for preservation, while a green line is plotted at
Z = 10 to indicate strong evidence. A \code{ggplot} object is produced and
can be edited outside of this function if desired. More information is given
in the documentation for \code{\link[WGCNA:modulePreservation]{WGCNA::modulePreservation()}}.
}
\examples{
\dontrun{
# Calculate Module Preservation
regions_disc <- modules_disc$regions
regions_rep <- modules_rep$regions
preservation <- getModulePreservation(methAdj_disc,
                                      regions_disc = regions_disc,
                                      meth_rep = methAdj_rep,
                                      regions_rep = regions_rep,
                                      corType = "pearson",
                                      file = "Module_Preservation_Stats.txt")

# Visualize Module Preservation
plotModulePreservation(preservation, file = "Module_Preservation_Plots.pdf")
}

}
\seealso{
\itemize{
\item \code{\link[=getModules]{getModules()}} to build a comethylation network and identify
modules of comethylated regions.
\item \code{\link[=getModulePreservation]{getModulePreservation()}} to calculate module preservation
statistics.
}
}
