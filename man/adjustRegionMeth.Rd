% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Examine_Region_Methylation.R
\name{adjustRegionMeth}
\alias{adjustRegionMeth}
\title{Adjust Methylation Data for Principal Components}
\usage{
adjustRegionMeth(
  meth,
  PCs,
  save = TRUE,
  file = "Adjusted_Region_Methylation.rds",
  verbose = TRUE
)
}
\arguments{
\item{meth}{A \code{numeric matrix}, where each row is a region and each
column is a sample. This is typically obtained from \code{\link[=getRegionMeth]{getRegionMeth()}}.}

\item{PCs}{A \code{numeric matrix}, where each row is a sample and each
column is a principal component. This is typically obtained from
\code{\link[=getPCs]{getPCs()}}.}

\item{save}{A \code{logical(1)} indicating whether to save the \code{matrix}.}

\item{file}{A \code{character(1)} giving the file name (.rds) for the saved
\code{matrix}.}

\item{verbose}{A \code{logical(1)} indicating whether messages should be
printed.}
}
\value{
A \code{numeric matrix}, where each row is a sample and each column
is a region.
}
\description{
\code{adjustRegionMeth()} adjusts region methylation data for the top
principal components, transposes it, and then saves it as a .rds file.
}
\details{
\code{adjustRegionMeth()} regresses out the top principal components
generated by \code{\link[=getPCs]{getPCs()}}. This is the same approach as taken by
\code{\link[sva:sva_network]{sva::sva_network()}}. More information on the function and approach is given
in the documentation and publications related to the \pkg{sva} package.
}
\examples{
\dontrun{

# Get Methylation Data
meth <- getRegionMeth(regions, bs = bs, file = "Region_Methylation.rds")

# Adjust Methylation Data for PCs
mod <- model.matrix(~1, data = pData(bs))
PCs <- getPCs(meth, mod = mod, file = "Top_Principal_Components.rds")
methAdj <- adjustRegionMeth(meth, PCs = PCs,
                            file = "Adjusted_Region_Methylation.rds")

# Compare Top PCs to Sample Traits
MEtraitCor <- getMEtraitCor(PCs, colData = colData, corType = "bicor",
                            file = "PC_Trait_Correlation_Stats.txt")
PCdendro <- getDendro(PCs, distance = "bicor")
PCtraitDendro <- getCor(PCs, y = colData, corType = "bicor", robustY = FALSE) \%>\%
        getDendro(transpose = TRUE)
plotMEtraitCor(PCtraitCor, moduleOrder = PCdendro$order,
               traitOrder = PCtraitDendro$order,
               file = "PC_Trait_Correlation_Heatmap.pdf")

# Assess Sample Similarity
getDendro(methAdj, distance = "euclidean") \%>\%
        plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))

# Select Soft Power Threshold
sft <- getSoftPower(methAdj, corType = "pearson", file = "Soft_Power.rds")
plotSoftPower(sft, file = "Soft_Power_Plots.pdf")

# Get Comethylation Modules
modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
                      corType = "pearson", file = "Modules.rds")
}

}
\seealso{
\itemize{
\item \code{\link[=getRegionMeth]{getRegionMeth()}} to extract region methylation values.
\item \code{\link[=getPCs]{getPCs()}} to calculate top principal components for region
methylation values.
\item \code{\link[=getMEtraitCor]{getMEtraitCor()}} to compare these top PCs to sample traits.
\item \code{\link[=getDendro]{getDendro()}} and \code{\link[=plotDendro]{plotDendro()}} to generate and visualize
dendrograms.
\item \code{\link[=getSoftPower]{getSoftPower()}} and \code{\link[=plotSoftPower]{plotSoftPower()}} to estimate the best
soft-thresholding power and visualize scale-free topology fit
and connectivity.
\item \code{\link[=getModules]{getModules()}} to build a comethylation network and identify
modules of comethylated regions.
}
}
