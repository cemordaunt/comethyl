% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Investigate_Functional_Enrichments.R
\name{listOntologies}
\alias{listOntologies}
\title{Get Ontologies Available in GREAT}
\usage{
listOntologies(
  genome = c("hg38", "hg19", "hg18", "mm10", "mm9", "danRer7"),
  version = c("4.0.4", "3.0.0", "2.0.2"),
  verbose = TRUE
)
}
\arguments{
\item{genome}{A \code{character(1)} giving the genome build. Supported genomes
include \code{hg38}, \code{hg19}, \code{hg18}, \code{mm10}, \code{mm9},
and \code{danRer7}. See \code{Details}.}

\item{version}{A \code{character(1)} specifying the version of GREAT. Possible
values include \code{4.0.4}, \code{3.0.0}, and \code{2.0.2}. Different
versions of GREAT support different genomes. See \code{Details}.}

\item{verbose}{A \code{logical(1)} indicating whether messages should be
printed.}
}
\value{
A \code{character vector}.
}
\description{
\code{listOntologies()} obtains a \code{character vector} of the available
ontologies for functional enrichment analysis by GREAT for a specified version
and genome build.
}
\details{
\code{listOntologies()} generates the possible ontologies to use for
functional enrichment analysis in \code{\link[=enrichModule]{enrichModule()}}. Supported ontologies may
change over time, so this function queries GREAT using the \pkg{rGREAT}
package to get the ones currently available.

GREAT supports different genomes depending on the version:
\describe{
\item{GREAT v4.0.4}{\code{hg38}, \code{hg19}, \code{mm10}, \code{mm9}}
\item{GREAT v3.0.0}{\code{hg19}, \code{mm10}, \code{mm9}, \code{danRer7}}
\item{GREAT v2.0.2}{\code{hg19}, \code{hg18}, \code{mm9}, \code{danRer7}}
}
}
\examples{
\dontrun{

# Get Comethylation Modules
modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
                      corType = "pearson", file = "Modules.rds")

# Annotate Modules
regionsAnno <- annotateModule(regions, module = "bisque4",
                              genome = "hg38",
                              file = "Annotated_bisque4_Module_Regions.txt")
geneList_bisque4 <- getGeneList(regionsAnno, module = "bisque4")

# Analyze Functional Enrichment
ontologies <- listOntologies("hg38", version = "4.0.4")
enrich_bisque4 <- enrichModule(regions, module = "bisque4",
                               genome = "hg38",
                               file = "bisque4_Module_Enrichment.txt")
plotEnrichment(enrich_bisque4,
               file = "bisque4_Module_Enrichment_Plot.pdf")
}

}
\seealso{
\itemize{
\item \code{\link[=getModules]{getModules()}} to build a comethylation network and identify
modules of comethylated regions.
\item \code{\link[=annotateModule]{annotateModule()}} and \code{\link[=getGeneList]{getGeneList()}} to annotate a set of
regions with genes and regulatory context and then extract
the gene symbols or IDs.
\item \code{\link[=enrichModule]{enrichModule()}} and \code{\link[=plotEnrichment]{plotEnrichment()}} to investigate
functional enrichment of module regions with GREAT.
}
}
