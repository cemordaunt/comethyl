#!/usr/bin/env Rscript

# Comethyl Analysis Pipeline ---------------------------------------------------
# Author: Charles Mordaunt
# Under Development

# Setup ####
.libPaths("/share/lasallelab/programs/comethyl/R_3.6")
AnnotationHub::setAnnotationHubOption("CACHE", value = "/share/lasallelab/programs/comethyl/R_3.6")
library(comethyl)

# Set Global Options ####
options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
enableWGCNAThreads(nThreads = 4)

# Read Bismark CpG Reports ####
colData <- read.xlsx("sample_info.xlsx", rowNames = TRUE)
bs <- getCpGs(colData, file = "Unfiltered_BSseq.rds")

# Examine CpG Totals at Different Cutoffs ####
CpGtotals <- getCpGtotals(bs, file = "CpG_Totals.txt")
plotCpGtotals(CpGtotals, file = "CpG_Totals.pdf")

# Filter BSobject ####
bs <- filterCpGs(bs, cov = 2, perSample = 0.75, file = "Filtered_BSseq.rds")

# Call Regions ####
regions <- getRegions(bs, file = "Unfiltered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99, file = "Unfiltered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99, file = "Unfiltered_SD_Plots.pdf")

# Examine Region Totals at Different Cutoffs ####
regionTotals <- getRegionTotals(regions, file = "Region_Totals.txt")
plotRegionTotals(regionTotals, file = "Region_Totals.pdf")

# Filter Regions ####
regions <- filterRegions(regions, covMin = 10, methSD = 0.05,
                         file = "Filtered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99, file = "Filtered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99, file = "Filtered_SD_Plots.pdf")

# Adjust Methylation Data for PCs ####
meth <- getRegionMeth(regions, bs = bs, file = "Region_Methylation.rds")
mod <- model.matrix(~1, data = pData(bs))
methAdj <- adjustRegionMeth(meth, mod = mod,
                            file = "Adjusted_Region_Methylation.rds")
getDendro(methAdj, distance = "euclidean") %>%
        plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))

# Select Soft Power Threshold ####
sft <- getSoftPower(methAdj, corType = "pearson", file = "Soft_Power.rds")
plotSoftPower(sft, file = "Soft_Power_Plots.pdf")

# Get Comethylation Modules ####
modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
                      corType = "pearson", file = "Modules.rds")
plotRegionDendro(modules, file = "Region_Dendrograms.pdf")
BED <- getModuleBED(modules$regions, file = "Modules.bed")

# Examine Correlations between Modules and Samples ####
MEs <- modules$MEs
moduleDendro <- getDendro(MEs, distance = "bicor")
plotDendro(moduleDendro, labelSize = 4, nBreaks = 5,
           file = "Module_ME_Dendrogram.pdf")
moduleCor <- getCor(MEs, corType = "bicor")
plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro,
            file = "Module_Correlation_Heatmap.pdf")

sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
plotDendro(sampleDendro, labelSize = 3, nBreaks = 5,
           file = "Sample_ME_Dendrogram.pdf")
sampleCor <- getCor(MEs, transpose = TRUE, corType = "bicor")
plotHeatmap(sampleCor, rowDendro = sampleDendro, colDendro = sampleDendro,
            file = "Sample_Correlation_Heatmap.pdf")

plotHeatmap(MEs, rowDendro = sampleDendro, colDendro = moduleDendro,
            legend.title = "Module\nEigennode",
            legend.position = c(0.37,0.89), file = "Sample_ME_Heatmap.pdf")

# Test Correlations between Module Eigennodes and Sample Traits ####
MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor",
                            file = "ME_Trait_Correlation_Stats.txt")
traitDendro <- getCor(MEs, y = colData, corType = "bicor", robustY = FALSE) %>%
        getDendro(transpose = TRUE)
plotDendro(traitDendro, labelSize = 3.5, expandY = c(0.65,0.08),
           file = "Trait_Dendrogram.pdf")
plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order,
               traitOrder = traitDendro$order,
               file = "ME_Trait_Correlation_Heatmap.pdf")
plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order,
               traitOrder = traitDendro$order, sigOnly = TRUE, star.size = 11,
               star.nudge_y = -0.27, legend.position = c(1.14, 0.745),
               colColorMargins = c(-1,5.1,0.5,10.47),
               file = "Sig_ME_Trait_Correlation_Heatmap.pdf", width = 7,
               height = 3.5)

# Explore Significant ME-Trait Correlations ####
# Plot Module Eigennodes vs Traits
plotMEtraitDot(MEs$bisque4, trait = colData$Diagnosis_ASD,
               traitCode = c("TD" = 0, "ASD" = 1),
               colors = c("TD" = "#3366CC", "ASD" = "#FF3366"), ylim = c(-0.2,0.2),
               xlab = "Diagnosis", ylab = "Bisque 4 Module Eigennode",
               file = "bisque4_ME_Diagnosis_Dotplot.pdf")
plotMEtraitScatter(MEs$paleturquoise, trait = colData$Gran, ylim = c(-0.15,0.15),
                   xlab = "Granulocytes", ylab = "Pale Turquoise Module Eigennode",
                   file = "paleturquoise_ME_Granulocytes_Scatterplot.pdf")
plotMEtraitScatter(MEs$paleturquoise, trait = colData$Bcell, ylim = c(-0.15,0.15),
                   xlab = "B-cells", ylab = "Pale Turquoise Module Eigennode",
                   file = "paleturquoise_ME_Bcells_Scatterplot.pdf")

# Plot Region Methylation vs Traits
regions <- modules$regions
plotMethTrait("bisque4", regions = regions, meth = meth,
              trait = colData$Diagnosis_ASD, traitCode = c("TD" = 0, "ASD" = 1),
              traitColors = c("TD" = "#3366CC", "ASD" = "#FF3366"),
              trait.legend.title = "Diagnosis",
              file = "bisque4_Module_Methylation_Diagnosis_Heatmap.pdf")
plotMethTrait("paleturquoise", regions = regions, meth = meth,
              trait = colData$Gran, expandY = 0.04, trait.legend.title = "Granulocytes",
              trait.legend.position = c(1.034,3.35),
              file = "paleturquoise_Module_Methylation_Granulocytes_Heatmap.pdf")
plotMethTrait("paleturquoise", regions = regions, meth = meth,
              trait = colData$Bcell, expandY = 0.04,
              trait.legend.title = "B-cells", trait.legend.position = c(1.004,3.35),
              file = "paleturquoise_Module_Methylation_Bcells_Heatmap.pdf")

# Annotate Modules ####
regionsAnno <- annotateModule(regions, module = c("bisque4", "paleturquoise"),
                              genome = "hg38",
                              file = "Annotated_bisque4_paleturquoise_Module_Regions.txt")
geneList_bisque4 <- getGeneList(regionsAnno, module = "bisque4")
geneList_paleturquoise <- getGeneList(regionsAnno, module = "paleturquoise")

# Analyze Functional Enrichment ####
ontologies <- listOntologies("hg38", version = "4.0.4")
enrich_bisque4 <- enrichModule(regions, module = "bisque4", genome = "hg38",
                               file = "bisque4_Module_Enrichment.txt")
plotEnrichment(enrich_bisque4, file = "bisque4_Module_Enrichment_Plot.pdf")
enrich_paleturquoise <- enrichModule(regions, module = "paleturquoise",
                                     genome = "hg38",
                                     file = "paleturquoise_Module_Enrichment.txt")
plotEnrichment(enrich_paleturquoise, axis.text.y.size = 14, width = 10,
               file = "paleturquoise_Module_Enrichment_Plot.pdf")

