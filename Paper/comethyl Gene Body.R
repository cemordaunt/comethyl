# comethyl Gene Body -------------------------------------------------------------------------------
# Charles Mordaunt
# 2/8/20

# Setup ####
setwd("~/Documents/Programming/comethyl/Testing/Gene Bodies")
.libPaths("/share/lasallelab/programs/comethyl/R_3.6")
AnnotationHub::setAnnotationHubOption("CACHE", value = "/share/lasallelab/programs/comethyl/R_3.6")
sapply(c("scales", "openxlsx", "rlist", "tidyverse", "ggdendro", "cowplot", "annotatr", "rtracklayer",
         "bsseq", "dmrseq", "WGCNA", "sva", "rGREAT", "R.devices", "biomaRt"), require, character.only = TRUE)
list.files("~/Documents/Programming/comethyl/R", full.names = TRUE) %>% lapply(source)

# Set Global Options ####
options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
enableWGCNAThreads(nThreads = 4)

# Load Previously Filtered BSseq Object ####
colData <- read.xlsx("sample_info.xlsx", rowNames = TRUE)
bs <- readRDS("comethyl_test/Filtered_BSseq.rds")

# Call Regions ####
regions <- getRegions(bs, annotation = "genes", genome = "hg38", file = "Unfiltered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99, file = "Unfiltered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99, file = "Unfiltered_SD_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.1, file = "Unfiltered_SD_Plots_Zoom.pdf")

# Examine Region Totals at Different Cutoffs ####
regionTotals <- getRegionTotals(regions, covMin = seq(0,100,10), methSD = seq(0,0.02,0.002),
                                file = "Region_Totals.txt") # Increase covMin range, decrease methSD range
plotRegionTotals(regionTotals, file = "Region_Totals.pdf")

# Filter Regions ####
# Try with same covMin filtering, but no methSD filter (only ~ 24k regions)
# covMin filter removes 670 regions, could increase this if desired
regions <- filterRegions(regions, covMin = 10, methSD = 0, file = "Filtered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99, file = "test.pdf") # "Filtered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99, file = "Filtered_SD_Plots.pdf")

# Adjust Methylation Data for PCs ####
meth <- getRegionMeth(regions, bs = bs, file = "Region_Methylation.rds")
mod <- model.matrix(~1, data = pData(bs))
methAdj <- adjustRegionMeth(meth, mod = mod, file = "Adjusted_Region_Methylation.rds") # Top 10 PCs
getDendro(methAdj, distance = "euclidean") %>% plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))
rm(bs)

# Select Soft Power Threshold ####
sft <- getSoftPower(methAdj, corType = "pearson", file = "Soft_Power.rds") # soft power = 12, took < 5 min
plotSoftPower(sft, file = "Soft_Power_Plots.pdf")

# Get Comethylation Modules ####
modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions, corType = "pearson",
                      file = "Modules.rds") # Took ~ 5 min
plotRegionDendro(modules, file = "Region_Dendrograms.pdf")
BED <- getModuleBED(modules$regions, file = "Modules.bed")

# Examine Correlations between Modules and Samples ####
MEs <- modules$MEs
moduleDendro <- getDendro(MEs, distance = "bicor")
plotDendro(moduleDendro, labelSize = 5, nBreaks = 5, file = "Module_ME_Dendrogram.pdf")
moduleCor <- getCor(MEs, corType = "bicor")
plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro,
            file = "Module_Correlation_Heatmap.pdf")

sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
plotDendro(sampleDendro, labelSize = 3, nBreaks = 5, file = "Sample_ME_Dendrogram.pdf")
sampleCor <- getCor(MEs, transpose = TRUE, corType = "bicor")
plotHeatmap(sampleCor, rowDendro = sampleDendro, colDendro = sampleDendro, file = "Sample_Correlation_Heatmap.pdf")

plotHeatmap(MEs, rowDendro = sampleDendro, colDendro = moduleDendro, legend.title = "Module\nEigennode",
            legend.position = c(0.37,0.89), file = "Sample_ME_Heatmap.pdf")

# Test Correlations between Module Eigennodes and Sample Traits ####
MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor", file = "ME_Trait_Correlation_Stats.txt")
traitDendro <- getCor(MEs, y = colData, corType = "bicor", robustY = FALSE) %>% getDendro(transpose = TRUE)
plotDendro(traitDendro, labelSize = 3.5, expandY = c(0.65,0.08), file = "Trait_Dendrogram.pdf")

plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order, traitOrder = traitDendro$order,
               file = "ME_Trait_Correlation_Heatmap.pdf")
plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order, traitOrder = traitDendro$order,
               sigOnly = TRUE, star.size = 11, star.nudge_y = -0.27, legend.position = c(1.14, 0.745),
               colColorMargins = c(-1,5.1,0.5,10.47), file = "Sig_ME_Trait_Correlation_Heatmap.pdf",
               width = 7, height = 3.5)

# Explore Significant ME-Trait Correlations ####
# Plot Module Eigennodes vs Traits
plotMEtraitDot(MEs$greenyellow, trait = colData$home_ownership, binwidth = 0.015, traitCode = c("No" = 0, "Yes" = 1),
               colors = c("No" = "#3366CC", "Yes" = "#FF3366"), ylim = c(-0.25,0.25), xlab = "Home Ownership",
               nBreaks = 5, ylab = "Green Yellow Module Eigennode", file = "greenyellow_ME_Home_Ownership_Dotplot.pdf")
plotMEtraitScatter(MEs$greenyellow, trait = colData$Gran, ylim = c(-0.25,0.25), xlab = "Granulocytes",
                   ylab = "Green Yellow Module Eigennode", file = "greenyellow_ME_Granulocytes_Scatterplot.pdf")
plotMEtraitScatter(MEs$greenyellow, trait = colData$CD4T, ylim = c(-0.25,0.25), xlab = "CD4+ T-Cells",
                   ylab = "Green Yellow Module Eigennode", file = "greenyellow_ME_CD4T_Scatterplot.pdf")
plotMEtraitScatter(MEs$greenyellow, trait = colData$CD8T, ylim = c(-0.25,0.25), xlab = "CD8+ T-Cells",
                   ylab = "Green Yellow Module Eigennode", file = "greenyellow_ME_CD8T_Scatterplot.pdf")

# Plot Region Methylation vs Traits
regions <- modules$regions
plotMethTrait("greenyellow", regions = regions, meth = meth, trait = colData$home_ownership,
              traitCode = c("No" = 0, "Yes" = 1), traitColors = c("No" = "#3366CC", "Yes" = "#FF3366"),
              trait.legend.title = "Home Ownership", trait.legend.position = c(1.05,4.39),
              file = "greenyellow_Module_Methylation_Home_Ownership_Heatmap.pdf")
plotMethTrait("greenyellow", regions = regions, meth = meth*100, trait = colData$Gran, expandY = 0.04,
              trait.legend.title = "Granulocytes", trait.legend.position = c(1.034,3.35),
              traitMargins = c(0,6,1,4.6), file = "greenyellow_Module_Methylation_Granulocytes_Heatmap.pdf")
plotMethTrait("greenyellow", regions = regions, meth = meth*100, trait = colData$CD4T, expandY = 0.04,
              trait.legend.title = "CD4+ T-Cells", trait.legend.position = c(1.038,3.35),
              traitMargins = c(0,6,1,4.6), file = "greenyellow_Module_Methylation_CD4T_Heatmap.pdf")
plotMethTrait("greenyellow", regions = regions, meth = meth*100, trait = colData$CD8T, expandY = 0.04,
              trait.legend.title = "CD8+ T-Cells", trait.legend.position = c(1.038,3.35),
              traitMargins = c(0,6,1,4.6), file = "greenyellow_Module_Methylation_CD8T_Heatmap.pdf")

# Annotate Modules ####
regionsAnno <- annotateModule(regions, module = "greenyellow", genome = "hg38",
                              file = "Annotated_greenyellow_Module_Regions.txt")
geneList <- getGeneList(regionsAnno, module = "greenyellow")
# Mapped with GREAT, includes overlapping and nearby, 14 regions, 26 genes
# [1] "AP3D1"    "ATOH1"    "CTIF"     "CUX1"     "DMRTB1"   "EGLN3"
# [7] "FAM184A"  "FOXN3"    "GLIS1"    "GRID2"    "GRM7"     "HSD17B12"
# [13] "IZUMO4"   "MCM9"     "NPAS3"    "NTN1"     "PMEPA1"   "RBMS3"
# [19] "SH2B2"    "STX8"     "TSHZ3"    "TTC17"    "TTC8"     "ZBP1"
# [25] "ZBTB7C"   "ZNF536"

# Analyze Functional Enrichment ####
enrichment <- enrichModule(regions, module = "greenyellow", genome = "hg38",
                           file = "greenyellow_Module_Enrichment.txt")
plotEnrichment(enrichment, file = "greenyellow_Module_Enrichment_Plot.pdf")

