# Comethyl CpG Clusters ----------------------------------------------------------------------------
# Charles Mordaunt
# 1/24/20

# Setup ####
setwd("~/Documents/Programming/comethyl/Testing/Default")
.libPaths("/share/lasallelab/programs/comethyl/R_3.6")
AnnotationHub::setAnnotationHubOption("CACHE", value = "/share/lasallelab/programs/comethyl/R_3.6")
sapply(c("scales", "openxlsx", "rlist", "tidyverse", "ggdendro", "cowplot", "annotatr", "rtracklayer",
         "bsseq", "dmrseq", "WGCNA", "sva", "rGREAT", "R.devices", "biomaRt"), require, character.only = TRUE)
list.files("~/Documents/Programming/comethyl/R", full.names = TRUE) %>% lapply(source)

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
regions <- filterRegions(regions, covMin = 10, methSD = 0.05, file = "Filtered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99, file = "Filtered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99, file = "Filtered_SD_Plots.pdf")

# Adjust Methylation Data for PCs ####
meth <- getRegionMeth(regions, bs = bs, file = "Region_Methylation.rds")
mod <- model.matrix(~1, data = pData(bs))
methAdj <- adjustRegionMeth(meth, mod = mod, file = "Adjusted_Region_Methylation.rds")
getDendro(methAdj, distance = "euclidean") %>% plotDendro(file = "Sample_Dendrogram.pdf", expandY = c(0.25,0.08))

# Select Soft Power Threshold ####
sft <- getSoftPower(methAdj, corType = "pearson", file = "Soft_Power.rds")
plotSoftPower(sft, file = "Soft_Power_Plots.pdf")

# Get Comethylation Modules ####
modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions, corType = "pearson", file = "Modules.rds")
plotRegionDendro(modules, file = "Region_Dendrograms.pdf")
BED <- getModuleBED(modules$regions, file = "Modules.bed")

# Examine Correlations between Modules and Samples ####
MEs <- modules$MEs
moduleDendro <- getDendro(MEs, distance = "bicor")
plotDendro(moduleDendro, labelSize = 4, nBreaks = 5, file = "Module_ME_Dendrogram.pdf")
moduleCor <- getCor(MEs, corType = "bicor")
plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro, file = "Module_Correlation_Heatmap.pdf")

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
plotMEtraitCor(MEtraitCor, moduleOrder = moduleDendro$order, traitOrder = traitDendro$order, sigOnly = TRUE, star.size = 11,
               star.nudge_y = -0.27, legend.position = c(1.14, 0.745), colColorMargins = c(-1,5.1,0.5,10.47),
               file = "Sig_ME_Trait_Correlation_Heatmap.pdf", width = 7, height = 3.5)

# Explore Significant ME-Trait Correlations ####
# Plot Module Eigennodes vs Traits
plotMEtraitDot(MEs$bisque4, trait = colData$Diagnosis_ASD, traitCode = c("TD" = 0, "ASD" = 1),
               colors = c("TD" = "#3366CC", "ASD" = "#FF3366"), ylim = c(-0.2,0.2), xlab = "Diagnosis",
               ylab = "Bisque 4 Module Eigennode", file = "bisque4_ME_Diagnosis_Dotplot.pdf")
plotMEtraitScatter(MEs$paleturquoise, trait = colData$Gran, ylim = c(-0.15,0.15), xlab = "Granulocytes",
                   ylab = "Pale Turquoise Module Eigennode", file = "paleturquoise_ME_Granulocytes_Scatterplot.pdf")
plotMEtraitScatter(MEs$paleturquoise, trait = colData$Bcell, ylim = c(-0.15,0.15), xlab = "B-cells",
                   ylab = "Pale Turquoise Module Eigennode", file = "paleturquoise_ME_Bcells_Scatterplot.pdf")

# Plot Region Methylation vs Traits
regions <- modules$regions
plotMethTrait("bisque4", regions = regions, meth = meth, trait = colData$Diagnosis_ASD, traitCode = c("TD" = 0, "ASD" = 1),
              traitColors = c("TD" = "#3366CC", "ASD" = "#FF3366"), trait.legend.title = "Diagnosis",
              file = "bisque4_Module_Methylation_Diagnosis_Heatmap.pdf")
plotMethTrait("paleturquoise", regions = regions, meth = meth, trait = colData$Gran, expandY = 0.04,
              trait.legend.title = "Granulocytes", trait.legend.position = c(1.034,3.35),
              file = "paleturquoise_Module_Methylation_Granulocytes_Heatmap.pdf")
plotMethTrait("paleturquoise", regions = regions, meth = meth, trait = colData$Bcell, expandY = 0.04,
              trait.legend.title = "B-cells", trait.legend.position = c(1.004,3.35),
              file = "paleturquoise_Module_Methylation_Bcells_Heatmap.pdf")

# Annotate Modules ####
regionsAnno <- annotateModule(regions, module = c("bisque4", "paleturquoise"), genome = "hg38",
                              file = "Annotated_bisque4_paleturquoise_Module_Regions.txt")
geneList_bisque4 <- getGeneList(regionsAnno, module = "bisque4")
geneList_paleturquoise <- getGeneList(regionsAnno, module = "paleturquoise")

# Analyze Functional Enrichment ####
ontologies <- listOntologies("hg38", version = "4.0.4")
enrich_bisque4 <- enrichModule(regions, module = "bisque4", genome = "hg38", file = "bisque4_Module_Enrichment.txt")
plotEnrichment(enrich_bisque4, file = "bisque4_Module_Enrichment_Plot.pdf")
enrich_paleturquoise <- enrichModule(regions, module = "paleturquoise", genome = "hg38",
                                     file = "paleturquoise_Module_Enrichment.txt")
plotEnrichment(enrich_paleturquoise, axis.text.y.size = 14, width = 10, file = "paleturquoise_Module_Enrichment_Plot.pdf")

# Circos Plot ####
library(circlize)
regions_col <- subset(regions, !module == "grey") %>% dplyr::select(chr, start, end, module)
pdf("Testing/Default/Module_Circos_Plot.pdf", width = 3.5, height = 3.5)
circos.par(gap.degree = 2, cell.padding = c(0.007,0,0.007,0), circle.margin = 0.00001)
circos.initializeWithIdeogram(species = "hg38", plotType = c("ideogram", "labels"))
circos.genomicTrack(regions_col, ylim = c(0, 1),
                    panel.fun = function(region, value, ...) {
                            circos.genomicRect(region, value, ytop = 1, ybottom = 0,
                                               border = unlist(value))
                    })
circos.clear()
dev.off()

# Module Region Counts ####
region_counts <- table(regions_col$module) %>% sort(decreasing = TRUE) %>% as.data.frame()
region_counts$Var1 <- as.character(region_counts$Var1) %>% factor(levels = rev(region_counts$Var1))
barplot <- ggplot(aes(x = Var1, y = Freq), data = region_counts) +
        geom_col(fill = "#132B43") +
        coord_flip() +
        scale_x_discrete(expand = expansion(c(0.02))) +
        scale_y_continuous(breaks = breaks_pretty(n = 4), expand = expansion(c(0.004, 0.03))) +
        theme_bw(base_size = 25) +
        theme(legend.position = "none", panel.grid.major = element_blank(),
              panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25),
              axis.ticks.y = element_blank(),
              panel.grid.minor = element_blank(), strip.background = element_blank(),
              axis.text.x = element_text(color = "black", size = 16),
              axis.text.y = element_blank(),
              axis.title.x = element_text(size = 20), axis.title.y = element_blank(),
              plot.margin = unit(c(1,1,0.5,0), "lines")) +
        ylab("Regions")
rowColors <- ggplot(data = data.frame(x = 0, y = 1:nrow(region_counts), color = rev(region_counts$Var1))) +
        geom_tile(aes(x = x, y = y, color = color, fill = color)) +
        scale_fill_identity(aesthetics = c("color", "fill")) +
        theme_void() +
        theme(legend.position = "none", plot.margin = unit(c(-0.15,-2,2.9,1), "lines"))
gg <- plot_grid(rowColors, barplot, ncol = 2, rel_widths = c(0.1, 1))
ggsave("Testing/Default/Module_Region_Counts.pdf", plot = gg, dpi = 600, width = 5, height = 7,
       units = "in")

