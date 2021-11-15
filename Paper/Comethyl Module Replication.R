# Comethyl Module Replication --------------------------------------------------
# Charles Mordaunt
# 11/14/21

# Setup ------------------------------------------------------------------------
# export PATH=/share/lasallelab/programs/R_OpenBLAS/install/R/bin:${PATH} # R with OpenBLAS
# module load R # Standard
# R --vanilla

# Load Packages ####
.libPaths("/share/lasallelab/programs/comethyl/R_3.6") # OpenBLAS
#.libPaths("/share/lasallelab/Charles/R_4.1") # Standard
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication")
library(tidyverse)
library(WGCNA)
library(comethyl)

# Set Global Options ####
options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
enableWGCNAThreads(nThreads = 2)

# Create and Filter BSseq Object -----------------------------------------------
# Filter Sample Trait Table ####
colData <- read.csv("MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv") %>%
        filter(Study == "MARBLES" & Sex == "M" & Platform == "HiSeq4000") %>%
        select(Sequencing_ID, Diagnosis_ASD = Diagnosis_Alg, home_ownership,
               Bcell:Gran) %>%
        mutate(Diagnosis_ASD = ifelse(Diagnosis_ASD == "ASD", yes = 1, no = 0),
               home_ownership = ifelse(home_ownership == 99, yes = NA,
                                       no = home_ownership)) %>%
        arrange(Sequencing_ID)
rownames(colData) <- colData$Sequencing_ID
colData <- select(colData, !Sequencing_ID)
write.table(colData, file = "Replication_Sample_Info.txt", sep = "\t",
            quote = FALSE)

# Read Bismark CpG Reports ####
bs <- getCpGs(colData, file = "Replication_Unfiltered_BSseq.rds") # 29.4M CpGs
bs_disc <- readRDS("../Discovery/comethyl_test/Filtered_BSseq.rds") # 20.7M CpGs
bs <- IRanges::subsetByOverlaps(bs, ranges = bs_disc) # 20.7M CpGs (filtered for CpGs in discovery BSseq)
bs <- filterCpGs(bs, cov = 1, perSample = 0.01, # 20.6M CpGs (lost 102K CpGs with no reads, 0.5%)
                 file = "Replication_Filtered_BSseq.rds")
rm(bs_disc)

# Call Regions and Build Network Using CpG Clusters ----------------------------
# Call Regions ####
regions_disc <- readRDS("Module_Replication_CpG/Discovery_CpG_Cluster_Modules.rds") %>%
        .[["regions"]] %>% select(chr, start, end, RegionID) # Use discovery regions
regions_disc <- with(regions_disc,
                     GenomicRanges::GRanges(seqnames = chr,
                                            ranges = IRanges::IRanges(start = start, end = end),
                                            name = RegionID))
regions <- getRegions(bs, custom = regions_disc, n = 1, save = FALSE)
regions <- select(regions, RegionID = name, chr:methSD)
write_tsv(regions, file = "Module_Replication_CpG/Unfiltered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99,
                file = "Module_Replication_CpG/Unfiltered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99,
            file = "Module_Replication_CpG/Unfiltered_SD_Plots.pdf")

# Examine Region Totals at Different Cutoffs ####
regionTotals <- getRegionTotals(regions,
                                file = "Module_Replication_CpG/Region_Totals.txt")
plotRegionTotals(regionTotals, file = "Module_Replication_CpG/Region_Totals.pdf")

# Filter Regions ####
regions <- filterRegions(regions, covMin = 5, methSD = 0.05, # 188,090 / 251,717 regions in discovery
                         file = "Module_Replication_CpG/Filtered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99,
                file = "Module_Replication_CpG/Filtered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99,
            file = "Module_Replication_CpG/Filtered_SD_Plots.pdf")

# Adjust Methylation Data for PCs ####
meth <- getRegionMeth(regions, bs = bs,
                      file = "Module_Replication_CpG/Region_Methylation.rds")
mod <- model.matrix(~1, data = bsseq::pData(bs))
methAdj <- adjustRegionMeth(meth, mod = mod, # Top 16 PCs
                            file = "Module_Replication_CpG/Adjusted_Region_Methylation.rds")
getDendro(methAdj, distance = "euclidean") %>%
        plotDendro(file = "Module_Replication_CpG/Sample_Dendrogram.pdf",
                   expandY = c(0.25, 0.08))

# Select Soft Power Threshold ####
# Previously:
# At 20, fit = 0.76, connectivity = 33 (5000 modules)
# At 25, fit = 0.96, connectivity = 9 (2000 modules), for discovery used 18, fit = 0.86, connectivity = 8
# Now: At 19, fit = 0.822 and mean connectivity = 34
# Try 25, fit = 0.97, connectivity = 6 (1293 modules)
sft <- getSoftPower(methAdj, powerVector = 1:30, corType = "pearson",
                    file = "Module_Replication_CpG/Soft_Power.rds")
plotSoftPower(sft, file = "Module_Replication_CpG/Soft_Power_Plots.pdf")

# Get Comethylation Modules ####
modules <- getModules(methAdj, power = 25, regions = regions, corType = "pearson",
                      nThreads = 2, file = "Module_Replication_CpG/Modules.rds")

# Examine Correlations between Modules and Samples ####
MEs <- modules$MEs
moduleDendro <- getDendro(MEs, distance = "bicor")
plotDendro(moduleDendro, labelSize = 2.5, nBreaks = 5, width = 48, height = 8,
           expandX = c(0.015, 0.015), expandY = c(0.15, 0.08),
           file = "CpG Module Replication/Module_ME_Dendrogram.pdf")
moduleCor <- getCor(MEs, corType = "bicor")
plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro,
            axis.text.size = 4, rowDendroMargins = c(-1.55, 0.95, -0.1, -0.5),
            colDendroMargins = c(0.95, -0.5, -0.5, 0.8),
            file = "CpG Module Replication/Module_Correlation_Heatmap.pdf",
            width = 22, height = 19)

sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
plotDendro(sampleDendro, labelSize = 3, nBreaks = 5,
           file = "CpG Module Replication/Sample_ME_Dendrogram.pdf")
sampleCor <- getCor(MEs, transpose = TRUE, corType = "bicor")
plotHeatmap(sampleCor, rowDendro = sampleDendro, colDendro = sampleDendro,
            file = "CpG Module Replication/Sample_Correlation_Heatmap.pdf")

plotHeatmap(MEs, rowDendro = sampleDendro, colDendro = moduleDendro,
            legend.title = "Module\nEigennode", legend.position = c(0.37,0.89),
            file = "CpG Module Replication/Sample_ME_Heatmap.pdf",
            width = 33, axis.text.size = 6)

# Call Regions and Build Network Using Gene Bodies -----------------------------
# Setup ####
enableWGCNAThreads(nThreads = 4)
bs <- readRDS("Replication_Filtered_BSseq.rds")

# Call Regions ####
regions_disc <- readRDS("Module_Replication_Gene/Discovery_Gene_Body_Modules.rds") %>%
        .[["regions"]] %>% select(chr, start, end, strand, RegionID, gene_id) # Use discovery regions
regions_disc <- with(regions_disc,
                     GenomicRanges::GRanges(seqnames = chr,
                                            ranges = IRanges::IRanges(start = start, end = end),
                                            strand = strand, name = RegionID, gene_id = gene_id))
regions <- getRegions(bs, custom = regions_disc, n = 1, save = FALSE)
regions <- select(regions, RegionID = name, chr:methSD, strand, gene_id)
write_tsv(regions, file = "Module_Replication_Gene/Unfiltered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99,
                file = "Module_Replication_Gene/Unfiltered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99,
            file = "Module_Replication_Gene/Unfiltered_SD_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.1,
            file = "Module_Replication_Gene/Unfiltered_SD_Plots_Zoom.pdf")

# Examine Region Totals at Different Cutoffs ####
regionTotals <- getRegionTotals(regions, covMin = seq(0, 100, 10),
                                methSD = seq(0, 0.02, 0.002),
                                file = "Module_Replication_Gene/Region_Totals.txt")
plotRegionTotals(regionTotals, file = "Module_Replication_Gene/Region_Totals.pdf")

# Filter Regions ####
regions <- filterRegions(regions, covMin = 10, methSD = 0,
                         file = "Module_Replication_Gene/Filtered_Regions.txt")
plotRegionStats(regions, maxQuantile = 0.99,
                file = "Module_Replication_Gene/Filtered_Region_Plots.pdf")
plotSDstats(regions, maxQuantile = 0.99,
            file = "Module_Replication_Gene/Filtered_SD_Plots.pdf")

# Adjust Methylation Data for PCs ####
meth <- getRegionMeth(regions, bs = bs,
                      file = "Module_Replication_Gene/Region_Methylation.rds")
mod <- model.matrix(~1, data = bsseq::pData(bs))
methAdj <- adjustRegionMeth(meth, mod = mod,
                            file = "Module_Replication_Gene/Adjusted_Region_Methylation.rds") # Top 9 PCs
getDendro(methAdj, distance = "euclidean") %>%
        plotDendro(file = "Module_Replication_Gene/Sample_Dendrogram.pdf",
                   expandY = c(0.25, 0.08))

# Select Soft Power Threshold ####
sft <- getSoftPower(methAdj, corType = "pearson",
                    file = "Module_Replication_Gene/Soft_Power.rds") # At 17, fit = 0.85 and mean connectivity = 3.3
plotSoftPower(sft, file = "Module_Replication_Gene/Soft_Power_Plots.pdf")

# Get Comethylation Modules ####
modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
                      corType = "pearson", file = "Module_Replication_Gene/Modules.rds")
plotRegionDendro(modules, file = "Module_Replication_Gene/Region_Dendrograms.pdf")
BED <- getModuleBED(modules$regions, file = "Module_Replication_Gene/Modules.bed")

# Examine Correlations between Modules and Samples ####
MEs <- modules$MEs
moduleDendro <- getDendro(MEs, distance = "bicor")
plotDendro(moduleDendro, labelSize = 4, nBreaks = 5,
           file = "Module_Replication_Gene/Module_ME_Dendrogram.pdf")
moduleCor <- getCor(MEs, corType = "bicor")
plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro,
            file = "Module_Replication_Gene/Module_Correlation_Heatmap.pdf")

sampleDendro <- getDendro(MEs, transpose = TRUE, distance = "bicor")
plotDendro(sampleDendro, labelSize = 3, nBreaks = 5,
           file = "Module_Replication_Gene/Sample_ME_Dendrogram.pdf")
sampleCor <- getCor(MEs, transpose = TRUE, corType = "bicor")
plotHeatmap(sampleCor, rowDendro = sampleDendro, colDendro = sampleDendro,
            file = "Module_Replication_Gene/Sample_Correlation_Heatmap.pdf")

plotHeatmap(MEs, rowDendro = sampleDendro, colDendro = moduleDendro,
            legend.title = "Module\nEigennode", legend.position = c(0.37,0.89),
            file = "Module_Replication_Gene/Sample_ME_Heatmap.pdf")

# Compare CpG Cluster Modules Between Discovery and Replication ----------------
# Setup ####
setwd("~/Documents/Programming/comethyl/Testing")
modules_disc <- readRDS("Default/Modules.rds")
methAdj_disc <- readRDS("Default/Adjusted_Region_Methylation.rds")
modules_rep <- readRDS("CpG Module Replication/Modules.rds")
methAdj_rep <- readRDS("CpG Module Replication/Adjusted_Region_Methylation.rds")

# Calculate Module Preservation ####
# Change columns from region IDs to genomic coordinates (region IDs don't correspond)
table(colnames(methAdj_disc) == modules_disc$regions$RegionID) # All TRUE
colnames(methAdj_disc) <- with(modules_disc$regions, paste(chr, start, end, sep = "_"))
table(colnames(methAdj_rep) == modules_rep$regions$RegionID) # All TRUE
colnames(methAdj_rep) <- with(modules_rep$regions, paste(chr, start, end, sep = "_"))
table(colnames(methAdj_disc) %in% colnames(methAdj_rep))
# FALSE   TRUE
# 63627 188090
table(colnames(methAdj_rep) %in% colnames(methAdj_disc))
# TRUE
# 188090

# Prepare multi-set inputs
multiData <- list(Discovery = list(data = methAdj_disc),
                  Replication = list(data = methAdj_rep))
multiColor <- list(Discovery = modules_disc$regions$module,
                   Replication = modules_rep$regions$module)
rm(methAdj_disc, methAdj_rep, modules_disc, modules_rep)

# Run modulePreservation
preservation <- modulePreservation(multiData, multiColor = multiColor,
                                   networkType = "signed", randomSeed = 5,
                                   goldName = "random", verbose = 3)
saveRDS(preservation, file = "CpG Module Replication/Module_Preservation.rds")

# Collect Results ####
categories <- c("quality", "preservation", "accuracy", "referenceSeparability",
                "testSeparability")
tables <- c("observed", "Z", "log.p")
pattern <- c("Z." = "", "Z" = "", "log.p." = "", "log.p" = "")
preservationStats <- lapply(categories, function(x){
        category <- lapply(tables, function(y){
                as.data.frame(preservation[[x]][[y]]$ref.Discovery$inColumnsAlsoPresentIn.Replication) %>%
                        rownames_to_column(var = "module") %>%
                        pivot_longer(!module & !moduleSize, names_to = "statistic",
                                     values_to = y) %>%
                        mutate(statistic = str_replace_all(statistic, pattern = fixed(pattern)))
        })
        names(category) <- tables
        category <- full_join(category$observed, y = category$Z,
                              by = c("module", "moduleSize", "statistic")) %>%
                full_join(y = category$log.p,
                          by = c("module", "moduleSize", "statistic")) %>%
                mutate(category = x)
})
preservationStats <- bind_rows(preservationStats) %>%
        filter(!statistic %in% c("medianRank.qual", "meanClusterCoeff.qual",
                                 "medianRank.pres", "medianRankConnectivity.pres",
                                 "medianRankDensity.pres", "cor.clusterCoeff",
                                 "meanClusterCoeff.pres", "coClustering")) %>%
        mutate(log.p = as.numeric(log.p)) %>%
        select(module, moduleSize, category, statistic:log.p) %>%
        arrange(module)
write_tsv(preservationStats,
          file = "CpG Module Replication/Module_Preservation_Statistics.txt")

# Summarize Results ####
# Which modules are good quality?
filter(preservationStats, statistic == "summary.qual" & Z > 10) %>% pull(module) # 52/55
filter(preservationStats, statistic == "summary.qual" & Z > 2 & Z <= 10) %>% pull(module) # ivory
filter(preservationStats, statistic == "summary.qual" & Z <= 2) %>% pull(module) # grey, random

# Which modules are preserved?
filter(preservationStats, statistic == "summary.pres" & Z > 10) %>% pull(module) # lightcyan
filter(preservationStats, statistic == "summary.pres" & Z > 2 & Z <= 10) %>% pull(module) # blue, floralwhite
filter(preservationStats, statistic == "summary.pres" & Z <= 2) %>% pull(module) # 52/55

# Which modules are accurate?
filter(preservationStats, statistic == "accuracy" & Z > 10) %>% pull(module) # "darkgreen"   "lightcyan"   "saddlebrown" "thistle1"    "white"
filter(preservationStats, statistic == "accuracy" & Z > 2 & Z <= 10) %>% pull(module) # 40/55
filter(preservationStats, statistic == "accuracy" & Z <= 2) %>% pull(module)
# [1] "black"          "brown4"         "darkslateblue"  "grey"           "ivory"
# [6] "palevioletred3"

# Which modules separate in the reference set?
filter(preservationStats, statistic == "separability.qual" & Z > 10) %>% pull(module) # none
filter(preservationStats, statistic == "separability.qual" & Z > 2 & Z <= 10) %>% pull(module)
# [1] "brown"          "cyan"           "darkgrey"       "darkolivegreen" "darkorange2"
# [6] "darkred"        "green"          "greenyellow"    "grey"           "ivory"
# [11] "magenta"        "mediumpurple3"  "paleturquoise"  "purple"         "salmon"
# [16] "sienna3"        "skyblue"        "steelblue"      "thistle2"       "turquoise"
# [21] "yellowgreen"
filter(preservationStats, statistic == "separability.qual" & Z <= 2) %>% pull(module) # 34/55

# Which modules separate in the test set?
filter(preservationStats, statistic == "separability.pres" & Z > 10) %>% pull(module) # none
filter(preservationStats, statistic == "separability.pres" & Z > 2 & Z <= 10) %>% pull(module) # none
filter(preservationStats, statistic == "separability.pres" & Z <= 2) %>% pull(module) # all

# Visualize Results ####
plot_data <- select(preservationStats, module, moduleSize, statistic, Z) %>%
        mutate(statistic = factor(statistic, levels = unique(statistic))) %>%
        filter(!module %in% c("grey", "random") & is.finite(Z))
scatterplot <- ggplot(plot_data) +
        geom_hline(yintercept = 2, lty = 2, color = "blue", size = 0.9) +
        geom_hline(yintercept = 10, lty = 2, color = "darkgreen", size = 0.9) +
        geom_point(aes(x = moduleSize, y = Z, color = module), size = 1.2) +
        facet_wrap(vars(statistic), nrow = 6, scales = "free_y") +
        scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
        scale_y_continuous(breaks = scales::breaks_pretty(n = 3),
                           expand = expansion(0.15)) +
        scale_color_identity() +
        theme_bw(base_size = 25) +
        theme(legend.position = "none", panel.grid.major = element_blank(),
              panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 0.9),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = 8.5, margin = unit(c(1,0,0.2,0), "lines")),
              axis.text = element_text(color = "black", size = 10),
              axis.title = element_text(size = 12),
              panel.spacing.x = unit(0.2, "lines"),
              panel.spacing.y = unit(-0.5, "lines"),
              plot.margin = unit(c(0,1,1,1), "lines")) +
        xlab("Module Size") +
        ylab("Z-score")
ggsave("CpG Module Replication/Module_Preservation_Statistics_Plot.pdf",
       plot = scatterplot, dpi = 600, width = 9, height = 9, units = "in")

# Compare Gene Modules Between Discovery and Replication -----------------------
# Setup ####
modules_disc <- readRDS("Gene Bodies/Modules.rds")
methAdj_disc <- readRDS("Gene Bodies/Adjusted_Region_Methylation.rds")
modules_rep <- readRDS("Gene Module Replication/Modules.rds")
methAdj_rep <- readRDS("Gene Module Replication/Adjusted_Region_Methylation.rds")

# Calculate Module Preservation ####
# Change columns from region IDs to gene IDs
table(colnames(methAdj_disc) == modules_disc$regions$RegionID) # All TRUE
colnames(methAdj_disc) <- paste("gene", modules_disc$regions$gene_id, sep = "_")
table(colnames(methAdj_rep) == modules_rep$regions$RegionID) # All TRUE
colnames(methAdj_rep) <- paste("gene", modules_rep$regions$gene_id, sep = "_")
table(colnames(methAdj_disc) %in% colnames(methAdj_rep))
# FALSE  TRUE
#   332 22841
table(colnames(methAdj_rep) %in% colnames(methAdj_disc))
# TRUE
# 22841

# Prepare multi-set inputs
multiData <- list(Discovery = list(data = methAdj_disc),
                  Replication = list(data = methAdj_rep))
multiColor <- list(Discovery = modules_disc$regions$module,
                   Replication = modules_rep$regions$module)
rm(methAdj_disc, methAdj_rep, modules_disc, modules_rep)

# Run modulePreservation
preservation <- modulePreservation(multiData, multiColor = multiColor,
                                   networkType = "signed", randomSeed = 5,
                                   goldName = "random", verbose = 3)
saveRDS(preservation, file = "Gene Module Replication/Module_Preservation.rds")

# Collect Results ####
categories <- c("quality", "preservation", "accuracy", "referenceSeparability",
                "testSeparability")
tables <- c("observed", "Z", "log.p")
pattern <- c("Z." = "", "Z" = "", "log.p." = "", "log.p" = "")
preservationStats <- lapply(categories, function(x){
        category <- lapply(tables, function(y){
                as.data.frame(preservation[[x]][[y]]$ref.Discovery$inColumnsAlsoPresentIn.Replication) %>%
                        rownames_to_column(var = "module") %>%
                        pivot_longer(!module & !moduleSize, names_to = "statistic",
                                     values_to = y) %>%
                        mutate(statistic = str_replace_all(statistic, pattern = fixed(pattern)))
        })
        names(category) <- tables
        category <- full_join(category$observed, y = category$Z,
                              by = c("module", "moduleSize", "statistic")) %>%
                full_join(y = category$log.p,
                          by = c("module", "moduleSize", "statistic")) %>%
                mutate(category = x)
})
preservationStats <- bind_rows(preservationStats) %>%
        filter(!statistic %in% c("medianRank.qual", "meanClusterCoeff.qual",
                                 "medianRank.pres", "medianRankConnectivity.pres",
                                 "medianRankDensity.pres", "cor.clusterCoeff",
                                 "meanClusterCoeff.pres", "coClustering")) %>%
        mutate(log.p = as.numeric(log.p)) %>%
        select(module, moduleSize, category, statistic:log.p) %>%
        arrange(module)
write_tsv(preservationStats,
          file = "Gene Module Replication/Module_Preservation_Statistics.txt")

# Summarize Results ####
# Which modules are good quality?
filter(preservationStats, statistic == "summary.qual" & Z > 10) %>% pull(module) # 15/17
filter(preservationStats, statistic == "summary.qual" & Z > 2 & Z <= 10) %>% pull(module) # None
filter(preservationStats, statistic == "summary.qual" & Z <= 2) %>% pull(module) # grey, random

# Which modules are preserved?
filter(preservationStats, statistic == "summary.pres" & Z > 10) %>% pull(module) # brown, pink
filter(preservationStats, statistic == "summary.pres" & Z > 2 & Z <= 10) %>% pull(module) # red
filter(preservationStats, statistic == "summary.pres" & Z <= 2) %>% pull(module) # 14/17

# Which modules are accurate?
filter(preservationStats, statistic == "accuracy" & Z > 10) %>% pull(module) # "brown" "grey"  "pink"
filter(preservationStats, statistic == "accuracy" & Z > 2 & Z <= 10) %>% pull(module)
# "black"       "cyan"        "greenyellow" "red"         "tan"
filter(preservationStats, statistic == "accuracy" & Z <= 2) %>% pull(module) # 8/16

# Which modules separate in the reference set?
filter(preservationStats, statistic == "separability.qual" & Z > 10) %>% pull(module) # None
filter(preservationStats, statistic == "separability.qual" & Z > 2 & Z <= 10) %>% pull(module) # brown, grey
filter(preservationStats, statistic == "separability.qual" & Z <= 2) %>% pull(module) # 15/17

# Which modules separate in the test set?
filter(preservationStats, statistic == "separability.pres" & Z > 10) %>% pull(module) # None
filter(preservationStats, statistic == "separability.pres" & Z > 2 & Z <= 10) %>% pull(module) # None
filter(preservationStats, statistic == "separability.pres" & Z <= 2) %>% pull(module) # All

# Visualize Results ####
plot_data <- select(preservationStats, module, moduleSize, statistic, Z) %>%
        mutate(statistic = factor(statistic, levels = unique(statistic))) %>%
        filter(!module %in% c("grey", "random") & is.finite(Z))
scatterplot <- ggplot(plot_data) +
        geom_hline(yintercept = 2, lty = 2, color = "blue", size = 0.9) +
        geom_hline(yintercept = 10, lty = 2, color = "darkgreen", size = 0.9) +
        geom_point(aes(x = moduleSize, y = Z, color = module), size = 1.5) +
        facet_wrap(vars(statistic), nrow = 6, scales = "free_y") +
        scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
        scale_y_continuous(breaks = scales::breaks_pretty(n = 3),
                           expand = expansion(0.15)) +
        scale_color_identity() +
        theme_bw(base_size = 25) +
        theme(legend.position = "none", panel.grid.major = element_blank(),
              panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 0.9),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = 8.5, margin = unit(c(1,0,0.2,0), "lines")),
              axis.text = element_text(color = "black", size = 10),
              axis.title = element_text(size = 12),
              panel.spacing.x = unit(0.2, "lines"),
              panel.spacing.y = unit(-0.5, "lines"),
              plot.margin = unit(c(0,1,1,1), "lines")) +
        xlab("Module Size") +
        ylab("Z-score")
ggsave("Gene Module Replication/Module_Preservation_Statistics_Plot.pdf",
       plot = scatterplot, dpi = 600, width = 9, height = 9, units = "in")

