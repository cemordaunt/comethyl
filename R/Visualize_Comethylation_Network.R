plotRegionDendro <- function(modules, save = TRUE, file = "Region_Dendrograms.pdf", width = 11, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotRegionDendro] Plotting region dendrograms and modules for each block")
        }
        blockColors <- lapply(modules$blockGenes, function(x) modules$colors[x])
        blockNames <- paste("Block ", 1:length(modules$dendrograms), " (", sapply(modules$blockGenes, length), " regions)", sep = "")
        if(save){
                if(verbose){
                        message("[plotRegionDendro] Saving plot as ", file)
                }
                pdf(file = file, width = width, height = height)
        }
        invisible(mapply(FUN = plotDendroAndColors, dendro = modules$dendrograms, colors = blockColors, main = blockNames,
                         MoreArgs = list(groupLabels = "Modules", dendroLabels = FALSE, marAll = c(1.5,5,3,1.5), saveMar = FALSE,
                                         cex.lab = 1.2, cex.colorLabels = 1.2, autoColorHeight = FALSE, lwd = 0.8,
                                         colorHeight = 0.15, cex.axis = 1, frame.plot = TRUE)))
        invisible(dev.off())
}

getModuleBED <- function(regions, grey = FALSE, save = TRUE, file = "Modules.bed", verbose = TRUE){
        if(!"module" %in% colnames(regions)){
                stop("[getModuleBED] Regions must have module annotation")
        }
        if(verbose){
                message("[getModuleBED] Creating bed file of regions annotated with identified modules")
        }
        if(!grey){
                if(verbose){
                        message("[getModuleBED] Excluding regions in grey (unassigned) module")
                }
                regions <- regions[!regions$module == "grey",]
        }
        regions$RegionID <- paste(regions$RegionID, regions$module, sep = "_")
        regions$rgb <- col2rgb(regions$module) %>% apply(2, paste, collapse = ",")
        BED <- cbind(regions[c("chr", "start", "end", "RegionID")], score = 0, strand = ".", thickStart = 0, thickEnd = 0,
                     rgb = regions$rgb)
        if(save){
                if(verbose){
                        message("[getModuleBED] Saving file as ", file)
                }
                name <- basename(file) %>% str_remove_all(pattern = ".bed")
                write(paste("track name='", name, "' description='", name, "' itemRgb='On'", sep = ""), file = file)
                write.table(BED, file = file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        }
        return(BED)
}

plotHeatmap <- function(x, rowDendro, colDendro, colors = blueWhiteRed(100, gamma = 0.3), limit = max(abs(x)),
                        axis.text.size = 8, legend.title = "Bicor", legend.title.size = 16, legend.text.size = 12,
                        legend.position = c(0.3,0.905), rowDendroMargins = c(-1.55,1,-0.1,-1.1),
                        colDendroMargins = c(1,-0.5,-1,0.8), rowColorMargins = c(-1.85,-1.5,0.55,1.7),
                        colColorMargins = c(-1.6,-0.85,1.8,0.55), save = TRUE, file = "Heatmap.pdf", width = 11,
                        height = 9.5, verbose = TRUE){
        if(verbose){
                message("[plotHeatmap] Plotting heatmap with dendrograms")
        }
        limits <- c(-limit, limit)
        x <- as.data.frame(x)
        rownames(x) <- str_remove_all(rownames(x), pattern = "ME")
        colnames(x) <- str_remove_all(colnames(x), pattern = "ME")
        rowModules <- sum(rownames(x) %in% colors()) == length(rownames(x))
        colModules <- sum(colnames(x) %in% colors()) == length(colnames(x))
        x$rowID <- factor(rownames(x), levels = rowDendro$labels[rev(rowDendro$order)])
        x <- reshape2::melt(x, id.vars = "rowID")
        x$variable <- factor(x$variable, levels = colDendro$labels[colDendro$order])
        hmMarginL <- ifelse(rowModules, yes = 2, no = -1)
        hmMarginB <- ifelse(colModules, yes = 2, no = -1)
        heatmap <- ggplot(data = x) +
                geom_tile(aes(x = variable, y = rowID, color = value, fill = value)) +
                scale_fill_gradientn(legend.title, colors = colors, limits = limits, aesthetics = c("color", "fill")) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_text(size = axis.text.size, color = "black", angle = 90, vjust = 0.5),
                      axis.text.y = element_text(size = axis.text.size, color = "black"),
                      axis.ticks = element_line(size = 0.8, color = "black"),
                      axis.title = element_blank(), legend.position = "none", panel.background = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25), panel.grid = element_blank(),
                      plot.background = element_blank(), plot.margin = unit(c(0,1,hmMarginB,hmMarginL), "lines"))
        legend <- get_legend(heatmap + theme(legend.position = legend.position, legend.background = element_blank(),
                                             legend.title = element_text(size = legend.title.size),
                                             legend.text = element_text(size = legend.text.size)))
        rowDendroPlot <- ggplot(data = dendro_data(rowDendro)$segments) +
                geom_segment(aes(x = -x, y = y, xend = -xend, yend = yend), lwd = 0.5, lineend = "square") +
                coord_flip() +
                theme_dendro() +
                theme(plot.margin = unit(rowDendroMargins, "lines"))
        colDendroPlot <- ggplot(data = dendro_data(colDendro)$segments) +
                geom_segment(aes(x = x, y = y, xend = xend, yend = yend), lwd = 0.5, lineend = "square") +
                theme_dendro() +
                theme(plot.margin = unit(colDendroMargins, "lines"))
        rowColors <- NULL
        colColors <- NULL
        if(rowModules){
                if(verbose){
                        message("[plotHeatmap] Using colors in row names for y-axis labels")
                }
                rowColors <- ggplot(data = data.frame(x = 0, y = 1:length(levels(x$rowID)), color = levels(x$rowID))) +
                        geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                        scale_fill_identity(aesthetics = c("color", "fill")) +
                        theme_void() +
                        theme(legend.position = "none", plot.margin = unit(rowColorMargins, "lines"))
                heatmap <- heatmap + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        }
        if(colModules){
                if(verbose){
                        message("[plotHeatmap] Using colors in column names for x-axis labels")
                }
                colColors <- ggplot(data = data.frame(x = 1:length(levels(x$variable)), y = 0, color = levels(x$variable))) +
                        geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                        scale_fill_identity(aesthetics = c("color", "fill")) +
                        theme_void() +
                        theme(legend.position = "none", plot.margin = unit(colColorMargins, "lines"))
                heatmap <- heatmap + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        }
        gg <- plot_grid(NULL, colDendroPlot, NULL, NULL, rowColors, heatmap, rowDendroPlot, legend, NULL, colColors,
                        NULL, NULL, nrow = 3, ncol = 4, rel_widths = c(0.045, 1, 0.15, 0.15),
                        rel_heights = c(0.15, 1, 0.045))
        if(save){
                if(verbose){
                        message("[plotHeatmap] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

getCor <- function(x, y = NULL, transpose = FALSE, corType = c("bicor", "pearson"), maxPOutliers = 0.1, robustY = TRUE,
                   verbose = TRUE){
        if(transpose){
                if(verbose){
                        message("[getCor] Transposing data")
                }
                x <- t(x)
        }
        corType <- match.arg(corType)
        if(verbose){
                message("[getCor] Calculating correlations using ", corType, " correlation")
        }
        if(corType == "bicor"){
                cor <- bicor(x, y = y, use = "pairwise.complete.obs", maxPOutliers = maxPOutliers,
                             robustY = robustY, pearsonFallback = "none")
        } else {
                if(corType == "pearson"){
                        cor <- WGCNA::cor(x, y = y, use = "pairwise.complete.obs")
                } else {
                        stop("[getCor] corType must be either bicor or pearson")
                }
        }
        return(cor)
}

