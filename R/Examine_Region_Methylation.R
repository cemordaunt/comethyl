getRegionMeth <- function(regions, bs, type = "raw", save = TRUE, file = "Region_Methylation.rds", verbose = TRUE){
        if(verbose){
                message("[getRegionMeth] Calculating region methylation from BSseq object")
        }
        meth <- getMeth(bs, regions = regions[,c("chr", "start", "end")], type = type, what = "perRegion")
        rownames(meth) <- regions$RegionID
        if(save){
                if(verbose){
                        message("[getRegionMeth] Saving file as ", file)
                }
                saveRDS(meth, file = file)
        }
        return(meth)
}

adjustRegionMeth <- function(meth, mod = matrix(1, nrow = ncol(meth), ncol = 1), save = TRUE,
                             file = "Adjusted_Region_Methylation.rds", verbose = TRUE){
        if(verbose){
                message("[adjustRegionMeth] Determining number of principal components to adjust for")
        }
        n.pc <- num.sv(meth, mod = mod, seed = 5)
        if(verbose){
                message("[adjustRegionMeth] Adjusting region methylation for the top ", n.pc, " principal components")
        }
        methAdj <- sva_network(meth, n.pc = n.pc) %>% t()
        if(save){
                if(verbose){
                        message("[adjustRegionMeth] Saving file as ", file)
                }
                saveRDS(methAdj, file = file)
        }
        return(methAdj)
}

getDendro <- function(x, transpose = FALSE, distance = c("euclidean", "pearson", "bicor"), maxPOutliers = 0.1,
                      verbose = TRUE){
        if(transpose){
                if(verbose){
                        message("[getDendro] Transposing data")
                }
                x <- t(x)
        }
        distance <- match.arg(distance)
        if(distance == "euclidean"){
                if(verbose){
                        message("[getDendro] Clustering by euclidean distance")
                }
                dist <- dist(x)
        } else {
                if(distance == "pearson"){
                        if(verbose){
                                message("[getDendro] Clustering with pearson correlation as the distance")
                        }
                        dist <- (1 - WGCNA::cor(x, use = "pairwise.complete.obs")) %>% as.dist()
                } else {
                        if(distance == "bicor"){
                                if(verbose){
                                        message("[getDendro] Clustering with bicor correlation as the distance")
                                }
                                dist <- (1 - bicor(x, maxPOutliers = maxPOutliers, use = "pairwise.complete.obs")) %>%
                                        as.dist()
                        } else {
                                stop("[getDendro] Error: Distance must be either euclidean, pearson, or bicor")
                        }
                }
        }
        dendro <- hclust(dist, method = "average")
        dendro$labels <- str_remove_all(dendro$labels, pattern = "ME")
        return(dendro)
}

plotDendro <- function(dendro, label = TRUE, labelSize = 2.5, expandX = c(0.03,0.03), expandY = c(0.3,0.08),
                       nBreaks = 4, save = TRUE, file = "Dendrogram.pdf", width = 11, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotDendro] Plotting dendrogram")
        }
        dendroPlot <- dendro_data(dendro)
        fix <- dendroPlot$segments$yend == 0
        dendroPlot$segments$yend[fix] <- dendroPlot$segments$y[fix] - max(dendroPlot$segments$y) * 0.05
        dendroPlot$labels$y <- dendroPlot$segments$yend[fix] - max(dendroPlot$segments$y) * 0.01
        dendroPlot$labels$label <- str_remove_all(dendroPlot$labels$label, pattern = "ME")
        gg <- ggplot()
        gg <- gg +
                geom_segment(data = dendroPlot$segments, aes(x = x, y = y, xend = xend, yend = yend), lwd = 0.3, lineend = "square") +
                scale_x_continuous(expand = expansion(mult = expandX)) +
                scale_y_continuous(expand = expansion(mult = expandY), breaks = breaks_pretty(n = nBreaks)) +
                ylab("Height") +
                theme_dendro() +
                theme(plot.margin = unit(c(1,1,0,1), "lines"),
                      panel.background = element_rect(color = "black", fill = "white", size = 1.1),
                      axis.ticks.y = element_line(), axis.text.y = element_text(size = 12, color = "black"),
                      axis.title.y = element_text(size = 16, angle = 90, vjust = 2))
        if(label){
                gg <- gg +
                        geom_text(data = dendroPlot$labels, aes(x = x, y = y, label = label), angle = 90, hjust = 1,
                                  size = labelSize)
        }
        if(save){
                if(verbose){
                        message("[plotDendro] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

