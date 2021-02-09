getMEtraitCor <- function(MEs, colData, corType = c("bicor", "pearson"), maxPOutliers = 0.1, robustY = FALSE,
                          adjMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
                          save = TRUE, file = "ME_Trait_Correlation_Stats.txt", verbose = TRUE){
        corType <- match.arg(corType)
        adjMethod <- match.arg(adjMethod)
        if(verbose){
                message("[getMEtraitCor] Testing associations between module eigennodes and sample traits using ",
                        corType, " correlation")
        }
        colData <- colData[rownames(MEs),]
        if(corType == "bicor"){
                cor <- bicorAndPvalue(x = MEs, y = colData, maxPOutliers = maxPOutliers, robustY = robustY,
                                      pearsonFallback = "none")
        } else {
                if(corType == "pearson"){
                        cor <- corAndPvalue(x = MEs, y = colData)
                } else {
                        stop("[getMEtraitCor] corType must be either bicor or pearson")
                }
        }
        stats <- list.rbind(cor) %>% as.data.frame()
        stats$module <- rownames(cor$p) %>% str_remove_all(pattern = "ME") %>% rep(length(cor)) %>%
                factor(levels = unique(.))
        stats$stat <- names(cor) %>% rep(each = nrow(cor$p)) %>% factor(levels = unique(.))
        stats <- reshape2::melt(stats, id.vars = c("module", "stat"), variable.name = "trait") %>%
                reshape2::dcast(formula = module + trait ~ stat, value.var = "value")
        if(verbose){
                message("[getMEtraitCor] Adjusting p-values using the ", adjMethod, " method")
        }
        stats$adj_p <- p.adjust(stats$p, method = adjMethod)
        if(corType == "bicor"){
                stats <- stats[,c("module", "trait", "nObs", "bicor", "Z", "t", "p", "adj_p")]
        } else {
                stats <- stats[,c("module", "trait", "nObs", "cor", "Z", "t", "p", "adj_p")]
        }
        if(save){
                if(verbose){
                        message("[getMEtraitCor] Saving file as ", file)
                }
                write.table(stats, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
        }
        return(stats)
}

plotMEtraitCor <- function(MEtraitCor, moduleOrder = 1:length(unique(MEtraitCor$module)),
                           traitOrder = 1:length(unique(MEtraitCor$trait)), sigOnly = FALSE, adj_p = 0.05, star.size = 8,
                           star.nudge_y = -0.38, colors = blueWhiteRed(100, gamma = 0.9), limit = NULL,
                           axis.text.size = 12, legend.position = c(1.08, 0.915), legend.text.size = 12,
                           legend.title.size = 16, colColorMargins = c(-0.7,4.21,1.2,11.07), save = TRUE,
                           file = "ME_Trait_Correlation_Heatmap.pdf", width = 11, height = 9.5, verbose = TRUE){
        if(verbose){
                message("[plotMEtraitCor] Plotting ME trait correlation heatmap")
        }
        MEtraitCor$module <- factor(MEtraitCor$module, levels = levels(MEtraitCor$module)[moduleOrder])
        MEtraitCor$trait <- factor(MEtraitCor$trait, levels = levels(MEtraitCor$trait)[rev(traitOrder)])
        MEtraitCor$Significant <- (MEtraitCor$adj_p < adj_p & !is.na(MEtraitCor$adj_p)) %>% factor(levels = c("TRUE", "FALSE"))
        if(sigOnly){
                sigModules <- MEtraitCor$module[MEtraitCor$Significant == "TRUE"] %>% unique() %>% as.character()
                sigTraits <- MEtraitCor$trait[MEtraitCor$Significant == "TRUE"] %>% unique() %>% as.character()
                MEtraitCor <- subset(MEtraitCor, module %in% sigModules & trait %in% sigTraits)
                MEtraitCor$module <- factor(MEtraitCor$module,
                                            levels = levels(MEtraitCor$module)[levels(MEtraitCor$module) %in% sigModules])
        }
        if("bicor" %in% colnames(MEtraitCor)){
                corType <- "bicor"
        } else {
                if("cor" %in% colnames(MEtraitCor)){
                        corType <- "cor"
                } else {
                        stop("[plotMEtraitCor] corType unknown, must be either bicor or cor")
                }
        }
        corData <- MEtraitCor[[corType]]
        if(is.null(limit)){
                limit <- max(abs(corData))
        }
        heatmap <- ggplot(data = MEtraitCor) +
                geom_tile(aes(x = module, y = trait, color = corData, fill = corData)) +
                geom_text(aes(x = module, y = trait, alpha = Significant), label = "*", color = "black",
                          size = star.size, nudge_y = star.nudge_y) +
                scale_fill_gradientn(str_to_title(corType), colors = colors, limits = c(-limit, limit),
                                     aesthetics = c("color", "fill")) +
                scale_x_discrete(expand = expansion(mult = 0.01)) +
                scale_y_discrete(expand = expansion(mult = 0.01)) +
                scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c("TRUE" = 1, "FALSE" = 0), guide = FALSE) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_blank(), axis.text.y = element_text(size = axis.text.size, color = "black"),
                      axis.ticks.x = element_blank(), axis.ticks.y = element_line(size = 0.8, color = "black"),
                      axis.title = element_blank(), legend.background = element_blank(), legend.position = legend.position,
                      legend.text = element_text(size = legend.text.size), legend.title = element_text(size = legend.title.size),
                      panel.background = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(), plot.background = element_blank(), plot.margin = unit(c(1,6,1,1), "lines"))
        colColors <- ggplot(data = data.frame(x = 1:length(levels(MEtraitCor$module)), y = 0,
                                              color = levels(MEtraitCor$module))) +
                geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                scale_fill_identity(aesthetics = c("color", "fill")) +
                theme_void() +
                theme(legend.position = "none", plot.margin = unit(colColorMargins, "lines"))
        gg <- plot_grid(heatmap, colColors, nrow = 2, rel_heights = c(1, 0.045))
        if(save){
                if(verbose){
                        message("[plotMEtraitCor] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

plotMEtraitDot <- function(ME, trait, traitCode = NULL, colors = NULL, fun.data = "median_hilow",
                           fun.args = list(conf.int = 0.5), binwidth = 0.01, stackratio = 1.4, dotsize = 0.85, ylim = NULL,
                           nBreaks = 4, axis.title.size = 20, axis.text.size = 16, xlab = "Trait", ylab = "Module Eigennode",
                           save = TRUE, file = "ME_Trait_Dotplot.pdf", width = 6, height = 6, verbose = TRUE){
        if(verbose){
                message("[plotMEtraitDot] Plotting module eigennode by categorical trait")
        }
        if(sum(is.na(trait)) > 1){
                message("[plotMEtraitDot] Removing NA values")
                ME <- ME[!is.na(trait)]
                trait <- trait[!is.na(trait)]
        }
        if(!is.null(traitCode)){
                trait <- names(traitCode)[match(trait, traitCode)] %>% factor(levels = names(traitCode))
        } else {
                trait <- as.factor(trait)
        }
        dotplot <- ggplot() +
                stat_summary(aes(x = trait, y = ME, group = trait), fun.data = fun.data, geom = "crossbar",
                             color = "black", size = 0.5, fun.args = fun.args) +
                geom_dotplot(aes(x = trait, y = ME, fill = trait, color = trait), binwidth = binwidth, binaxis = "y",
                             stackdir = "center", position = "dodge", stackratio = stackratio, dotsize = dotsize) +
                coord_cartesian(ylim = ylim) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                theme_bw(base_size = 25) +
                theme(legend.position = "none", panel.grid.major = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25),
                      panel.grid.minor = element_blank(), strip.background = element_blank(),
                      axis.text = element_text(color = "black", size = axis.text.size),
                      axis.title = element_text(size = axis.title.size), plot.margin = unit(c(1,1,1,1), "lines")) +
                xlab(xlab) +
                ylab(ylab)
        if(!is.null(colors)){
                dotplot <- dotplot +
                        scale_color_manual(breaks = names(colors), values = colors, aesthetics = c("color", "fill"))
        }
        if(verbose){
                message("[plotMEtraitDot] Saving file as ", file)
        }
        ggsave(file, plot = dotplot, dpi = 600, width = width, height = height, units = "in")
}

plotMEtraitScatter <- function(ME, trait, color = "#132B43", xlim = NULL, ylim = NULL, nBreaks = 4, point.size = 2.5,
                               axis.title.size = 20, axis.text.size = 16, xlab = "Trait", ylab = "Module Eigennode",
                               save = TRUE, file = "ME_Trait_Scatterplot.pdf", width = 6, height = 6, verbose = TRUE){
        if(verbose){
                message("[plotMEtraitScatter] Plotting module eigennode by continuous trait")
        }
        scatterplot <- ggplot() +
                geom_smooth(aes(x = trait, y = ME), method = MASS::rlm, formula = y ~ x, color = "#56B1F7", fill = "#336A98") +
                geom_point(aes(x = trait, y = ME), color = color, size = point.size) +
                coord_cartesian(xlim = xlim, ylim = ylim) +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                theme_bw(base_size = 25) +
                theme(legend.position = "none", panel.grid.major = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25),
                      panel.grid.minor = element_blank(), strip.background = element_blank(),
                      axis.text = element_text(color = "black", size = axis.text.size),
                      axis.title = element_text(size = axis.title.size), plot.margin = unit(c(1,1,1,1), "lines")) +
                xlab(xlab) +
                ylab(ylab)
        if(verbose){
                message("[plotMEtraitScatter] Saving file as ", file)
        }
        ggsave(file, plot = scatterplot, dpi = 600, width = width, height = height, units = "in")
}

plotMethTrait <- function(module, regions, meth, trait, discrete = NULL, traitCode = NULL, traitColors = NULL,
                          heatmapColors = blueWhiteRed(100, gamma = 0.3), limit = NULL, expandY = 0.05, axis.text.size = 11,
                          heatmap.legend.position = c(1.1,0.743), trait.legend.position = c(1.017,4.39),
                          heatmap.legend.title = "Relative\nMethylation (%)", trait.legend.title = "Trait",
                          legend.text.size = 11, legend.title.size = 14, heatmapMargins = c(1,8,0,1),
                          traitMargins = c(0,6,1,5.15), save = TRUE, file = "Module_Methylation_Trait_Heatmap.pdf",
                          width = 11, height = 4, verbose = TRUE){
        if(verbose){
                message("[plotMethTrait] Plotting ", module, " module region methylation by ", trait.legend.title)
        }
        if(sum(is.na(trait)) > 1){
                message("[plotMethTrait] Removing NA values")
                meth <- meth[,!is.na(trait)]
                trait <- trait[!is.na(trait)]
        }
        RegionIDs <- rev(regions$RegionID[regions$module == module])
        if(is.null(discrete)){
                discrete <- length(unique(trait)) <= 5
        }
        if(discrete){
                if(!is.null(traitCode)){
                        trait <- names(traitCode)[match(trait, traitCode)] %>% factor(levels = names(traitCode))
                } else {
                        trait <- as.factor(trait)
                }
        }
        meth <- meth[RegionIDs,order(trait)]
        meth <- (meth - DelayedMatrixStats::rowMeans2(meth)) %>% reshape2::melt()
        if(is.null(limit)){
                limit <- max(abs(meth$value))
        }
        heatmap <- ggplot(data = meth) +
                geom_tile(aes(x = Var2, y = Var1, color = value, fill = value)) +
                scale_fill_gradientn(heatmap.legend.title, colors = heatmapColors, limits = c(-limit,limit),
                                     aesthetics = c("color", "fill")) +
                scale_y_discrete(expand = expansion(mult = expandY)) +
                theme_bw(base_size = 24) +
                theme(axis.text.x = element_blank(), axis.text.y = element_text(size = axis.text.size, color = "black"),
                      axis.ticks.x = element_blank(), axis.ticks.y = element_line(size = 0.8, color = "black"),
                      axis.title = element_blank(), legend.background = element_blank(),
                      legend.position = heatmap.legend.position, legend.text = element_text(size = legend.text.size),
                      legend.title = element_text(size = legend.title.size), panel.background = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25), panel.grid = element_blank(),
                      plot.background = element_blank(), plot.margin = unit(heatmapMargins, "lines"))
        colColors <- ggplot(data = data.frame(x = 1:length(trait), y = 0, color = sort(trait))) +
                geom_tile(aes(x = x, y = y, color = color, fill = color)) +
                theme_void() +
                theme(legend.position = trait.legend.position, legend.text = element_text(size = legend.text.size),
                      legend.title = element_text(size = legend.title.size), plot.margin = unit(traitMargins, "lines"))
        if(discrete){
                if(!is.null(traitColors)){
                        colColors <- colColors +
                                scale_color_manual(trait.legend.title, breaks = names(traitColors), values = traitColors,
                                                   aesthetics = c("color", "fill"))
                } else {
                        colColors <- colColors +
                                scale_color_discrete(trait.legend.title, aesthetics = c("color", "fill"))
                }
        } else {
                colColors <- colColors +
                        scale_color_continuous(trait.legend.title, aesthetics = c("color", "fill"))
        }
        gg <- plot_grid(heatmap, colColors, nrow = 2, rel_heights = c(1,0.15))
        if(save){
                if(verbose){
                        message("[plotMethTrait] Saving file as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
}

