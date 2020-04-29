getSoftPower <- function(meth, powerVector = 1:20, corType = c("pearson", "bicor"), maxPOutliers = 0.1,
                         RsquaredCut = 0.8, blockSize = 40000, gcInterval = blockSize - 1, save = TRUE,
                         file = "Soft_Power.rds", verbose = TRUE){
        corType <- match.arg(corType)
        if(verbose){
                message("[getSoftPower] Analyzing scale-free topology with ", corType,
                        " correlation to determine best soft-thresholding power")
                verboseNum <- 1
        } else {
                verboseNum <- 0
        }
        if(corType == "pearson"){
                sft <- pickSoftThreshold(meth, RsquaredCut = RsquaredCut, powerVector = powerVector, networkType = "signed",
                                         moreNetworkConcepts = TRUE, corFnc = "cor", blockSize = blockSize,
                                         gcInterval = gcInterval, verbose = verboseNum)
        } else {
                if(corType == "bicor"){
                        sft <- pickSoftThreshold(meth, RsquaredCut = RsquaredCut, powerVector = powerVector,
                                                 networkType = "signed", moreNetworkConcepts = TRUE, corFnc = "bicor",
                                                 corOptions = list(maxPOutliers = maxPOutliers), blockSize = blockSize,
                                                 gcInterval = gcInterval, verbose = verboseNum)
                } else {
                        stop("[getSoftPower] Error: corType must be either pearson or bicor")
                }
        }
        if(is.na(sft$powerEstimate)){
                fit <- -sign(sft$fitIndices[,"slope"]) * sft$fitIndices[,"SFT.R.sq"]
                sft$powerEstimate <- sft$fitIndices$Power[fit == max(fit)]
        }
        if(verbose){
                message("[getSoftPower] At soft power threshold = ", sft$powerEstimate,
                        ", fit = ", round(sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == sft$powerEstimate], 3),
                        " and mean connectivity = ", round(sft$fitIndices$mean.k.[sft$fitIndices$Power == sft$powerEstimate], 1))
        }
        if(save){
                if(verbose){
                        message("[getSoftPower] Saving file as ", file)
                }
                saveRDS(sft, file = file)
        }
        return(sft)
}

plotSoftPower <- function(sft, pointCol = "#132B43", lineCol = "red", nBreaks = 4, save = TRUE,
                          file = "Soft_Power_Plots.pdf", width = 8.5, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotSoftPower] Plotting scale-free topology fit and mean connectivity by soft power threshold")
        }
        fitIndices <- data.frame(power = sft$fitIndices$Power,
                                 fit = -sign(sft$fitIndices[,"slope"]) * sft$fitIndices[,"SFT.R.sq"],
                                 log10_meanConnectivity = log10(sft$fitIndices$mean.k.),
                                 powerEstimate = sft$powerEstimate) %>%
                reshape2::melt(id.vars = c("power", "powerEstimate"))
        powerEstimateY <- min(0, fitIndices$value)
        gg <- ggplot(data = fitIndices)
        gg <- gg +
                geom_vline(aes(xintercept = powerEstimate), color = lineCol) +
                geom_text(aes(x = powerEstimate, y = powerEstimateY, label = powerEstimate), color = lineCol,
                          nudge_x = -1) +
                geom_point(aes(x = power, y = value), color = pointCol, size = 1.2) +
                facet_wrap(vars(variable), nrow = 1, ncol = 2, scales = "free_y", strip.position = "left") +
                xlab("Soft Power Threshold") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                expand_limits(x = 0, y = c(0,1)) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 12, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"), axis.title.x = element_text(size = 16),
                      axis.title.y = element_blank(), legend.position = "none",
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(), panel.spacing.x = unit(0.3, "lines"),
                      panel.spacing.y = unit(0.8, "lines"), plot.margin = unit(c(1,1,0.7,0.2), "lines"),
                      strip.background = element_blank(), strip.placement = "outside",
                      strip.switch.pad.wrap = unit(0, "lines"), strip.text.x = element_text(size = 16))
        if(save){
                if(verbose){
                        message("[plotSoftPower] Saving plots as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

getModules <- function(meth, power, regions, maxBlockSize = 40000, corType = c("pearson", "bicor"),
                       maxPOutliers = 0.1, deepSplit = 4, minModuleSize = 10, mergeCutHeight = 0.1, nThreads = 4, save = TRUE,
                       file = "Modules.rds", verbose = TRUE){
        if(is.null(power)){
                stop("[getModules] You must select a soft power threshold")
        }
        if(is.null(regions)){
                stop("[getModules] You must include regions")
        }
        corType <- match.arg(corType)
        if(!corType %in% c("pearson", "bicor")){
                stop("[getModules] corType must be either pearson or bicor")
        }
        if(verbose){
                message("[getModules] Constructing network and detecting modules in blocks using ", corType, " correlation")
                verboseNum <- 3
        } else {
                verboseNum <- 0
        }
        modules <- blockwiseModules(meth, checkMissingData = FALSE, maxBlockSize = maxBlockSize, corType = corType,
                                    maxPOutliers = maxPOutliers, power = power, networkType = "signed", TOMtype = "signed",
                                    deepSplit = deepSplit, minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight,
                                    nThreads = nThreads, verbose = verboseNum)
        if(verbose){
                message("[getModules] Assigning modules and calculating module membership using ", corType, " correlation")
        }
        colnames(modules$MEs) <- str_remove_all(colnames(modules$MEs), pattern = "ME")
        if(corType == "pearson"){
                membership <- WGCNA::cor(x = meth, y = modules$MEs, use = "pairwise.complete.obs", nThreads = nThreads)
        } else {
                membership <- bicor(x = meth, y = modules$MEs, use = "pairwise.complete.obs", maxPOutliers = maxPOutliers,
                                    nThreads = nThreads)
        }
        regions$module <- modules$colors[match(regions$RegionID, names(modules$colors))]
        regions <- lapply(unique(regions$module), function(x){
                regions <- regions[regions$module == x,]
                regions$membership <- membership[regions$RegionID,x]
                regions$hubRegion <- regions$membership == max(regions$membership)
                return(regions)
        })
        regions <- list.rbind(regions) %>% .[order(as.integer(str_remove_all(.$RegionID, pattern = "Region_"))),]
        modules$regions <- regions
        if(save){
                if(verbose){
                        message("[getModules] Saving modules as ", file)
                }
                saveRDS(modules, file = file)
        }
        return(modules)
}

