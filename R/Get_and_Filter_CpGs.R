getCpGs <- function(colData, path = getwd(), pattern = "*CpG_report.txt.gz",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), BPPARAM = MulticoreParam(10),
                    save = TRUE, file = "Unfiltered_BSseq.rds", verbose = TRUE){
        if(verbose){
                message("[getCpGs] Loading CpG-level data")
        }
        files <- list.files(path, pattern = pattern) %>% .[pmatch(rownames(colData), table = .)]
        loci <- bsseq:::.readBismarkAsFWGRanges(files[1]) %>% chrSelectBSseq(seqnames = chroms)
        bs <- read.bismark(files, loci = loci, colData = colData, BPPARAM = BPPARAM, verbose = verbose)
        if(verbose){
                message("[getCpGs] Final BSseq Object:")
                print(bs)
        }
        if(save){
                if(verbose){
                        message("[getCpGs] Saving file as ", file)
                }
                saveRDS(bs, file = file)
        }
        return(bs)
}

getCpGtotals <- function(bs, cov = seq(0,10,1), perSample = seq(0.5,1,0.05), save = TRUE, file = "CpG_Totals.txt",
                         verbose = TRUE){
        if(verbose){
                message("[getCpGtotals] Calculating CpG totals at specified cov and perSample cutoffs")
        }
        nSample <- (perSample * ncol(bs)) %>% ceiling()
        bsCov <- getCoverage(bs)
        nCpGs <- sapply(cov, FUN = .nCpGsByCov, bsCov = bsCov, nSample = nSample)
        nCpGs <- as.integer(nCpGs)
        perCpGs <- (nCpGs / length(bs))
        CpGtotals <- data.frame(cov = rep(cov, each = length(perSample)), perSample = rep(perSample, times = length(cov)),
                                nSample = rep(nSample, times = length(cov)), nCpGs_M = nCpGs / 10^6, perCpGs = perCpGs)
        if(save){
                if(verbose){
                        message("[getCpGtotals] Saving file as ", file)
                }
                write.table(CpGtotals, file = file, sep = "\t", row.names = FALSE)
        }
        return(CpGtotals)
}

plotCpGtotals <- function(CpGtotals, nBreaks = 4, legend.position = c(1.08,0.73), save = TRUE,
                          file = "CpG_Totals.pdf", width = 11, height = 4.25, verbose = TRUE){
        if(verbose){
                message("[plotCpGtotals] Plotting CpG totals")
        }
        CpGtotals$perSample <- CpGtotals$perSample * 100
        gg <- ggplot(data = CpGtotals)
        gg <- gg +
                geom_line(aes(x = perSample, y = nCpGs_M, group = cov, color = cov)) +
                geom_text(data = subset(CpGtotals, perSample == min(perSample)),
                          aes(x = perSample, y = nCpGs_M, group = cov, color = cov, label = cov),
                          size = 4.5, check_overlap = TRUE, nudge_x = -0.5, hjust = 1) +
                xlab("Samples (%) Cutoff") +
                ylab("Total CpGs (Millions)") +
                scale_x_continuous(breaks = breaks_pretty(n = nBreaks + 1), expand = expansion(mult = c(0.05, 0.03))) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks)) +
                scale_color_gradient("Coverage\nCutoff", breaks = breaks_pretty(n = nBreaks - 1)) +
                theme_bw(base_size = 24) +
                theme(axis.text = element_text(size = 14, color = "black"),
                      axis.ticks = element_line(size = 1.25, color = "black"), axis.title = element_text(size = 18),
                      legend.background = element_blank(), legend.position = legend.position,
                      legend.title = element_text(size = 18), legend.text = element_text(size = 14),
                      panel.border = element_rect(color = "black", size = 1.25),
                      panel.grid = element_blank(), plot.margin = unit(c(1,7,0.7,0.7), "lines"))
        if(save){
                if(verbose){
                        message("[plotCpGtotals] Saving plot as ", file)
                }
                ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
        }
        return(gg)
}

filterCpGs <- function(bs, cov = 2, perSample = 0.75, save = TRUE, file = "Filtered_BSseq.rds", verbose = TRUE){
        if(verbose){
                message("[getCpGs] Filtering CpG-level data for loci with at least ", cov, " reads in at least ",
                        perSample * 100, "% of samples")
        }
        covSample <- (getCoverage(bs) >= cov) %>% DelayedMatrixStats::rowSums2()
        bs <- bs[covSample >= (perSample * ncol(bs)),]
        if(verbose){
                message("[getCpGs] Final BSseq Object:")
                print(bs)
        }
        if(save){
                if(verbose){
                        message("[getCpGs] Saving file as ", file)
                }
                saveRDS(bs, file = file)
        }
        return(bs)
}

.nCpGsBySample <- function(covSample, nSample){
        n <- length(covSample[covSample >= nSample])
        return(n)
}

.nCpGsByCov <- function(bsCov, cov, nSample){
        covSample <- (bsCov >= cov) %>% DelayedMatrixStats::rowSums2()
        nCpGs <- sapply(nSample, FUN = .nCpGsBySample, covSample = covSample)
        return(nCpGs)
}
