listOntologies <- function(genome = c("hg38", "hg19", "hg18", "mm10", "mm9", "danRer7"),
                           version = c("4.0.4", "3.0.0", "2.0.2"), verbose =  TRUE){
        genome <- match.arg(genome)
        version <- match.arg(version)
        if((version == "4.0.4" & genome %in% c("hg18", "danRer7")) | (version == "3.0.0" & genome %in% c("hg38", "hg18")) |
           (version == "2.0.2" & genome %in% c("hg38", "mm10"))){
                stop("[listOntologies] The ", genome, " genome assembly is not supported for GREAT v", version)
        }
        if(verbose){
                message("[listOntologies] Getting available ontologies for GREAT v", version, " with the ", genome,
                        " genome assembly")
        }
        ontologies <- GRanges("chr1", ranges = IRanges(1, end = 2)) %>%
                submitGreatJob(species = genome, request_interval = 0, version = version) %>%
                availableOntologies()
        return(ontologies)
}

enrichModule <- function(regions, module = NULL, genome = c("hg38", "hg19", "hg18", "mm10", "mm9", "danRer7"),
                         includeCuratedRegDoms = FALSE, rule = c("basalPlusExt", "twoClosest", "oneClosest"),
                         adv_upstream = 5, adv_downstream = 1, adv_span = 1000, adv_twoDistance = 1000,
                         adv_oneDistance = 1000, version = c("4.0.4", "3.0.0", "2.0.2"),
                         ontologies = c("GO Molecular Function", "GO Biological Process", "GO Cellular Component",
                                        "Mouse Phenotype", "Human Phenotype"), min_background_region_hits = 5,
                         adjMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
                         min_region_hits = 2, pvalue_threshold = 0.01, save = TRUE, file = "Module_Enrichment_Results.txt",
                         verbose =  TRUE){
        if(!is.null(module)){
                if(verbose){
                        message("[enrichModule] Analyzing functional enrichment for regions in ",
                                paste(module, collapse = ", "), " module(s)")
                }
                moduleRegions <- regions[regions$module %in% module,]
        } else {
                if(verbose){
                        message("[enrichModule] Analyzing functional enrichment for all regions assigned to a module")
                }
                moduleRegions <- regions[!regions$module == "grey",]
        }
        genome <- match.arg(genome)
        rule <- match.arg(rule)
        version <- match.arg(version)
        if((version == "4.0.4" & genome %in% c("hg18", "danRer7")) | (version == "3.0.0" & genome %in% c("hg38", "hg18")) |
           (version == "2.0.2" & genome %in% c("hg38", "mm10"))){
                stop("[enrichModule] The ", genome, " genome assembly is not supported for GREAT v", version)
        }
        adjMethod <- match.arg(adjMethod)
        if(verbose){
                message("[enrichModule] Using GREAT v", version, " with the ", genome, " genome assembly and the ", rule,
                        " rule")
        }
        if(genome == "danRer7"){
                ontologies <- ontologies[!ontologies %in% c("Mouse Phenotype", "Human Phenotype")]
        }
        gr <- with(moduleRegions, GRanges(chr, ranges = IRanges(start, end = end)))
        bg <- with(regions, GRanges(chr, ranges = IRanges(start, end = end)))
        job <- submitGreatJob(gr, bg = bg, species = genome, includeCuratedRegDoms = includeCuratedRegDoms, rule = rule,
                              adv_upstream = adv_upstream, adv_downstream = adv_downstream, adv_span = adv_span,
                              adv_twoDistance = adv_twoDistance, adv_oneDistance = adv_oneDistance, request_interval = 0,
                              version = version)
        if(verbose){
                message("[enrichModule] Getting results and adjusting p-values using the ", adjMethod, " method")
        }
        enrichTables <- getEnrichmentTables(job, ontology = ontologies)
        results <- list.rbind(enrichTables)
        rownames(results) <- 1:nrow(results)
        pattern <- c("hyper_" = "", "id" = "ID", "name" = "term", "total_regions" = "background_region_hits",
                     "expected" = "expected_region_hits", "foreground_region_hits" = "region_hits",
                     "region_set_coverage" = "region_coverage", "foreground_gene_hits" = "gene_hits",
                     "total_genes_annotated" = "term_genes", "raw_pvalue" = "p")
        colnames(results) <- colnames(results) %>% str_to_lower() %>% str_replace_all(pattern = pattern)
        results$ontology <- rep(names(enrichTables), sapply(enrichTables, nrow))
        results <- subset(results, background_region_hits >= min_background_region_hits)
        results$log_p <- -log10(results$p)
        results$adj_p <- p.adjust(results$p, method = adjMethod)
        results$log_adj_p <- -log10(results$adj_p)
        results <- subset(results, region_hits >= min_region_hits & p < pvalue_threshold)
        if(verbose){
                message("[enrichModule] Extracting genes for enriched terms")
        }
        results$genes <- mapply(function(x,y){
                suppressGraphics(plotRegionGeneAssociationGraphs(job, type = 1, ontology = x, termID = y, request_interval = 0,
                                                                 max_tries = 10000, verbose = FALSE)) %>% .$gene %>%
                        unique() %>% sort() %>% paste(collapse = ", ")
        }, x = results$ontology, y = results$ID, USE.NAMES = FALSE)
        results <- results[order(results$p), c("ID", "term", "ontology", "background_region_hits", "expected_region_hits",
                                               "region_hits", "fold_enrichment", "region_coverage", "term_region_coverage",
                                               "gene_hits", "genes", "background_gene_hits", "term_genes", "p", "log_p",
                                               "adj_p", "log_adj_p")]
        if(save){
                if(verbose){
                        message("[enrichModule] Saving file as ", file)
                }
                write.table(results, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
        }
        return(results)
}

plotEnrichment <- function(enrichment, nTerms = 15, fill = "#132B43", xlim = NULL, nBreaks = 4, axis.title.x.size = 20,
                           axis.text.x.size = 16, axis.text.y.size = 16, save = TRUE, file = "Module_Enrichment_Plot.pdf",
                           width = 8, height = 6, verbose = TRUE){
        if(verbose){
                message("[plotEnrichment] Plotting module enrichments from GREAT")
        }
        enrichment <- enrichment[order(enrichment$log_p, decreasing = TRUE),]
        if(nrow(enrichment) > nTerms){
                enrichment <- enrichment[1:nTerms,]
        }
        enrichment$term <- factor(enrichment$term, levels = rev(unique(enrichment$term)))
        scatterplot <- ggplot() +
                geom_col(aes(x = term, y = log_p), data = enrichment, fill = fill) +
                coord_flip(ylim = xlim) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks), expand = expansion(c(0.004, 0.03))) +
                theme_bw(base_size = 25) +
                theme(legend.position = "none", panel.grid.major = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25),
                      panel.grid.minor = element_blank(), strip.background = element_blank(),
                      axis.text.x = element_text(color = "black", size = axis.text.x.size),
                      axis.text.y = element_text(color = "black", size = axis.text.y.size),
                      axis.title.x = element_text(size = axis.title.x.size), axis.title.y = element_blank(),
                      plot.margin = unit(c(1,1,0.5,1), "lines")) +
                ylab(expression(-log[10]*(italic(p)-value)))
        if(verbose){
                message("[plotEnrichment] Saving file as ", file)
        }
        ggsave(file, plot = scatterplot, dpi = 600, width = width, height = height, units = "in")
}

