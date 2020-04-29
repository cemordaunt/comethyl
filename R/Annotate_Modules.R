annotateModule <- function(regions, module = NULL, grey = FALSE, genome = c("hg38", "hg19", "hg18", "mm10", "mm9", "danRer7"),
                           includeCuratedRegDoms = FALSE, rule = c("basalPlusExt", "twoClosest", "oneClosest"),
                           adv_upstream = 5, adv_downstream = 1, adv_span = 1000, adv_twoDistance = 1000,
                           adv_oneDistance = 1000, version = c("4.0.4", "3.0.0", "2.0.2"), save = TRUE,
                           file = "Annotated_Module_Regions.txt", verbose =  TRUE){
        if(!is.null(module)){
                if(verbose){
                        message("[annotateModule] Filtering for regions in ", paste(module, collapse = ", "), " module(s)")
                }
                regions <- regions[regions$module %in% module,]
        }
        if(!grey){
                if(verbose){
                        message("[annotateModule] Excluding regions in grey (unassigned) module")
                }
                regions <- regions[!regions$module == "grey",]
        }
        genome <- match.arg(genome)
        rule <- match.arg(rule)
        version <- match.arg(version)
        if((version == "4.0.4" & genome %in% c("hg18", "danRer7")) | (version == "3.0.0" & genome %in% c("hg38", "hg18")) |
           (version == "2.0.2" & genome %in% c("hg38", "mm10"))){
                stop("[annotateModule] The ", genome, " genome assembly is not supported for GREAT v", version)
        }
        if(verbose){
                message("[annotateModule] Using the ", genome, " genome assembly for annotations")
                message("[annotateModule] Adding genes to regions using GREAT v", version, " with the ", rule, " rule")
        }
        GR_regions <- with(regions, GRanges(chr, ranges = IRanges(start, end = end), RegionID = RegionID))
        job <- submitGreatJob(GR_regions, species = genome, includeCuratedRegDoms = includeCuratedRegDoms, rule = rule,
                              adv_upstream = adv_upstream, adv_downstream = adv_downstream, adv_span = adv_span,
                              adv_twoDistance = adv_twoDistance, adv_oneDistance = adv_oneDistance, request_interval = 0,
                              version = version)
        regions_genes <- suppressGraphics(plotRegionGeneAssociationGraphs(job, type = 1, request_interval = 0)) %>%
                as.data.frame() %>%
                merge(x = regions[,c("RegionID", "chr", "start", "end")],
                      y = .[,c("seqnames", "start", "end", "gene", "distTSS")], by.x = c("chr", "start", "end"),
                      by.y = c("seqnames", "start", "end"), all = TRUE, sort = FALSE)
        regions_genes$RegionID <- factor(regions_genes$RegionID, levels = unique(regions_genes$RegionID))
        colnames(regions_genes) <- str_replace_all(colnames(regions_genes), pattern = c("gene" = "gene_symbol",
                                                                                        "distTSS" = "distance_to_TSS"))
        if(verbose){
                message("[annotateModule] Adding gene info from BioMart")
        }
        dataset <- str_sub(genome, start = 1, end = 2) %>% switch(hg = "hsapiens_gene_ensembl", mm = "mmusculus_gene_ensembl",
                                                                  da = "drerio_gene_ensembl")
        ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset)
        genes_annotated <- suppressMessages(getBM(attributes = c("external_gene_name", "description", "ensembl_gene_id", "entrezgene_id"),
                                                  filters = c("external_gene_name"), values = regions_genes$gene, mart = ensembl,
                                                  useCache = FALSE))
        colnames(genes_annotated) <- colnames(genes_annotated) %>%
                str_replace_all(pattern = c("external_gene_name" = "gene_symbol", "description" = "gene_description",
                                            "ensembl_gene_id" = "gene_ensemblID", "entrezgene_id" = "gene_entrezID"))
        genes_annotated$gene_description <- str_split_fixed(genes_annotated$gene_description,
                                                            pattern = fixed(" ["), n = 2)[,1] %>%
                str_remove_all(",")
        regions_annotated <- merge(x = regions_genes, y = genes_annotated, by = "gene_symbol", all.x = TRUE, all.y = FALSE,
                                   sort = FALSE) %>%
                unique() %>%
                aggregate(formula = cbind(gene_symbol, distance_to_TSS, gene_description, gene_ensemblID, gene_entrezID) ~ RegionID,
                          data = ., FUN = function(x) paste(x, collapse = " | "), simplify = TRUE, drop = FALSE) %>%
                merge(x = regions, y = ., by = "RegionID", all.x = TRUE, all.y = FALSE, sort = FALSE)
        if(!genome %in% c("hg18", "danRer7")){
                if(verbose){
                        message("[annotateModule] Getting gene context from annotatr")
                }
                annotations <- paste(genome, c("basicgenes", "genes_intergenic", "enhancers_fantom"), sep = "_")
                regions_GeneReg <- suppressWarnings(suppressMessages(build_annotations(genome = genome, annotations = annotations))) %>%
                        GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
                        annotate_regions(regions = GR_regions, annotations = ., ignore.strand = TRUE, quiet = TRUE) %>%
                        as.data.frame()
                colnames(regions_GeneReg)[colnames(regions_GeneReg) == "annot.type"] <- "gene_context"
                pattern <- rep("", times = 5)
                names(pattern) <- c(genome, "genes", "s", "fantom", "_")
                regions_GeneReg$gene_context <- str_replace_all(regions_GeneReg$gene_context, pattern = pattern)
                regions_annotated <- aggregate(formula = gene_context ~ RegionID, data = regions_GeneReg,
                                               FUN = function(x) paste(unique(x), collapse = ", "), simplify = TRUE) %>%
                        merge(x = regions_annotated, y = ., by = "RegionID", all.x = TRUE, all.y = FALSE, sort = FALSE)
                if(verbose){
                        message("[annotateModule] Getting CpG context from annotatr")
                }
                regions_CpGs <- suppressMessages(build_annotations(genome = genome, annotations = paste(genome, "cpgs", sep = "_"))) %>%
                        GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
                        annotate_regions(regions = GR_regions, annotations = ., ignore.strand = TRUE, quiet = TRUE) %>%
                        as.data.frame()
                colnames(regions_CpGs)[colnames(regions_CpGs) == "annot.type"] <- "CpG_context"
                pattern <- c("island", "shore", "shelf", "open sea")
                names(pattern) <- paste(genome, "cpg", c("islands", "shores", "shelves", "inter"), sep = "_")
                regions_CpGs$CpG_context <- str_replace_all(regions_CpGs$CpG_context, pattern = pattern)
                regions_annotated <- aggregate(formula = CpG_context ~ RegionID, data = regions_CpGs,
                                               FUN = function(x) paste(unique(x), collapse = ", "), simplify = TRUE) %>%
                        merge(x = regions_annotated, y = ., by = "RegionID", all.x = TRUE, all.y = FALSE, sort = FALSE)
        } else {
                if(verbose){
                        message("[annotateModule] Gene and CpG context not supported by annotatr for hg18 and danRer7")
                }
        }
        regions_annotated <- regions_annotated[order(regions_annotated$module, as.integer(str_remove_all(regions_annotated$RegionID,
                                                                                                         pattern = "Region_"))),]
        if(save){
                if(verbose){
                        message("[annotateModule] Saving file as ", file)
                }
                write.table(regions_annotated, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
        }
        return(regions_annotated)
}

getGeneList <- function(regions_annotated, module = NULL,type = c("symbol", "description", "ensemblID", "entrezID"),
                        verbose = TRUE){
        if(!is.null(module)){
                if(verbose){
                        message("[getGeneList] Filtering for regions in ", paste(module, collapse = ", "), " module(s)")
                }
                regions_annotated <- regions_annotated[regions_annotated$module %in% module,]
        }
        type <- match.arg(type)
        if(verbose){
                message("[getGeneList] Getting gene ", type, "s")
        }
        geneList <- str_split(regions_annotated[,paste("gene", type, sep = "_")], pattern = fixed(" | ")) %>%
                as_vector() %>% unique() %>% sort()
        return(geneList)
}

