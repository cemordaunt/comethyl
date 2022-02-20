#' Annotate Module Regions
#'
#' \code{annotateModule()} takes a \code{data.frame} of regions with module
#' assignments, annotates them with genes using GREAT, adds additional gene
#' information with Ensembl BioMart, provides gene regulatory context with
#' annotatr, and then saves this as a .txt file. Support is provided for several
#' genomes, including \code{hg38}, \code{hg19}, \code{hg18}, \code{mm10},
#' \code{mm9}, and \code{danRer7}.
#'
#' \code{regions} can be filtered for one or more modules of interest, or the
#' grey (unassigned) module can be excluded. Gene annotation is performed by the
#' \pkg{rGREAT} package, which allows for different annotation rules and versions
#' of GREAT. The default \code{basalPlusExt} annotation rule associates a gene
#' with a region if the region is within the basal regulatory domain of the gene
#' (5 kb upstream and 1 kb downstream of the TSS) or if it is within 1 Mb
#' upstream or downstream of the TSS and not in the basal regulatory domain of
#' another gene. Other rules include \code{twoClosest} and \code{oneClosest},
#' which effectively assign the two nearest genes or one nearest genes,
#' respectively. See [rGREAT::submitGreatJob()] for more details.
#'
#' GREAT supports different genomes depending on the version:
#' \describe{
#'         \item{GREAT v4.0.4}{\code{hg38}, \code{hg19}, \code{mm10}, \code{mm9}}
#'         \item{GREAT v3.0.0}{\code{hg19}, \code{mm10}, \code{mm9}, \code{danRer7}}
#'         \item{GREAT v2.0.2}{\code{hg19}, \code{hg18}, \code{mm9}, \code{danRer7}}
#' }
#'
#' Gene information is provided by the \pkg{biomaRt} package, which adds the gene
#' description along with Ensembl and NCBI Entrez gene IDs. Regulatory context
#' is added by the \pkg{annotatr} package. This provides positional context of
#' the region relative to nearby genes, enhancers, and CpG islands. Note that
#' \pkg{annotatr} does not support the \code{hg18} or \code{danRer7} genomes.
#'
#' @param regions A \code{data.frame} of regions with module assignments,
#'         typically obtained from [getModules()].
#' @param module A \code{character} giving the name of one or more modules to
#'         annotate. If null, all modules will be annotated.
#' @param grey A \code{logical(1)} specifying whether or not to include the grey
#'         (unassigned) module.
#' @param genome A \code{character(1)} giving the genome build for the regions.
#'         Supported genomes include \code{hg38}, \code{hg19}, \code{hg18},
#'         \code{mm10}, \code{mm9}, and \code{danRer7}. See \code{Details}.
#' @param includeCuratedRegDoms A \code{logical(1)} indicating whether to include
#'         curated regulatory domains for GREAT gene annotation.
#' @param rule A \code{character(1)} specifying the rule used by GREAT for gene
#'         annotation. Possible values include \code{basalPlusExt},
#'         \code{twoClosest}, and \code{oneClosest}. See [rGREAT::submitGreatJob()]
#'         for more details for this and the next six arguments.
#' @param adv_upstream A \code{numeric(1)} giving the distance upstream of the
#'         TSS (in kb) to define a basal regulatory domain in the
#'         \code{basalPlusExt} rule.
#' @param adv_downstream A \code{numeric(1)} giving the distance downstream of
#'         the TSS (in kb) to define a basal regulatory domain in the
#'         \code{basalPlusExt} rule.
#' @param adv_span A \code{numeric(1)} specifying the distance upstream and
#'         downstream of the TSS (in kb) to define the maximum extension of the
#'         regulatory domain in the \code{basalPlusExt} rule.
#' @param adv_twoDistance A \code{numeric(1)} giving the distance upstream and
#'         downstream of the TSS (in kb) to define the maximum extension of the
#'         regulatory domain in the \code{twoClosest} rule.
#' @param adv_oneDistance A \code{numeric(1)} giving the distance upstream and
#'         downstream of the TSS (in kb) to define the maximum extension of the
#'         regulatory domain in the \code{oneClosest} rule.
#' @param version A \code{character(1)} specifying the version of GREAT to use
#'         for gene annotation. Possible values include \code{4.0.4},
#'         \code{3.0.0}, and \code{2.0.2}. Different versions of GREAT support
#'         different genomes. See \code{Details}.
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{data.frame}.
#' @param file A \code{character(1)} giving the file name (.txt) for the saved
#'         \code{data.frame}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{data.frame} adding gene and regulatory annotations to the
#'         regions.
#'
#' @seealso \itemize{
#'         \item [getModules()] to build a comethylation network and identify
#'                 modules of comethylated regions.
#'         \item [getGeneList()] to extract a list of genes or IDs from the
#'                 annotated regions.
#'         \item [listOntologies()], [enrichModule()], and [plotEnrichment()] to
#'                 investigate functional enrichment of module regions with GREAT.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Annotate Modules
#' regionsAnno <- annotateModule(regions, module = "bisque4",
#'                               genome = "hg38",
#'                               file = "Annotated_bisque4_Module_Regions.txt")
#' geneList_bisque4 <- getGeneList(regionsAnno, module = "bisque4")
#'
#' # Analyze Functional Enrichment
#' ontologies <- listOntologies("hg38", version = "4.0.4")
#' enrich_bisque4 <- enrichModule(regions, module = "bisque4",
#'                                genome = "hg38",
#'                                file = "bisque4_Module_Enrichment.txt")
#' plotEnrichment(enrich_bisque4,
#'                file = "bisque4_Module_Enrichment_Plot.pdf")
#' }
#'
#' @export
#'
#' @import GenomicRanges
#' @import rGREAT
#' @import stringr
#' @import biomaRt
#' @import utils
#'
#' @importFrom magrittr %>%

annotateModule <- function(regions, module = NULL, grey = FALSE,
                           genome = c("hg38", "hg19", "hg18", "mm10", "mm9",
                                      "danRer7"),
                           includeCuratedRegDoms = FALSE,
                           rule = c("basalPlusExt", "twoClosest", "oneClosest"),
                           adv_upstream = 5, adv_downstream = 1, adv_span = 1000,
                           adv_twoDistance = 1000, adv_oneDistance = 1000,
                           version = c("4.0.4", "3.0.0", "2.0.2"), save = TRUE,
                           file = "Annotated_Module_Regions.txt", verbose = TRUE){
        if(!is.null(module)){
                if(verbose){
                        message("[annotateModule] Filtering for regions in ",
                                paste(module, collapse = ", "), " module(s)")
                }
                regions <- regions[regions$module %in% module,]
        }
        if("grey" %in% regions$module & !grey){
                if(verbose){
                        message("[annotateModule] Excluding regions in grey (unassigned) module")
                }
                regions <- regions[!regions$module == "grey",]
        }
        genome <- match.arg(genome)
        rule <- match.arg(rule)
        version <- match.arg(version)
        if((version == "4.0.4" & genome %in% c("hg18", "danRer7")) |
           (version == "3.0.0" & genome %in% c("hg38", "hg18")) |
           (version == "2.0.2" & genome %in% c("hg38", "mm10"))){
                stop("[annotateModule] The ", genome,
                     " genome assembly is not supported for GREAT v", version)
        }
        if(verbose){
                message("[annotateModule] Using the ", genome,
                        " genome assembly for annotations")
                message("[annotateModule] Adding genes to regions using GREAT v",
                        version, " with the ", rule, " rule")
        }
        GR_regions <- with(regions, GRanges(chr,
                                            ranges = IRanges::IRanges(start, end = end),
                                            RegionID = RegionID))
        job <- submitGreatJob(GR_regions, species = genome,
                              includeCuratedRegDoms = includeCuratedRegDoms,
                              rule = rule, adv_upstream = adv_upstream,
                              adv_downstream = adv_downstream,
                              adv_span = adv_span,
                              adv_twoDistance = adv_twoDistance,
                              adv_oneDistance = adv_oneDistance,
                              request_interval = 0, version = version)
        regions_genes <- R.devices::suppressGraphics(
                plotRegionGeneAssociationGraphs(job, type = 1, request_interval = 0)) %>%
                as.data.frame() %>%
                merge(x = regions[,c("RegionID", "chr", "start", "end")],
                      y = .[,c("seqnames", "start", "end", "gene", "distTSS")],
                      by.x = c("chr", "start", "end"),
                      by.y = c("seqnames", "start", "end"), all = TRUE,
                      sort = FALSE)
        regions_genes$RegionID <- factor(regions_genes$RegionID,
                                         levels = unique(regions_genes$RegionID))
        colnames(regions_genes) <- str_replace_all(colnames(regions_genes),
                                                   pattern = c("gene" = "gene_symbol",
                                                               "distTSS" = "distance_to_TSS"))
        if(verbose){
                message("[annotateModule] Adding gene info from BioMart")
        }
        dataset <- str_sub(genome, start = 1, end = 2) %>%
                switch(hg = "hsapiens_gene_ensembl", mm = "mmusculus_gene_ensembl",
                       da = "drerio_gene_ensembl")
        ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset)
        genes_annotated <- suppressMessages(
                getBM(attributes = c("external_gene_name", "description",
                                     "ensembl_gene_id", "entrezgene_id"),
                      filters = c("external_gene_name"),
                      values = regions_genes$gene,
                      mart = ensembl, useCache = FALSE))
        colnames(genes_annotated) <- colnames(genes_annotated) %>%
                str_replace_all(pattern = c("external_gene_name" = "gene_symbol",
                                            "description" = "gene_description",
                                            "ensembl_gene_id" = "gene_ensemblID",
                                            "entrezgene_id" = "gene_entrezID"))
        genes_annotated$gene_description <- str_split_fixed(genes_annotated$gene_description,
                                                            pattern = fixed(" ["), n = 2)[,1] %>%
                str_remove_all(",")
        regions_annotated <- merge(x = regions_genes, y = genes_annotated,
                                   by = "gene_symbol", all.x = TRUE,
                                   all.y = FALSE, sort = FALSE) %>%
                unique() %>%
                stats::aggregate(formula = cbind(gene_symbol, distance_to_TSS,
                                                 gene_description, gene_ensemblID,
                                                 gene_entrezID) ~ RegionID,
                                 data = .,
                                 FUN = function(x) paste(x, collapse = " | "),
                                 simplify = TRUE, drop = FALSE) %>%
                merge(x = regions, y = ., by = "RegionID", all.x = TRUE,
                      all.y = FALSE, sort = FALSE)
        if(!genome %in% c("hg18", "danRer7")){
                if(verbose){
                        message("[annotateModule] Getting gene context from annotatr")
                }
                annotations <- paste(genome,
                                     c("basicgenes", "genes_intergenic", "enhancers_fantom"),
                                     sep = "_")
                regions_GeneReg <- suppressWarnings(suppressMessages(
                        annotatr::build_annotations(genome = genome,
                                                    annotations = annotations))) %>%
                        GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
                        annotatr::annotate_regions(regions = GR_regions,
                                                   annotations = .,
                                                   ignore.strand = TRUE,
                                                   quiet = TRUE) %>%
                        as.data.frame()
                colnames(regions_GeneReg)[colnames(regions_GeneReg) == "annot.type"] <- "gene_context"
                pattern <- rep("", times = 5)
                names(pattern) <- c(genome, "genes", "s", "fantom", "_")
                regions_GeneReg$gene_context <- str_replace_all(regions_GeneReg$gene_context,
                                                                pattern = pattern)
                regions_annotated <- stats::aggregate(formula = gene_context ~ RegionID,
                                                      data = regions_GeneReg,
                                                      FUN = function(x) paste(unique(x),
                                                                              collapse = ", "),
                                                      simplify = TRUE) %>%
                        merge(x = regions_annotated, y = ., by = "RegionID",
                              all.x = TRUE, all.y = FALSE, sort = FALSE)
                if(verbose){
                        message("[annotateModule] Getting CpG context from annotatr")
                }
                regions_CpGs <- suppressMessages(
                        annotatr::build_annotations(genome = genome,
                                                    annotations = paste(genome, "cpgs",
                                                                        sep = "_"))) %>%
                        GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
                        annotatr::annotate_regions(regions = GR_regions,
                                                   annotations = .,
                                                   ignore.strand = TRUE,
                                                   quiet = TRUE) %>%
                        as.data.frame()
                colnames(regions_CpGs)[colnames(regions_CpGs) == "annot.type"] <- "CpG_context"
                pattern <- c("island", "shore", "shelf", "open sea")
                names(pattern) <- paste(genome, "cpg", c("islands", "shores",
                                                         "shelves", "inter"),
                                        sep = "_")
                regions_CpGs$CpG_context <- str_replace_all(regions_CpGs$CpG_context,
                                                            pattern = pattern)
                regions_annotated <- stats::aggregate(formula = CpG_context ~ RegionID,
                                                      data = regions_CpGs,
                                                      FUN = function(x) paste(unique(x),
                                                                              collapse = ", "),
                                                      simplify = TRUE) %>%
                        merge(x = regions_annotated, y = ., by = "RegionID",
                              all.x = TRUE, all.y = FALSE, sort = FALSE)
        } else {
                if(verbose){
                        message("[annotateModule] Gene and CpG context not supported by annotatr for hg18 and danRer7")
                }
        }
        regions_annotated <- regions_annotated[order(regions_annotated$module,
                                                     as.integer(str_remove_all(regions_annotated$RegionID,
                                                                               pattern = "Region_"))),]
        if(save){
                if(verbose){
                        message("[annotateModule] Saving file as ", file)
                }
                write.table(regions_annotated, file = file, sep = "\t",
                            quote = FALSE, row.names = FALSE)
        }
        return(regions_annotated)
}

#' Extract a Gene List from Annotated Regions
#'
#' \code{getGeneList()} takes a \code{data.frame} of regions annotated with gene
#' information and extracts a \code{vector} of unique gene symbols, descriptions,
#' or IDs.
#'
#' \code{getGeneList()} is designed to be used in combination with
#' [annotateModule()]. \code{regions} can be filtered for one or more modules of
#' interest. Values that can be extracted include gene \code{symbol},
#' \code{description}, \code{ensemblID} and \code{entrezID}.
#'
#' @param regions_annotated A \code{data.frame} of regions with gene annotations,
#'         typically produced by [annotateModule()].
#' @param module A \code{character} giving the name of one or more modules to
#'         include. If null, all modules will be included.
#' @param type A \code{character(1)} specifying the type of gene information to
#'         extract. Possible values include \code{symbol}, \code{description},
#'         \code{ensemblID} and \code{entrezID}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{vector} of unique values.
#'
#' @seealso \itemize{
#'         \item [getModules()] to build a comethylation network and identify
#'                 modules of comethylated regions.
#'         \item [annotateModule()] to annotate a set of regions with genes and
#'                 regulatory context.
#'         \item [listOntologies()], [enrichModule()], and [plotEnrichment()] to
#'                 investigate functional enrichment of module regions with GREAT.
#' }
#'
#' @examples \dontrun{
#'
#' # Get Comethylation Modules
#' modules <- getModules(methAdj, power = sft$powerEstimate, regions = regions,
#'                       corType = "pearson", file = "Modules.rds")
#'
#' # Annotate Modules
#' regionsAnno <- annotateModule(regions, module = "bisque4",
#'                               genome = "hg38",
#'                               file = "Annotated_bisque4_Module_Regions.txt")
#' geneList_bisque4 <- getGeneList(regionsAnno, module = "bisque4")
#'
#' # Analyze Functional Enrichment
#' ontologies <- listOntologies("hg38", version = "4.0.4")
#' enrich_bisque4 <- enrichModule(regions, module = "bisque4",
#'                                genome = "hg38",
#'                                file = "bisque4_Module_Enrichment.txt")
#' plotEnrichment(enrich_bisque4,
#'                file = "bisque4_Module_Enrichment_Plot.pdf")
#' }
#'
#' @export
#'
#' @import stringr
#'
#' @importFrom magrittr %>%

getGeneList <- function(regions_annotated, module = NULL,
                        type = c("symbol", "description", "ensemblID", "entrezID"),
                        verbose = TRUE){
        if(!is.null(module)){
                if(verbose){
                        message("[getGeneList] Filtering for regions in ",
                                paste(module, collapse = ", "), " module(s)")
                }
                regions_annotated <- regions_annotated[regions_annotated$module %in% module,]
        }
        type <- match.arg(type)
        if(verbose){
                message("[getGeneList] Getting gene ", type, "s")
        }
        geneList <- str_split(regions_annotated[,paste("gene", type, sep = "_")],
                              pattern = fixed(" | ")) %>%
                purrr::as_vector() %>% unique() %>% sort()
        return(geneList)
}

