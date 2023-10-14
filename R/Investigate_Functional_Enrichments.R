#' Get Ontologies Available in GREAT
#'
#' \code{listOntologies()} obtains a \code{character vector} of the available
#' ontologies for functional enrichment analysis by GREAT for a specified version
#' and genome build.
#'
#' \code{listOntologies()} generates the possible ontologies to use for
#' functional enrichment analysis in [enrichModule()]. Supported ontologies may
#' change over time, so this function queries GREAT using the \pkg{rGREAT}
#' package to get the ones currently available.
#'
#' GREAT supports different genomes depending on the version:
#' \describe{
#'         \item{GREAT v4.0.4}{\code{hg38}, \code{hg19}, \code{mm10}, \code{mm9}}
#'         \item{GREAT v3.0.0}{\code{hg19}, \code{mm10}, \code{mm9}, \code{danRer7}}
#'         \item{GREAT v2.0.2}{\code{hg19}, \code{hg18}, \code{mm9}, \code{danRer7}}
#' }
#'
#' @param genome A \code{character(1)} giving the genome build. Supported genomes
#'         include \code{hg38}, \code{hg19}, \code{hg18}, \code{mm10}, \code{mm9},
#'         and \code{danRer7}. See \code{Details}.
#' @param version A \code{character(1)} specifying the version of GREAT. Possible
#'         values include \code{4.0.4}, \code{3.0.0}, and \code{2.0.2}. Different
#'         versions of GREAT support different genomes. See \code{Details}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{character vector}.
#'
#' @seealso \itemize{
#'         \item [getModules()] to build a comethylation network and identify
#'                 modules of comethylated regions.
#'         \item [annotateModule()] and [getGeneList()] to annotate a set of
#'                 regions with genes and regulatory context and then extract
#'                 the gene symbols or IDs.
#'         \item [enrichModule()] and [plotEnrichment()] to investigate
#'                 functional enrichment of module regions with GREAT.
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
#'
#' @importFrom magrittr %>%

listOntologies <- function(genome = c("hg38", "hg19", "hg18", "mm10", "mm9",
                                      "danRer7"),
                           version = c("4.0.4", "3.0.0", "2.0.2"),
                           verbose =  TRUE){
        genome <- match.arg(genome)
        version <- match.arg(version)
        if((version == "4.0.4" & genome %in% c("hg18", "danRer7")) |
           (version == "3.0.0" & genome %in% c("hg38", "hg18")) |
           (version == "2.0.2" & genome %in% c("hg38", "mm10"))){
                stop("[listOntologies] The ", genome,
                     " genome assembly is not supported for GREAT v", version)
        }
        if(verbose){
                message("[listOntologies] Getting available ontologies for GREAT v",
                        version, " with the ", genome, " genome assembly")
        }
        ontologies <- GRanges("chr1", ranges = IRanges::IRanges(1, end = 2)) %>%
                submitGreatJob(species = genome, request_interval = 0,
                               version = version) %>%
                availableOntologies()
        return(ontologies)
}

#' Analyze Module Functional Enrichment with GREAT
#'
#' \code{enrichModule()} takes a \code{data.frame} of regions with module
#' assignments, filters it for the module(s) of interest, and submits it to
#' GREAT for functional enrichment analysis compared to gene sets in the
#' specified ontologies. The results are then processed and saved as a .txt file.
#'
#' Submission to GREAT is performed by the \pkg{rGREAT} package, which allows
#' for different annotation rules and versions of GREAT. The default
#' \code{basalPlusExt} annotation rule associates a gene with a region if the
#' region is within the basal regulatory domain of the gene (5 kb upstream and 1
#' kb downstream of the TSS) or if it is within 1 Mb upstream or downstream of
#' the TSS and not in the basal regulatory domain of another gene. Other rules
#' include \code{twoClosest} and \code{oneClosest}, which effectively assign the
#' two nearest genes or one nearest genes, respectively. See
#' [rGREAT::submitGreatJob()] for more details.
#'
#' GREAT supports different genomes depending on the version:
#' \describe{
#'         \item{GREAT v4.0.4}{\code{hg38}, \code{hg19}, \code{mm10}, \code{mm9}}
#'         \item{GREAT v3.0.0}{\code{hg19}, \code{mm10}, \code{mm9}, \code{danRer7}}
#'         \item{GREAT v2.0.2}{\code{hg19}, \code{hg18}, \code{mm9}, \code{danRer7}}
#' }
#'
#' Functional enrichment analysis is performed for regions in the module(s) of
#' interest relative to the background set of all regions in \code{regions},
#' including the grey (unassigned) module. These regions are compared for overlap
#' with the regulatory domains of genes annotated to functional gene sets of the
#' ontologies of interest. The default ontologies include
#' \code{GO Molecular Function}, \code{GO Biological Process},
#' \code{GO Cellular Component}, \code{Mouse Phenotype}, and
#' \code{Human Phenotype}. Initial enrichment results are filtered for terms
#' with a minimum number of overlaps with the background set of regions,
#' p-values are adjusted for multiple comparisons using the specified method,
#' and then the results are filtered again for a minimum number of overlaps with
#' the module(s) of interest and a p-value threshold. Finally, gene symbols are
#' obtained for the significant gene sets, and the result is saved as a .txt
#' file.
#'
#' @param regions A \code{data.frame} of regions with module assignments,
#'         typically obtained from [getModules()].
#' @param module A \code{character} giving the name of one or more modules to
#'         analyze functional enrichment. If null, all modules will be analyzed,
#'         except the grey (unassigned) module.
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
#' @param ontologies A \code{character} giving the ontology databases to query.
#'         Default ontologies include \code{GO Molecular Function},
#'         \code{GO Biological Process}, \code{GO Cellular Component},
#'         \code{Mouse Phenotype}, and \code{Human Phenotype}. All possible
#'         ontologies can be obtained with [listOntologies()].
#' @param min_background_region_hits A \code{numeric(1)} giving the minimum
#'         number of overlaps of gene set regulatory domains with background
#'         regions to include that gene set in the results. This affects the
#'         results used when adjusting p-values for multiple comparisons.
#' @param adjMethod A \code{character(1)} specifying the method for adjusting
#'         p-values, Potential values include \code{fdr}, \code{holm},
#'         \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH},
#'         \code{BY}, and \code{none}.
#' @param min_region_hits A \code{numeric(1)} giving the minimum number of
#'         overlaps of gene set regulatory domains with module regions to
#'         include that gene set in the results.
#' @param pvalue_threshold A \code{numeric(1)} giving the maximum p-value
#'         for enrichment to include a gene set in the results.
#' @param save A \code{logical(1)} indicating whether to save the
#'         \code{data.frame}.
#' @param file A \code{character(1)} giving the file name (.txt) for the saved
#'         \code{data.frame}.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{data.frame} of functional enrichment results.
#'
#' @seealso \itemize{
#'         \item [getModules()] to build a comethylation network and identify
#'                 modules of comethylated regions.
#'         \item [annotateModule()] and [getGeneList()] to annotate a set of
#'                 regions with genes and regulatory context and then extract
#'                 the gene symbols or IDs.
#'         \item [listOntologies()] and [plotEnrichment()] to investigate
#'                 functional enrichment of module regions with GREAT.
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
#' @import utils
#'
#' @importFrom magrittr %>%

enrichModule <- function(regions, module = NULL,
                         genome = c("hg38", "hg19", "hg18", "mm10", "mm9",
                                    "danRer7"),
                         includeCuratedRegDoms = FALSE,
                         rule = c("basalPlusExt", "twoClosest", "oneClosest"),
                         adv_upstream = 5, adv_downstream = 1, adv_span = 1000,
                         adv_twoDistance = 1000, adv_oneDistance = 1000,
                         version = c("4.0.4", "3.0.0", "2.0.2"),
                         ontologies = c("GO Molecular Function",
                                        "GO Biological Process",
                                        "GO Cellular Component",
                                        "Mouse Phenotype", "Human Phenotype"),
                         min_background_region_hits = 5,
                         adjMethod = c("fdr", "holm", "hochberg", "hommel",
                                       "bonferroni", "BH", "BY", "none"),
                         min_region_hits = 2, pvalue_threshold = 0.01,
                         save = TRUE, file = "Module_Enrichment_Results.txt",
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
        if((version == "4.0.4" & genome %in% c("hg18", "danRer7")) |
           (version == "3.0.0" & genome %in% c("hg38", "hg18")) |
           (version == "2.0.2" & genome %in% c("hg38", "mm10"))){
                stop("[enrichModule] The ", genome,
                     " genome assembly is not supported for GREAT v", version)
        }
        adjMethod <- match.arg(adjMethod)
        if(verbose){
                message("[enrichModule] Using GREAT v", version, " with the ",
                        genome, " genome assembly and the ", rule, " rule")
        }
        if(genome == "danRer7"){
                ontologies <- ontologies[!ontologies %in% c("Mouse Phenotype",
                                                            "Human Phenotype")]
        }
        gr <- with(moduleRegions,
                   GRanges(chr, ranges = IRanges::IRanges(start, end = end)))
        bg <- with(regions,
                   GRanges(chr, ranges = IRanges::IRanges(start, end = end)))
        job <- submitGreatJob(gr, bg = bg, species = genome,
                              includeCuratedRegDoms = includeCuratedRegDoms,
                              rule = rule, adv_upstream = adv_upstream,
                              adv_downstream = adv_downstream,
                              adv_span = adv_span,
                              adv_twoDistance = adv_twoDistance,
                              adv_oneDistance = adv_oneDistance,
                              request_interval = 0, version = version)
        if(verbose){
                message("[enrichModule] Getting results and adjusting p-values using the ",
                        adjMethod, " method")
        }
        enrichTables <- getEnrichmentTables(job, ontology = ontologies)
        results <- rlist::list.rbind(enrichTables)
        rownames(results) <- 1:nrow(results)
        pattern <- c("hyper_" = "", "id" = "ID", "name" = "term",
                     "total_regions" = "background_region_hits",
                     "expected" = "expected_region_hits",
                     "foreground_region_hits" = "region_hits",
                     "region_set_coverage" = "region_coverage",
                     "foreground_gene_hits" = "gene_hits",
                     "total_genes_annotated" = "term_genes", "raw_pvalue" = "p")
        colnames(results) <- colnames(results) %>% str_to_lower() %>%
                str_replace_all(pattern = pattern)
        results$ontology <- rep(names(enrichTables), sapply(enrichTables, nrow))
        results <- subset(results, background_region_hits >= min_background_region_hits)
        results$log_p <- -log10(results$p)
        results$adj_p <- stats::p.adjust(results$p, method = adjMethod)
        results$log_adj_p <- -log10(results$adj_p)
        results <- subset(results, region_hits >= min_region_hits &
                                  p < pvalue_threshold)
        if(verbose){
                message("[enrichModule] Extracting genes for enriched terms")
        }
        results$genes <- mapply(function(x,y){
                R.devices::suppressGraphics(
                        plotRegionGeneAssociationGraphs(job, type = 1,
                                                        ontology = x, termID = y,
                                                        request_interval = 0,
                                                        max_tries = 10000,
                                                        verbose = FALSE)) %>%
                        .$gene %>% unique() %>% sort() %>% paste(collapse = ", ")
        }, x = results$ontology, y = results$ID, USE.NAMES = FALSE)
        results <- results[order(results$p), c("ID", "term", "ontology",
                                               "background_region_hits",
                                               "expected_region_hits",
                                               "region_hits", "fold_enrichment",
                                               "region_coverage",
                                               "term_region_coverage",
                                               "gene_hits", "genes",
                                               "background_gene_hits",
                                               "term_genes", "p", "log_p",
                                               "adj_p", "log_adj_p")]
        if(save){
                if(verbose){
                        message("[enrichModule] Saving file as ", file)
                }
                write.table(results, file = file, quote = FALSE, sep = "\t",
                            row.names = FALSE)
        }
        return(results)
}

#' Plot Functional Enrichment Results
#'
#' \code{plotEnrichment()} takes a \code{data.frame} of enrichment results from
#' [enrichModule()], plots the log p-values in a bar plot, and saves it as a .pdf.
#'
#' \code{plotEnrichment()} is designed to be used in combination with
#' [enrichModule()]. The top 15 gene sets are plotted by default, but this can
#' be expanded if needed. A \code{ggplot} object is produced and can be edited
#' outside of this function if desired.
#'
#' @param enrichment A \code{data.frame} of functional enrichment results,
#'         typically obtained from [enrichModule()].
#' @param nTerms A \code{numeric(1)} specifying the number of terms to include in
#'         the plot.
#' @param fill A \code{character(1)} giving the color of the bars.
#' @param xlim A \code{numeric(2)} specifying the limits of the x-axis.
#' @param nBreaks A \code{numeric(1)} indicating the number of breaks to use in
#'         the x-axis.
#' @param axis.title.x.size A \code{numeric(1)} with the size of the x-axis title.
#' @param axis.text.x.size A \code{numeric(1)} giving the size of the x-axis text.
#' @param axis.text.y.size A \code{numeric(1)} giving the size of the y-axis text.
#' @param save A \code{logical(1)} indicating whether to save the plot.
#' @param file A \code{character(1)} giving the file name (.pdf) for the saved
#'         plot.
#' @param width A \code{numeric(1)} specifying the width in inches of the saved
#'         plot.
#' @param height A \code{numeric(1)} specifying the height in inches of the
#'         saved plot.
#' @param verbose A \code{logical(1)} indicating whether messages should be
#'         printed.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \itemize{
#'         \item [getModules()] to build a comethylation network and identify
#'                 modules of comethylated regions.
#'         \item [annotateModule()] and [getGeneList()] to annotate a set of
#'                 regions with genes and regulatory context and then extract
#'                 the gene symbols or IDs.
#'         \item [listOntologies()] and [enrichModule()], to investigate
#'                 functional enrichment of module regions with GREAT.
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
#' @import ggplot2
#'
#' @importFrom scales breaks_pretty

plotEnrichment <- function(enrichment, nTerms = 15, fill = "#132B43",
                           xlim = NULL, nBreaks = 4, axis.title.x.size = 20,
                           axis.text.x.size = 16, axis.text.y.size = 16,
                           save = TRUE, file = "Module_Enrichment_Plot.pdf",
                           width = 8, height = 6, verbose = TRUE){
        if(verbose){
                message("[plotEnrichment] Plotting module enrichment from GREAT")
        }
        enrichment <- enrichment[order(enrichment$log_p, decreasing = TRUE),]
        if(nrow(enrichment) > nTerms){
                enrichment <- enrichment[1:nTerms,]
        }
        enrichment$term <- factor(enrichment$term,
                                  levels = rev(unique(enrichment$term)))
        barplot <- ggplot() +
                geom_col(aes(x = term, y = log_p), data = enrichment,
                         fill = fill) +
                coord_flip(ylim = xlim) +
                scale_y_continuous(breaks = breaks_pretty(n = nBreaks),
                                   expand = expansion(c(0.004, 0.03))) +
                theme_bw(base_size = 25) +
                theme(legend.position = "none", panel.grid.major = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25),
                      axis.ticks = element_line(size = 1.25),
                      panel.grid.minor = element_blank(),
                      strip.background = element_blank(),
                      axis.text.x = element_text(color = "black",
                                                 size = axis.text.x.size),
                      axis.text.y = element_text(color = "black",
                                                 size = axis.text.y.size),
                      axis.title.x = element_text(size = axis.title.x.size),
                      axis.title.y = element_blank(),
                      plot.margin = unit(c(1,1,0.5,1), "lines")) +
                ylab(expression(-log[10]*(italic(p)-value)))
        if(save){
                if(verbose){
                        message("[plotEnrichment] Saving plot as ", file)
                }
                ggsave(filename = file, plot = barplot, dpi = 600, width = width,
                       height = height, units = "in")
        }
        return(barplot)
}

