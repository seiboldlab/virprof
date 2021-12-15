#!/usr/bin/env Rscript
logfile <- file(snakemake@log[[1]], open="wt")
sink(logfile)
sink(logfile, type="message")

library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)

if (FALSE) {

    library(methods)
    Snakemake <- setClass(
        "Snakemake",
        slots = c(
            input = "list",
            output = "list",
            params = "list",
            wildcards = "list",
            threads = "numeric",
            log = "list",
            resources = "list",
            config = "list",
            rule = "character",
            bench_iteration = "numeric",
            scriptdir = "character",
            source = "function"
        )
    )
    snakemake <- Snakemake(
        input = list(
            "gtf" = 'references/hg38/ALL.gtf'
        ),
        output = list(
            "se" = "output.se.rds"
        ),
        params = list(
            "single_gene_only" = FALSE
        ),
        threads = 4
    )
    
}

message("Using ", snakemake@threads, " threads")
options(MulticoreParam=MulticoreParam(workers=snakemake@threads))

message("Loading GTF (", snakemake@input$gtf, ")")
txdb <- makeTxDbFromGFF(snakemake@input$gtf)

message("Extracting Exonic Parts")
message("Single Gene Only set to ",snakemake@params$single_gene_only) 
exons <- exonicParts(txdb, linked.to.single.gene.only = snakemake@params$single_gene_only)

bam <-BamFile(
    file = snakemake@input$bam,
    yieldSize = 100000,
    obeyQname = TRUE  # sorted by query name name (read id)
)

message("Summarizing Overlaps")
overlaps <- summarizeOverlaps(
    exons, bam,
    mode="Union",
    singleEnd=FALSE,
    ignore.strand=FALSE,
    preprocess.reads = invertStrand,
    inter.feature=FALSE,
    fragments=TRUE
)

message("Writing SummarizedExperiment to ", snakemake@output$se)
saveRDS(overlaps, snakemake@output$se)

message("Done")
