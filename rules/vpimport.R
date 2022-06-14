#!/usr/bin/env Rscript

#' We expect to be called from snakemake script directive, so having
#' `snakemake` object with `snakemake@input` etc containing paths.


#' We also need to redirect our output to log ourselves...
logfile <- file(snakemake@log[[1]], open="wt")
sink(logfile)
sink(logfile, type="message")

snakemake@source("../R/gene_plot.R")
library(lubridate)

if (snakemake@params$task == "create") {
    message("Creating VP object")
    vp <- VirProfFromCSV(
        calls_fname=snakemake@input$calls,
        hits_fname=snakemake@input$hits,
        features_fname=snakemake@input$features
    )
    message("Placing contigs")
    vp <- place_contigs(vp)
    message("Adding contig coverage depths")
    bamlist <- paste(snakemake@input$bam, collapse=",")
    vp <- coverage_depth(vp, bamlist, scaffold=FALSE)
    message("Saving RDS")
    saveRDS(vp, snakemake@output$rds)
} else if (snakemake@params$task == "scaffold_depth") {
    message("Loading VP Object from ", snakemake@input$rds)
    vp <- readRDS(snakemake@input$rds)
    message("Adding scaffold coverage depths")
    bamlist <- paste(snakemake@input$bam, collapse=",")
    vp <- coverage_depth(vp, bamlist, scaffold=TRUE)
    message("Saving RDS to ", snakemake@output$rds)
    saveRDS(vp, snakemake@output$rds)
} else if (snakemake@params$task ==  "combine") {
    vp <- new("VirProf")
    total <- length(snakemake@input$rds)
    for (n in seq_along(snakemake@input$rds)) {
        fname <- snakemake@input$rds[n]
        message(now(), " Loading #", n, "/", total, ": ", fname)
        vp2 <- readRDS(fname)
        message(now(), " merging...")
        vp <- BiocGenerics::combine(vp, vp2)
        if (n %% 50 == 0) {
            message(now(), " Collecting garbage...")
            gc()
            message(now(), " Dedup'ing intermediate")
            vp <- unique(vp)
            message(now(), " Collecting garbage...")
            gc()
        }
    }
    message(now(), " Collecting garbage...")
    gc()
    message(now(), " Dedup'ing final")
    vp <- unique(vp)
    message(now(), " Saving RDS to ", snakemake@output$rds)
    saveRDS(vp, snakemake@output$rds)
}

message(now(), " DONE")

