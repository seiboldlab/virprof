#!/usr/bin/env Rscript

#' We expect to be called from snakemake script directive, so having
#' `snakemake` object with `snakemake@input` etc containing paths.

#' We also need to redirect our output to log ourselves...
logfile <- file(snakemake@log[[1]], open="wt")
sink(logfile)
sink(logfile, type="message")

#' Test data for interactive debugging
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
        input = list('test.ref_hg38gtest.qc.quant_salmon_sa/A1.salmon/quant.sf',
                     'test.ref_hg38gtest.qc.quant_salmon_sa/A2.salmon/quant.sf',
                     'test.ref_hg38gtest.qc.quant_salmon_sa/B.salmon/quant.sf',
                     'test.ref_hg38gtest.qc.quant_salmon_sa/C.salmon/quant.sf',
                     'test.ref_hg38gtest.qc.quant_salmon_sa/D.salmon/quant.sf',
                     'references/hg38gtest/ALL.gtf',
                     'test/qiime_mapping.tsv',
                     "counts" = c('test.ref_hg38gtest.qc.quant_salmon_sa/A1.salmon/quant.sf',
                                  'test.ref_hg38gtest.qc.quant_salmon_sa/A2.salmon/quant.sf',
                                  'test.ref_hg38gtest.qc.quant_salmon_sa/B.salmon/quant.sf',
                                  'test.ref_hg38gtest.qc.quant_salmon_sa/C.salmon/quant.sf',
                                  'test.ref_hg38gtest.qc.quant_salmon_sa/D.salmon/quant.sf'),
                     "gtf" = 'references/hg38gtest/ALL.gtf',
                     "meta" = 'test/qiime_mapping.tsv'),
        output = list('test.ref_hg38gtest.qc.quant_salmon_sa.group_ALL.tximport_vp/ALL.gene_counts.rds',
                      'test.ref_hg38gtest.qc.quant_salmon_sa.group_ALL.tximport_vp/ALL.tx_counts.rds',
                      "counts" = 'test.ref_hg38gtest.qc.quant_salmon_sa.group_ALL.tximport_vp/ALL.gene_counts.rds',
                      "transcripts" = 'test.ref_hg38gtest.qc.quant_salmon_sa.group_ALL.tximport_vp/ALL.tx_counts.rds'),
        params = list('input_type', 'version',
                      "input_type" = 'Salmon',
                      "version" = '0.1.0'),
        wildcards = list('test.ref_hg38gtest.qc.quant_salmon_sa.group_ALL.', 'ALL',
                         "_YMP_DIR" = 'test.ref_hg38gtest.qc.quant_salmon_sa.group_ALL.',
                         "target" = 'ALL'),
        threads = 2,
        log = list('test.ref_hg38gtest.qc.quant_salmon_sa.group_ALL.tximport_vp/ALL.log'),
        resources = list('tmpdir', 'mem', 'walltime', 'mem_mb', 'mem_gb',
                         "tmpdir" = '/tmp',
                         "mem" = 8589934592,
                         "walltime" = '0-23:59:59',
                         "mem_mb" = 8192,
                         "mem_gb" = 8),
        config = list(),
        rule = 'vp_tximport_salmon',
        bench_iteration = as.numeric(NA),
        scriptdir = '/Seibold/home/pruessee/projects/virprof/projects/virprof/rules',
        source = function(...){
            wd <- getwd()
            setwd(snakemake@scriptdir)
            source(...)
            setwd(wd)
        }
    )
}


read_json_nolist <- function(fname, ...) {
    tmp<-jsonlite::read_json(fname, ...)
    tmp[sapply(tmp, function(x) length(x) == 1)]
}


message("Importing ", snakemake@params$input_type, " data into R")

message("1. ----------- Loading packages ----------")
library(tximport)
library(readr)
library(GenomicFeatures)
library(rtracklayer)
library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(purrr)
library(jsonlite)

metadata <- list(
    date = date(),
    virprof_version = snakemake@params$version,
    pipeline = snakemake@params$label
)

message("2. ----------- Loading files ----------")
message("2.1. ----------- Loading Sample Sheet ----------")
message("Filename = ", snakemake@input$meta)
samples <- read.csv(snakemake@input$meta, sep="\t")
samples <- Filter(function(x)!all(is.na(x)), samples)  # remove all-NA columns

message("2.2. ----------- Loading GTF ----------")
message("Filename = ", snakemake@input$gtf)
gr <- rtracklayer::import.gff(snakemake@input$gtf)

if (snakemake@params$input_type == "Salmon") {
    message("3.1. ----------- Loading quant.sf files ----------")
    files <- snakemake@input$counts
    names(files) <- gsub(".salmon", "", basename(dirname(snakemake@input$counts)))
    files <- files[order(names(files))]
    txi <- tximport(files, type="salmon", txOut=TRUE)

    message("3.2. ----------- Loading meta_info.json files ----------")
    salmon_all_meta <-
        files %>%
        gsub("/quant.sf", "/aux_info/meta_info.json", .) %>%
        purrr::map_df(read_json_nolist, simplifyVector = TRUE, .id = "idcolumn")

    # Extract data varying per sample
    extra_coldata <- salmon_all_meta %>%
        select(where(~ length(unique(.x)) != 1)) %>%
        select(-start_time, -end_time)

    must_be_identical <- c(
        "index_decoy_seq_hash",
        "index_decoy_name_hash",
        "num_decoy_targets",
        "index_seq_hash",
        "index_name_hash",
        "num_valid_targets",
        "seq_bias_correct",
        "gc_bias_correct",
        "salmon_version"
    )

    if (length(intersect(colnames(extra_coldata), must_be_identical)) > 0) {
        errorfn = paste0(logfile, ".error.csv")
        message("Samples were run with multiple references or varried parameters.",
                "Refusing to aggregate.",
                "Writing coldata to ", errorfn)
        readr::write_csv(extra_coldata, errorfn)
        stop()
    }

    # Extract data constant across dataset
    metadata$salmon <- salmon_all_meta %>%
        summarize(
            across(where(~ length(unique(.x)) == 1), ~ unique(.x)[1])
        ) %>%
        as.list()
} else if (snakemake@params$input_type == "RSEM") {
    files <- snakemake@input$transcripts
    names(files) <- gsub(".isoforms.results", "", basename(files))
    txi <- tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)
    extra_coldata <- data.frame(idcolumn = character())
} else if (snakemake@params$input_type == "ExonSE") {
    files <- snakemake@input$counts
    names(files) <- gsub(".exon.se.rds", "", basename(files))
    sel <- vector("list", length(files))
    for (i in seq_along(files)) {
        message("Reading ", files[[i]], " ...")
        sel[[i]] <- readRDS(files[[i]])
    }
    extra_coldata <- data.frame(idcolumn = character())
    se <- do.call("cbind", sel)
    colnames(se) <- names(files)
}

message("4. ----------- Assembling SummarizedExperiment ----------")

message("4.1. ----------- Preparing colData (sample sheet) -----------")
idcolumn <- names(which(sapply(samples, function(x) all(sort(x)==names(files)))))
if (length(idcolumn) == 0) {
    stop("The sample sheet columns and file names didn't match up. Something is wrong. Bailing out.")
}

coldata <- samples %>%
    group_by(across(all_of(idcolumn))) %>%
    summarize(
        across(where(~ length(unique(.x)) == 1), ~ unique(.x)[1]),
        across(where(~ length(unique(.x)) > 1), ~ paste(as.character(.x), collapse=";")),
        num_units=n(),
        .groups="drop"
    ) %>%
    arrange(across(all_of(idcolumn))) %>%
    left_join(
        extra_coldata,
        by = set_names("idcolumn", idcolumn[[1]])
    )

if (snakemake@params$input_type == "ExonSE") {
    stopifnot(all(colnames(se) == coldata[idcolumn[[1]]]))
    colData(se) <- as(coldata, "DataFrame")
    message("5. ----------- Writing RDS with exon se object ----------")
    message("Filename = ", snakemake@output$counts)
    saveRDS(se, snakemake@output$counts)
    saveRDS(se, snakemake@output$transcripts)
    } else {

message("4.2. ----------- Preparing rowData (gene sheet) ----------")
txmeta <- mcols(gr)[mcols(gr)$type=="transcript", ]  # only transcript rows
txmeta <- subset(txmeta, select = -type)
rownames(txmeta) <- txmeta$transcript_id  # set names
txmeta <- txmeta[rownames(txi$counts), ]  # only rows for which we have counts
txmeta <- Filter(function(x)!all(is.na(x)), txmeta)  # remove all-NA columns

message("4.3. ----------- Creating object ----------")
se <- SummarizedExperiment(
    assays = txi[c("counts", "abundance", "length")],
    rowData = txmeta,
    colData = coldata,
    metadata = c(
        metadata,
        list(
            countsFromAbundance = txi$countsFromAbundance  # should be no
        )
    )
)

message("5. ----------- Writing RDS with transcript se object ----------")
message("Filename = ", snakemake@output$transcripts)
saveRDS(se, snakemake@output$transcripts)

if (snakemake@params$input_type == "Salmon") {
    message("6. ----------- Summarizing transcript counts to gene counts ----------")
    txi_genes <- summarizeToGene(txi, txmeta[,c("transcript_id", "gene_id")])
} else if (snakemake@params$input_type == "RSEM") {
    gene_files <- snakemake@input$counts
    names(gene_files) <- gsub(".genes.results", "", basename(gene_files))
    txi_genes <- tximport(gene_files, type = "rsem", txIn = FALSE, txOut = FALSE)
   
    ## Something inside of tximport seems to reset the log sink on the
    ## second call. Resetting it here:
    sink(logfile)
    sink(logfile, type="message")
}

message("7. ----------- Assembling SummarizedExperiment ----------")
gmeta <-  mcols(gr)[mcols(gr)$type=="gene", ]  # only transcript rows
gmeta <- subset(gmeta, select = -type)
rownames(gmeta) <- gmeta$gene_id  # set names
gmeta <- gmeta[rownames(txi_genes$counts), ]  # only rows for which we have counts
gmeta <- Filter(function(x)!all(is.na(x)), gmeta)  # remove all-NA columns

gse <- SummarizedExperiment(
    assays = txi_genes[c("counts", "abundance", "length")],
    colData = coldata,
    rowData = gmeta,
    metadata = c(
        metadata,
        list(
            countsFromAbundance = txi_genes$countsFromAbundance  # should be no
        )
    )
)

message("Rounding counts to keep DESeq2 happy")
assay(gse) <- round(assay(gse))
mode(assay(gse)) <- "integer"

## Rename length assay IFF we are having counts, not TPM
## (not sure if otherwise is possible with Salmon, but since this is
## checked inside of deseq/tximeta, let's do check here as well).
if (snakemake@params$input_type == "Salmon") {
    if (txi_genes$countsFromAbundance == "no") {
        message("Renaming length assay to avgTxLength so DESeq2 will use for size estimation")
        assayNames(gse)[assayNames(gse) == "length"] <- "avgTxLength"
    }
}

message("8. ----------- Writing RDS with gene se object ----------")
message("Filename = ", snakemake@output$transcripts)
saveRDS(gse, snakemake@output$counts)

}
message("done")
