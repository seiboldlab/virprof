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
library(fs)
library(stringr)
library(tidyr)

metadata <- list(
    date = date(),
    virprof_version = snakemake@params$version,
    pipeline = snakemake@params$label
)

message("2. ----------- Loading files ----------")
message("2.1. ----------- Loading Sample Sheet ----------")
message("Filename = ", snakemake@input$meta)
samples <- read_tsv(snakemake@input$meta)
samples <- Filter(function(x)!all(is.na(x)), samples)  # remove all-NA columns

message("2.2. ----------- Loading GTF ----------")
message("Filename = ", snakemake@input$gtf)
gr <- rtracklayer::import.gff(snakemake@input$gtf)

message("2.3. ----------- Loading MultiQC Report Data ----------")

# FastQC
for (n in c(1,2)) {
    suffix = c("", "_1")[n]
    round = c("raw", "trimmed")[n]
    fastqc_fn <- path(
        snakemake@input$multiqc,
        str_glue("multiqc_fastqc{suffix}.txt")
    )
    if (fs::file_exists(fastqc_fn)) {
        df <- read_tsv(fastqc_fn, show_col_types = FALSE) %>%
            # Sequence length can be `95` or `95-151`, split here
            tidyr::separate(
                col = `Sequence length`,
                into = c("read_len_min", "read_len_max"),
                convert = TRUE,
                fill = "right"
            ) %>%
            {print(.); .} %>%
            transmute(
                sample = sub(".R[12]$", "", Sample),
                mate = if_else(str_ends(Sample, "R2"), "R2", "R1"),
                num_reads = `Total Sequences`,
                read_len_min,
                read_len_max = if_else(
                    is.na(read_len_max), read_len_min, read_len_max
                ),
                read_len_avg = avg_sequence_length,
                pct_gc = `%GC`,
                # This is the percentage of unique reads (dedup/total)
                pct_unique = total_deduplicated_percentage,
                trimmed = round == "trimmed"
            )
        if (is.null(metadata$fastqc)) {
            metadata$fastqc <- df
        } else {
            metadata$fastqc <- bind_rows(metadata$fastqc, df)
        }
    }
}

if (snakemake@params$input_type == "Salmon") {
    message("3.1. ---------- Checking for empty samples --------")
    files <- snakemake@input$counts
    names(files) <- gsub(".salmon", "", basename(dirname(files)))
    files <- files[order(names(files))]

    no_data <- sapply(files, function(fn) {
        nrow(read_tsv(fn, n_max = 1, guess_max = 1, col_types = "cdddd")) == 0
    })
    if (any(no_data)) {
        message("Warning: excluded ", length(which(no_data)), " empty samples:")
        message("  ", paste(names(files[no_data]), collapse = ", "))
        metadata$excluded_empty_samples <- names(files[no_data])
        files <- files[!no_data]
    }

    message("3.2. ----------- Loading quant.sf files ----------")
    txi <- tximport(files, type="salmon", txOut=TRUE)

    message("3.3. ----------- Loading meta_info.json files ----------")
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

# Extract the columns used with the active grouping to identify the samples:
idcolumns <- names(which(sapply(samples, function(x) {
    n_distinct(x) == length(x) && # column must be unique
        all(names(files) %in% as.character(x)) # must identify each sample
})))

if (length(idcolumns) == 0) {
    message("The sample sheet columns and file names didn't match up. Something is wrong. Bailing out.")
    message("samples:")
    print(as.data.frame(samples))
    message("files:")
    print(as.data.frame(files))
    stop("Unable to determine id columns")
}

coldata <- samples %>%
    # Make sure all idcolumns are character type (we don't want these numeric)
    mutate(across(all_of(idcolumns), as.character)) %>%
    # Remove samples (really, results) filtered out above
    filter(.data[[idcolumns[1]]] %in% names(files)) %>%
    # Group sample sheet to match active result grouping
    group_by(across(all_of(idcolumns))) %>%
    summarize(
        # Columns with group unique values get that value assigned
        across(where(~ length(unique(.x)) == 1), ~ unique(.x)[1]),
        # Columns with multiple values get them semi colon separated
        across(where(~ length(unique(.x)) > 1), ~ paste(as.character(.x), collapse=";")),
        # Also attach the count of units grouped
        num_units=n(),
        .groups="drop"
    ) %>%
    arrange(across(all_of(idcolumns))) %>%
    # Merge in the extra_coldata gathered above
    left_join(
        extra_coldata,
        by = set_names("idcolumn", idcolumns[[1]])
    )

if (snakemake@params$input_type == "ExonSE") {
    stopifnot(all(colnames(se) == coldata[idcolumns[[1]]]))
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

    message("8. ----------- Writing RDS with metadata object ----------")
    message("Filename = ", snakemake@output$transcripts)
    saveRDS(metadata(gse), snakemake@output$stats)
}
message("done")
