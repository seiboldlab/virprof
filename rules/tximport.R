#!/usr/bin/env Rscript

library(cli)

# This creates a snakemake object such as is passed to us from
# snakemake. Useful for debugging interactively.
if (!exists("snakemake")) {
    cli_h1("0. Parsing commandline options (no snakemake object found)")
    library(optparse)
    selfarg <- grep("--file", commandArgs(trailingOnly = FALSE), value = TRUE)
    self <- sub("--file=", "", selfarg)
    sourcedir <- dirname(self)

    parser <- OptionParser(
        option_list = list(
            make_option("--in-samplesheet", help = "PRJ/qiime_mapping.tsv"),
            make_option("--in-multiqc", help = "PRJ.ref_hg38g.qc.quant_salmon_sa.group_ALL.qc_multiqc/multiqc_report_data"),
            make_option("--in-gtf", help = "references/hg38g/ALL.gtf"),
            make_option("--out-counts", default = "out.counts.rds"),
            make_option("--out-transcripts", default ="out.txcounts.rds"),
            make_option("--out-stats", default = "out.stats.rds"),
            make_option("--out-log"),
            make_option("--input-type", default = "Salmon"),
            make_option("--set-version", default = "0.0"),
            make_option("--set-label", default = "unknown"),
            make_option("--threads", default = 4)
        ),
        usage = "usage: %prog [options] count_file ...",
        description = paste(
            sep="\n",
            "Creates unified SummarizedExperiment RDS",
            "",
            "Example:",
            "%prog \\",
            "  --in-sample-sheet PROJ/qiime_mapping.tsv \\",
            "  --in-multiqc PROJ.ref_hg38g.qc.quant_salmon_sa.group_ALL.qc_multiqc/multiqc_report_data \\",
            "  --in-gtf references/hg38g/ALL.gtf",
            "  PROJ.ref_hg38g.qc.quant_salmon_sa/*/quant.sf"
        )
    )
    opt <- parse_args2(parser, args = commandArgs(trailingOnly = TRUE))

    Snakemake <- methods::setClass(
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
            counts = opt$args,
            meta = opt$options$in_samplesheet,
            multiqc = opt$options$in_multiqc,
            gtf = opt$options$in_gtf
        ),
        output = list(
            counts = opt$options$out_counts,
            transcripts = opt$options$out_transcripts,
            stats = opt$options$out_stats
        ),
        log = list(opt$options$out_log),
        threads = opt$options$threads,
        params = list(
            version = opt$options$set_version,
            label = opt$options$set_label,
            input_type = opt$options$input_type
        ),
        scriptdir = sourcedir,
        source = function(...) {
            wd <- getwd()
            setwd(snakemake@scriptdir)
            source(...)
            setwd(wd)
        }
    )

    if (length(snakemake@input$counts) == 0) {
        cli_abort("Missing input count files")
    }
    if (any(!fs::file_exists(snakemake@input$counts))) {
        cli_abort("Some input count files are missing")
    }
    if (!fs::file_exists(snakemake@input$meta)) {
        cli_abort("Sample sheet is missing")
    }
    if (!fs::file_exists(snakemake@input$multiqc)) {
        cli_abort("MultiQC data folder is missing")
    }
    if (!fs::file_exists(snakemake@input$gtf)) {
        cli_abort("GTF file is missing")
    }
}

#' We need to redirect our output to log if running from snakemake...

if (!is.null(snakemake@log[[1]])) {
    logfile <- file(snakemake@log[[1]], open = "wt")
    sink(logfile)
    sink(logfile, type = "message")
}

cli_alert("Importing {snakemake@params$input_type} data into R")

cli_h1("1. Loading packages")
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
library(tibble)
library(DESeq2)
library(lubridate)
library(parallel)
library(future)
library(future.apply)
library(furrr)

# Load faster saveRDS
snakemake@source("_rds.R")

# Set threads
plan(tweak(multicore, workers = snakemake@threads))

get_idcols <- function(table, ids, allow_short = FALSE) {
    ids <- unique(ids)
    if (length(ids) > nrow(table) && !allow_short) {
        rlang::abort("Sampe sheet has fewer entries than unique ids to match")
    }
    res <- names(which(sapply(table, function(x) {
        (
            # column must be unique:
            n_distinct(x) == length(x)
            # and must match ids
            && all(as.character(x) %in% as.character(ids))
        )
    })))
    if (length(res) == 0) {
        print("Sample Sheet:")
        print(as.data.frame(table))
        print("IDs:")
        print(ids)
        rlang::abort(paste(
            "Can't find column in sample sheet matching IDs.",
            "Check whether {project}/qiime_mapping.csv is up to date"
        ))
    }
    res
}



is_header_only_csv <- function(files, col_types) {
    future_sapply(files, function(fn) {
        nrow(read_tsv(
            fn, n_max = 1, guess_max = 1,
            col_types = "cdddd", lazy = FALSE,
            progress = FALSE
        )) == 0
    })
}

# First, run tximport on threads chunks in parallel as parsing
# the files takes quite a while (seconds per file, which adds when
# you have thousands of samples).
parallel_load_txi <- function(files, threads) {
    # from parallel package - splits list into ncl sub-lists with even
    # length; modified to return no empty lists if |x|<ncl
    split_list <- function(x, ncl) {
        ncl <- min(length(x), ncl)
        lapply(splitIndices(length(x), ncl), function(i) x[i])
    }

    txi_parts <- future_lapply(
        split_list(files, threads),
        tximport,
        type = "salmon",
        txOut = TRUE
    )

    # Now merge the three assays
    txi <- list()
    for (assay in c("abundance", "counts", "length")) {
        txi[[assay]] <- do.call(
            cbind,
            lapply(txi_parts, function(x) x[[assay]])
        )
    }
    txi$countsFromAbundance <- txi_parts[[1]]$countsFromAbundance
    return(txi)
}

txi_add_zero_samples <- function(txi, samples) {
    # Add fake data (zero, avg length) for the failed samples
    if (length(samples) == 0) {
        return(txi)
    }
    cols <- length(samples)
    rows <- nrow(txi$abundance)
    zeroes <- matrix(rep(0, cols * rows), ncol = cols)
    colnames(zeroes) <- samples
    lengths <- matrix(rep(rowMeans(txi$length), cols), ncol = cols)
    colnames(lengths) <- samples
    txi$counts <- cbind(txi$counts, zeroes)
    txi$abundance <- cbind(txi$abundance, zeroes)
    txi$length <- cbind(txi$length, lengths)
    return(txi)
}

salmon_load_meta <- function(files, consistency_check = TRUE,
                             logfile_name = NULL) {
    res <- list()
    read_json_nolist <- function(fname, ...) {
        tmp <- jsonlite::read_json(fname, ...)
        tmp[sapply(tmp, function(x) length(x) == 1)]
    }
    res$all <- files %>%
        gsub("/quant.sf", "/aux_info/meta_info.json", .) %>%
        future_map_dfr(
            read_json_nolist,
            simplifyVector = TRUE,
            .id = "idcolumn"
        )
    res$varying <- res$all %>%
        select(-start_time, -end_time) %>%
        select(where(~ length(unique(.x)) != 1))
    res$constant <- res$all %>%
        summarize(
            across(where(~ length(unique(.x)) == 1), ~ unique(.x)[1])
        ) %>%
        as.list()
    if (consistency_check) {
        must_be_constant <- c(
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
        pcols = intersect(colnames(res$varying), must_be_constant)
        if (length(pcols) > 0) {
            if (!is.null(logfile_name)) {
                errorfn <- paste0(logfile_name, ".error.csv")
                message("Writing salmon ERROR information to ", errorfn)
                message(str_glue("Column(s) {pcols} vary but should not"))
                readr::write_csv(res$all, errorfn)
            }
            rlang::abort(
                "Samples were run with multiple references or ",
                "varied parameters. Refusing to aggregate."
            )
        }
    }
    res
}

message("2. ----------- Loading files ----------")
message("2.1. ----------- Loading Sample Sheet ----------")
message("Filename = ", snakemake@input$meta)

sample_sheet <- read_tsv(snakemake@input$meta, show_col_types = FALSE) %>%
    dplyr::select(where(~!all(is.na(.))))

metadata <- list(
    virprof_version = snakemake@params$version,
    pipeline = snakemake@params$label,
    date = now(),
    sample_sheet = sample_sheet
)

message("2.2. ----------- Loading MultiQC Report Data ----------")

multiqc_data_file <- fs::path(snakemake@input$multiqc, "multiqc_data.json")
multiqc <- jsonlite::fromJSON(multiqc_data_file)

extract_counts <- function(col, json) {
    counts <- sapply(
        names(json),
        \(id) as.integer(json[[id]]["Total Sequences"])
    )
    res <- enframe(counts, name="ids", value = col) %>%
        separate("ids", c("ids", "mate"), sep = "\\.", fill="right")
    err <- group_by(res, ids) %>% filter(n_distinct(.data[[col]]) != 1)
    if (nrow(err)) {
        cli_alert("Found errors importing {col}:")
        print(n=1000, err)
        cli_abort("Failed to import read counts from multiqc fastqc report")
    }
    res <- res %>% group_by(ids) %>%
        summarize("{col}":=unique(.data[[col]]), .groups="drop")
}
fastqc_trimmed <- extract_counts(
    "trimmed_reads", multiqc$report_saved_raw_data$multiqc_fastqc
)
fastqc_raw <- extract_counts(
    "raw_reads", multiqc$report_saved_raw_data$multiqc_fastqc_1
)
fastqc <- full_join(fastqc_raw, fastqc_trimmed, by = "ids")
fastqc_idcols <- get_idcols(sample_sheet, fastqc$ids, allow_short = TRUE)
cli_alert_info("FastQC data identified by '{fastqc_idcols}'")

sample_sheet <- sample_sheet %>%
    left_join(fastqc, by = set_names("ids", fastqc_idcols[1]))
err <- filter(sample_sheet, is.na(trimmed_reads) | is.na(raw_reads))
if (nrow(err) > 0) {
    cli_alert("Read counts missing for some samples:")
    print(n=1000, err)
    cli_abort("Incomplete read count information from multiqc fastqc report")
}


message("2.3. ----------- Loading GTF ----------")
message("Filename = ", snakemake@input$gtf)
gr <- rtracklayer::import.gff(snakemake@input$gtf)

if (snakemake@params$input_type == "Salmon") {
    files <- snakemake@input$counts
    names(files) <- gsub(".salmon", "", basename(dirname(files)))
    files <- files[order(names(files))]

    message("3.1. ---------- Checking for failed samples --------")
    no_data <- is_header_only_csv(files, col_types = "cdddd")
    metadata$failed_samples <- if (any(no_data)) {
        message("Warning: excluded ", length(which(no_data)),
                " failed samples:")
        message("  ", paste(names(files[no_data]), collapse = ", "))
        names(files[no_data])
    } else {
        c()
    }

    message("3.2. ----------- Loading meta_info.json files ----------")
    # This will abort if there were mismatches in version or reference
    salmon_meta <- salmon_load_meta(files[!no_data],
                                    logfile = snakemake@log[[1]])
    extra_coldata <- salmon_meta$varying
    metadata$salmon <- salmon_meta$constant

    message("3.3. ----------- Loading quant.sf files ----------")
    txi <- parallel_load_txi(files[!no_data], snakemake@threads)
    txi <- txi_add_zero_samples(txi, metadata$failed_samples)
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

# Extract the sample sheet columns used to identify the samples. This
# could be "unit", or c("unit", "sample") or something else that
# grouping was active on when the counts where generated. We cannot
# currently handle ids that require multiple columns to be combined.
#
# This needs `names(files)` to be the sample names extracted from Salmon
# or RSEM above. Columns must be unique and containing those names.
idcolumns <- get_idcols(sample_sheet, names(files))
message("Samples in output identified by: ",
        paste(idcolumns, collapse = ", "))

coldata <- sample_sheet %>%
    # Make sure all idcolumns are character type (we don't want these
    # numeric)
    mutate(across(all_of(idcolumns), as.character)) %>%
    # Group sample sheet to match active result grouping
    group_by(across(all_of(idcolumns))) %>%
    summarize(
        # Columns with group unique values get that value assigned
        across(where(~ length(unique(.x)) == 1), ~ unique(.x)[1]),
        # Columns with multiple values get them semi colon separated
        across(where(~ length(unique(.x)) > 1),
               ~ paste(as.character(.x), collapse = ";")),
        # Also attach the count of units grouped
        num_units = n(),
        .groups = "drop"
    ) %>%
    arrange(across(all_of(idcolumns)))

if ("idcolumn" %in% colnames(extra_coldata)) {
    # Merge in the extra_coldata gathered above
    coldata <- coldata %>% left_join(
        extra_coldata,
        by = set_names("idcolumn", idcolumns[[1]])
    )
}

message("Coldata:")
print(coldata)

# Sort array columns to match coldata
for (assay in c("counts", "abundance", "length")) {
    message("Sorting ", assay)
    txi[[assay]] <- txi[[assay]][, coldata[[idcolumns[1]]], drop = FALSE]
}

if (snakemake@params$input_type == "ExonSE") {
    stopifnot(all(colnames(se) == coldata[idcolumns[[1]]]))
    colData(se) <- as(coldata, "DataFrame")
    message("5. ----------- Writing RDS with exon se object ----------")
    message("Filename = ", snakemake@output$counts)
    saveRDS(se, snakemake@output$counts)
    saveRDS(se, snakemake@output$transcripts)
} else {
    message("4.2. ----------- Preparing rowData (gene sheet) ----------")
    # only transcript rows
    txmeta <- mcols(gr)[mcols(gr)$type == "transcript", ]
    txmeta <- subset(txmeta, select = -type)
    rownames(txmeta) <- txmeta$transcript_id
    # only rows for which we have counts
    txmeta <- txmeta[rownames(txi$counts), ]
    # remove all-NA columns
    txmeta <- Filter(function(x) !all(is.na(x)), txmeta)

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
        message("6. ----------- Summarizing transcript counts to gene counts",
                " ----------")
        txi_genes <- summarizeToGene(
            txi, txmeta[, c("transcript_id", "gene_id")]
        )
    } else if (snakemake@params$input_type == "RSEM") {
        gene_files <- snakemake@input$counts
        names(gene_files) <- gsub(".genes.results", "", basename(gene_files))
        txi_genes <- tximport(gene_files, type = "rsem",
                              txIn = FALSE, txOut = FALSE)

        if (!is.null(logfile)) {
            ## Something inside of tximport seems to reset the log sink on the
            ## second call. Resetting it here:
            sink(logfile)
            sink(logfile, type = "message")
        }
    }

    message("7. ----------- Assembling SummarizedExperiment ----------")
    # only transcript rows:
    gmeta <-  mcols(gr)[mcols(gr)$type == "gene", ]
    gmeta <- subset(gmeta, select = -type)
    rownames(gmeta) <- gmeta$gene_id
    # only rows for which we have counts:
    gmeta <- gmeta[rownames(txi_genes$counts), ]
    # remove all-NA columns:
    gmeta <- Filter(function(x) !all(is.na(x)), gmeta)

    gse <- SummarizedExperiment(
        assays = txi_genes[c("counts", "abundance", "length")],
        colData = coldata,
        rowData = gmeta,
        metadata = c(
            metadata,
            list(
                # (should be 'no')
                countsFromAbundance = txi_genes$countsFromAbundance
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
            message("Renaming length assay to avgTxLength so DESeq2",
                    " will use for size estimation")
            assayNames(gse)[assayNames(gse) == "length"] <- "avgTxLength"
        }
    }

    message("8. ----------- Collecting mapping metadata ----------")

    mito_genes <-
        rowData(gse) %>%
        as_tibble() %>%
        filter(
            str_detect(gene_name, "^MT-")
        ) %>%
        pull(gene_id)

    ribo_genes <-
        rowData(gse) %>%
        as_tibble() %>%
        filter(
            !gene_id %in% c(mito_genes),
            gene_type %in% c("rRNA", "rRNA_pseudogene")
        ) %>%
        pull(gene_id)

    non_coding_genes <-
        rowData(gse) %>%
        as_tibble() %>%
        filter(
            !gene_id %in% c(mito_genes, ribo_genes),
            gene_type != "protein_coding"
        ) %>%
        pull(gene_id)

    coding_genes <- rowData(gse) %>%
        as_tibble() %>%
        filter(
            !gene_id %in% c(mito_genes, ribo_genes, non_coding_genes)
        ) %>%
        pull(gene_id)

    mapping <- colSums(assay(gse)) %>%
        enframe(name = idcolumns[[1]], value = "count_total") %>%
        left_join(
            colSums(assay(gse)[mito_genes, , drop = FALSE]) %>%
                enframe(name = idcolumns[[1]], value = "count_mito"),
            by = idcolumns[[1]]
        ) %>%
        left_join(
            colSums(assay(gse)[ribo_genes, , drop = FALSE]) %>%
                enframe(name = idcolumns[[1]], value = "count_ribo"),
            by = idcolumns[[1]]
        ) %>%
        left_join(
            colSums(assay(gse)[non_coding_genes, , drop = FALSE]) %>%
                enframe(name = idcolumns[[1]], value = "count_noncoding"),
            by = idcolumns[[1]]
        ) %>%
        left_join(
            colSums(assay(gse)[coding_genes, , drop = FALSE] > 0) %>%
                enframe(name = idcolumns[[1]], value = "num_expr_genes"),
            by = idcolumns[[1]]
        ) %>%
        transmute(
            across(all_of(idcolumns[[1]])),
            mapped_read_count = count_total,
            pct_mito = count_mito / count_total * 100,
            pct_ribo = count_ribo / count_total * 100,
            pct_noncoding = count_noncoding / count_total * 100,
            num_expr_genes
        )

    message("8. ----------- Try Generating PCA data ----------")

    tryCatch({
        dds <- DESeqDataSet(gse[coding_genes, , drop = FALSE], design = ~ 1)
        dds <- dds[, colSums(counts(dds)) > 0, drop = FALSE]
        dds <- dds[rowSums(counts(dds)) > 0, , drop = FALSE]
        # Using poscounts type size factor estimation so we don't fail
        # if there is no gene without zeros.
        dds <- estimateSizeFactors(dds, type = "poscounts")
        vsd <- DESeq2::varianceStabilizingTransformation(dds)
        pca <- plotPCA(vsd, intgroup = idcolumns[[1]], returnData = TRUE) %>%
            select(all_of(idcolumns[[1]]), PC1, PC2)
        mapping <- left_join(mapping, pca, by = idcolumns[[1]])
        message("Success")
    }, error = function(err) {
        print(err)
        message("Failed - Going on")
    })
    metadata(gse)$mapping <- mapping

    message("9. ----------- Writing RDS with gene se object ----------")
    message("Filename = ", snakemake@output$transcripts)
    saveRDS(gse, snakemake@output$counts)

    message("10. ----------- Writing RDS with metadata object ----------")
    message("Filename = ", snakemake@output$stats)
    metadata(gse)$coldata <- colData(gse)
    saveRDS(metadata(gse), snakemake@output$stats)
}
message("done")
