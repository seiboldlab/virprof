#!/usr/bin/env Rscript

#' We expect to be called from snakemake script directive, so having
#' `snakemake` object with `snakemake@input` etc containing paths.

#' We also need to redirect our output to log ourselves...
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)
sink(logfile, type = "message")

read_json_nolist <- function(fname, ...) {
    tmp <- jsonlite::read_json(fname, ...)
    tmp[sapply(tmp, function(x) length(x) == 1)]
}

#' from parallel package - splits list into ncl sub-lists with even
#' length; modified to return no empty lists if |x|<ncl
#'
split_list <- function(x, ncl) {
    ncl <- min(length(x), ncl)
    lapply(splitIndices(length(x), ncl), function(i) x[i])
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
library(tibble)
library(DESeq2)
library(lubridate)
library(parallel)
library(future)
library(future.apply)
library(furrr)

if (FALSE) {
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
    setwd("/Seibold/tmp/pipeline/work")
    project <- "gala"
    snakemake <- Snakemake(
        scriptdir = "/Seibold/tmp/pipeline/work/virprof/rules",
        source = function(...) {
            wd <- getwd()
            setwd(snakemake@scriptdir)
            source(...)
            setwd(wd)
        }
    )
    snakemake@input$counts <- fs::dir_ls(
        path = paste0(project, ".ref_hg38g.qc.quant_salmon_sa"),
        glob = "*/quant.sf",
        recurse = TRUE
    )
    snakemake@input$meta <- file.path(project, "qiime_mapping.tsv")
    snakemake@input$multiqc <- paste0(
        project,
        ".ref_hg38g.qc.quant_salmon_sa.group_ALL.qc_multiqc/multiqc_report_data"
    )
    snakemake@input$gtf <- "references/hg38g/ALL.gtf"
    snakemake@threads <- 16
    snakemake@params$version <- "0.0.0"
    snakemake@params$label <- "testing manually"
    snakemake@params$input_type <- "salmon"
}

# Load faster saveRDS
snakemake@source("_rds.R")

# Set threads
plan(tweak(multicore, workers = snakemake@threads))

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

# FastQC
metadata$fastqc <- NULL
for (n in c(1,2)) {
    suffix <- c("", "_1")[n]
    round <- c("raw", "trimmed")[n]
    fastqc_fn <- fs::path(
        snakemake@input$multiqc,
        str_glue("multiqc_fastqc{suffix}.txt")
    )
    if (fs::file_exists(fastqc_fn)) {
        message("  Loading ", round, " read fastqc data")
        df <- read_tsv(fastqc_fn, show_col_types = FALSE) %>%
            # Sequence length can be `95` or `95-151`, split here
            tidyr::separate(
                col = `Sequence length`,
                into = c("read_len_min", "read_len_max"),
                convert = TRUE,
                fill = "right"
            ) %>%
            transmute(
                # filename is just sample.R[12].fq.gz, ignoring
                fastqc_id = sub(".R[12]$", "", Sample),
                mate = if_else(str_ends(Sample, "R2"), "R2", "R1"),
                num_reads = `Total Sequences`,
                read_len_min,
                read_len_max = ifelse(
                    is.na(read_len_max), read_len_min, read_len_max
                ),
                read_len_avg = avg_sequence_length,
                pct_gc = `%GC`,
                # This is the percentage of unique reads (dedup/total)
                pct_unique = total_deduplicated_percentage,
                trimmed = round == "trimmed"
            )
        # Find columns identifying the fastqc files in sample sheet
        fastqc_idcolumns <- names(which(sapply(sample_sheet, function(x) {
            n_distinct(x) == length(x) && # column must be unique
                # must identify each sample
                all(df$fastqc_id %in% as.character(x))
        })))
        if (length(fastqc_idcolumns) == 0) {
            rlang::abort(paste(
                "Can't find columns identifying files.",
                "Check whether {project}/qiime_mapping.csv is up to date"
            ))
        }
        message("    FastQC data in ", round, " identified by: ",
                paste(fastqc_idcolumns, collapse = ", "))
        # Extract paths
        #
        # FIXME: This is a bad hack. We should get the detected fwd
        #        and rev read from YMP somehow. Also, this won't work
        #        for SRR sources.
        fastq_id_to_path <- sample_sheet %>%
            pivot_longer(
                cols = where(~any(str_detect(., "(fastq|fq).gz"),
                                  # must return T/F, no NA to where:
                                  na.rm = TRUE )),
                values_to = "fastq_file_path",
                names_to = "path_col"
            ) %>%
            mutate(
                # Deduce mate from file name
                mate = if_else(
                    str_detect(basename(fastq_file_path), "(_|\\.)R1"),
                    "R1", "R2"
                )
            ) %>%
            # Remove anything we can't uniquely match
            group_by(across(fastqc_idcolumns[1]), mate) %>%
            filter(n() == 1) %>%
            ungroup()
        # Merge this into the fastqc df
        df <- df %>%
            left_join(
                fastq_id_to_path,
                by = c(fastqc_id = fastqc_idcolumns[1], "mate")
            )
        if (any(is.na(df$fastq_file_path))) {
            message(
                "WARNING: ",
                "Failed to identify fastq file paths for all samples"
            )
            message("---- BEGIN fastqc data w/o file path ----")
            print(filter(df, is.na(fastq_file_path)))
            message("---- END fastqc data w/o file path ----")
        }
        # Rename sample column to it's name from the sample sheet
        df <- df %>% dplyr::rename("{fastqc_idcolumns[1]}" := fastqc_id)
        if (is.null(metadata$fastqc)) {
            metadata$fastqc <- df
        } else {
            metadata$fastqc <- bind_rows(metadata$fastqc, df)
        }
    }
}

message("2.3. ----------- Loading GTF ----------")
message("Filename = ", snakemake@input$gtf)
gr <- rtracklayer::import.gff(snakemake@input$gtf)

if (snakemake@params$input_type == "Salmon") {
    message("3.1. ---------- Checking for failed samples --------")
    files <- snakemake@input$counts
    names(files) <- gsub(".salmon", "", basename(dirname(files)))
    files <- files[order(names(files))]

    no_data <- future_sapply(files, function(fn) {
        nrow(read_tsv(
            fn, n_max = 1, guess_max = 1,
            col_types = "cdddd", lazy = FALSE,
            progress = FALSE
        )) == 0
    })

    if (any(no_data)) {
        message("Warning: excluded ", length(which(no_data)),
                " failed samples:")
        message("  ", paste(names(files[no_data]), collapse = ", "))
        metadata$failed_samples <- names(files[no_data])
    }

    message("3.2. ----------- Loading quant.sf files ----------")

    # First, run tximport on threads chunks in parallel as parsing
    # the files takes quite a while (seconds per file, which adds when
    # you have thousands of samples).
    txi_parts <- future_lapply(
        split_list(files[!no_data], snakemake@threads),
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
    rm(txi_parts)

    # Add fake data (zero, avg length) for the failed samples
    if (length(metadata$failed_samples) > 0) {
        cols <- length(metadata$failed_samples)
        rows <- nrow(txi$abundance)
        zeroes <- matrix(rep(0, cols * rows), ncol = cols)
        colnames(zeroes) <- metadata$failed_samples
        lengths <- matrix(rep(rowMeans(txi$length), cols), ncol = cols)
        colnames(lengths) <- metadata$failed_samples
        txi$counts <- cbind(txi$counts, zeroes)
        txi$abundance <- cbind(txi$abundance, zeroes)
        txi$length <- cbind(txi$length, lengths)
    }

    message("3.3. ----------- Loading meta_info.json files ----------")
    salmon_all_meta <-
        files[!no_data] %>%
        gsub("/quant.sf", "/aux_info/meta_info.json", .) %>%
        future_map_dfr(read_json_nolist, simplifyVector = TRUE,
                      .id = "idcolumn")

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
        errorfn <- paste0(logfile, ".error.csv")
        message("Samples were run with multiple references or ",
                "varied parameters. Refusing to aggregate.",
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

# Extract the sample sheet columns used to identify the samples. This
# could be "unit", or c("unit", "sample") or something else that
# grouping was active on when the counts where generated. We cannot
# currently handle ids that require multiple columns to be combined.
#
# This needs `names(files)` to be the sample names extracted from Salmon
# or RSEM above. Columns must be unique and containing those names.
idcolumns <- names(which(sapply(sample_sheet, function(x) {
    n_distinct(x) == length(x) && # column must be unique
        all(names(files) %in% as.character(x)) # must identify each sample
})))

if (length(idcolumns) == 0) {
    message("The sample sheet columns and file names didn't match up.",
            " Something is wrong. Bailing out.")
    message("samples:")
    print(as.data.frame(sample_sheet))
    message("files:")
    print(as.data.frame(files))
    stop("Unable to determine id columns")
}

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
    arrange(across(all_of(idcolumns))) %>%
    # Merge in the extra_coldata gathered above
    left_join(
        extra_coldata,
        by = set_names("idcolumn", idcolumns[[1]])
    )

# Sort array columns to match coldata
for (assay in c("counts", "abundance", "length")) {
    message("Sorting ", assay)
    txi[[assay]] <- txi[[assay]][, coldata[[ idcolumns[1] ]] ]
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

        ## Something inside of tximport seems to reset the log sink on the
        ## second call. Resetting it here:
        sink(logfile)
        sink(logfile, type = "message")
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
            colSums(assay(gse)[mito_genes, ]) %>%
                enframe(name = idcolumns[[1]], value = "count_mito"),
            by = idcolumns[[1]]
        ) %>%
        left_join(
            colSums(assay(gse)[ribo_genes, ]) %>%
                enframe(name = idcolumns[[1]], value = "count_ribo"),
            by = idcolumns[[1]]
        ) %>%
        left_join(
            colSums(assay(gse)[non_coding_genes, ]) %>%
                enframe(name = idcolumns[[1]], value = "count_noncoding"),
            by = idcolumns[[1]]
        ) %>%
        left_join(
            colSums(assay(gse)[coding_genes, ] > 0) %>%
                enframe(name = idcolumns[[1]], value = "num_expr_genes"),
            by = idcolumns[[1]]
        ) %>%
        transmute(
            across(all_of(idcolumns)),
            mapped_read_count = count_total,
            pct_mito = count_mito / count_total * 100,
            pct_ribo = count_ribo / count_total * 100,
            pct_noncoding = count_noncoding / count_total * 100,
            num_expr_genes
        )

    message("8. ----------- Try Generating PCA data ----------")

    tryCatch({
        dds <- DESeqDataSet(gse[coding_genes, ], design = ~ 1)
        dds <- dds[, colSums(counts(dds)) > 0]
        dds <- dds[rowSums(counts(dds)) > 0, ]
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
