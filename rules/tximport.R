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
library(tibble)
library(DESeq2)
library(lubridate)

message("2. ----------- Loading files ----------")
message("2.1. ----------- Loading Sample Sheet ----------")
message("Filename = ", snakemake@input$meta)

sample_sheet <- read_tsv(snakemake@input$meta, show_col_types = FALSE) %>%
    select(where(~!all(is.na(.))))

metadata <- list(
    virprof_version = snakemake@params$version,
    pipeline = snakemake@params$label,
    date = now(),
    sample_sheet = sample_sheet
)

message("2.2. ----------- Loading MultiQC Report Data ----------")

# FastQC
for (n in c(1,2)) {
    suffix <- c("", "_1")[n]
    round <- c("raw", "trimmed")[n]
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
        message("FastQC data in ", round, " identified by: ",
                paste(fastqc_idcolumns, collapse = ", "))
        # Extract paths
        fastq_id_to_path <- sample_sheet %>%
            pivot_longer(
                cols = where(~any(str_detect(., "(fastq|fq).gz"),
                           # must return T/F, no NA to where:
                                  na.rm = TRUE )),
                values_to = "fastq_file_path"
            ) %>%
            mutate(
                # Deduce mate from file name
                mate = if_else(str_detect(fastq_file_path, "R1"),
                                     "R1", "R2")
            ) %>%
            # Remove anything we can't uniquely match
            group_by(across(fastqc_idcolumns[1]), mate) %>%
            filter(n() == 1) %>%
            ungroup()
        # Merge this into the fastqc df
        df <- df %>%
            left_join(
                fastq_id_to_path,
                by = c(fastqc_id = fastqc_idcolumns[[1]], "mate")
            )
        # Rename sample column to it's name from the sample sheet
        df <- df %>% rename("{fastqc_idcolumns[1]}" := fastqc_id)

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
    message("3.1. ---------- Checking for empty samples --------")
    files <- snakemake@input$counts
    names(files) <- gsub(".salmon", "", basename(dirname(files)))
    files <- files[order(names(files))]

    no_data <- sapply(files, function(fn) {
        nrow(read_tsv(fn, n_max = 1, guess_max = 1, col_types = "cdddd")) == 0
    })
    if (any(no_data)) {
        message("Warning: excluded ", length(which(no_data)),
                " empty samples:")
        message("  ", paste(names(files[no_data]), collapse = ", "))
        metadata$empty_samples <- names(files[no_data])
    }

    message("3.2. ----------- Loading quant.sf files ----------")
    txi <- tximport(files[!no_data], type="salmon", txOut=TRUE)

    message("3.3. ----------- Loading meta_info.json files ----------")
    salmon_all_meta <-
        files[!no_data] %>%
        gsub("/quant.sf", "/aux_info/meta_info.json", .) %>%
        purrr::map_df(read_json_nolist, simplifyVector = TRUE,
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
        errorfn = paste0(logfile, ".error.csv")
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

    if (length(metadata$empty_samples) > 0) {
        txi$abundance <- rep(0, nrow(txi$abundance)) %>%
            matrix(ncol = length(metadata$empty_samples)) %>%
            set_colnames(metadata$empty_samples) %>%
            cbind(txi$abundance, .)
        txi$abundance <- txi$abundance[,names(files)]
        txi$counts <- rep(0, nrow(txi$counts)) %>%
            matrix(ncol = length(metadata$empty_samples)) %>%
            set_colnames(metadata$empty_samples) %>%
            cbind(txi$counts, .)
        txi$counts <- txi$counts[,names(files)]
        txi$length <- rep(rowMeans(txi$length),
                          length(metadata$empty_samples)) %>%
            matrix(ncol = length(metadata$empty_samples)) %>%
            set_colnames(metadata$empty_samples) %>%
            cbind(txi$length, .)
        txi$length <- txi$length[,names(files)]
    }
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
               ~ paste(as.character(.x), collapse=";")),
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

    message("8. ----------- Collecting mapping metadata ----------")

    mito_genes <- rowData(gse) %>%
        as_tibble() %>%
        filter(str_detect(gene_name, "^MT-")) %>%
        pull(gene_id)

    ribo_genes <-
        rowData(gse) %>%
        as_tibble() %>%
        filter(gene_type %in% c("rRNA", "rRNA_pseudogene")) %>%
        pull(gene_id)

    non_coding_genes <-
        rowData(gse) %>%
        as_tibble() %>%
        filter(
            gene_type != "protein_coding",
            !gene_id %in% c(mito_genes, ribo_genes)
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
            colSums(assay(gse)[mito_genes,]) %>%
                enframe(name = idcolumns[[1]], value = "count_mito"),
            by = idcolumns[[1]]
        ) %>%
        left_join(
            colSums(assay(gse)[ribo_genes,]) %>%
                enframe(name = idcolumns[[1]], value = "count_ribo"),
            by = idcolumns[[1]]
        ) %>%
        left_join(
            colSums(assay(gse)[non_coding_genes,]) %>%
                enframe(name = idcolumns[[1]], value = "count_noncoding"),
            by = idcolumns[[1]]
        ) %>%
        left_join(
            colSums(assay(gse)[coding_genes,] > 0) %>%
                enframe(name = idcolumns[[1]], value = "num_expr_genes"),
            by = idcolumns[[1]]
        ) %>%
        transmute(
            across(all_of(idcolumns[[1]])),
            pct_mito = count_mito / count_total * 100,
            pct_ribo = count_ribo / count_total * 100,
            pct_noncoding = count_noncoding / count_total * 100,
            num_expr_genes
        )

    message("8. ----------- Try Generating PCA data ----------")

    tryCatch({
        dds <- DESeqDataSet(gse[coding_genes,], design = ~ 1)
        dds <- dds[ , colSums(counts(dds)) > 0 ]
        dds <- dds[ rowSums(counts(dds)) > 0, ]
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
