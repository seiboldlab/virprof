#!/usr/bin/env Rscript
library(here)
library(optparse)
library(magrittr)
suppressPackageStartupMessages(library(dplyr))
library(readr)
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
library(stringr)
library(openxlsx)
library(fuzzyjoin)
library(tibble)

source(here("R", "virushostdb.R"))

#' List of terms marking known respiratory viruses
respiratory_viruses <- c(
    "Rhinovirus",
    "Coronavirus",
    "Influenza",
    "Metapneumovirus",
    "Enterovirus",
    "Orthopneumovirus",
    "Parainfluenza",
    "Respirovirus",
    "Rubulavirus",
    "Respiratory"
)


#' Write Excel workbook with column widths adjusted
write_xlsx <- function(sheets, file,
                       colWidths = "auto", maxWidth = 80, ...) {
    message("Writing XLSX...")
    message("... output: ", file)
    mw <- getOption("openxslx.maxWidth")
    options(openxlsx.maxWidth = maxWidth)
    write.xlsx(sheets, file, colWidths = colWidths, ...)
    options(openxlsx.maxWidth = mw)
    message("... done")
}


field_name_map <- list(
    "Sample" = "sample",
    "Frequent Words" = "words",
    "Log E-Value" = "log_evalue",
    "Log E-Value (min)" = "min_log_evalue",
    "Log E-Values (min)" = "min_log_evalues",
    "# Contigs" = "n_frag",
    "Reference Accession" = "sacc",
    "Reference Accessions" = "saccs",
    "Reference Title" = "stitle",
    "Reference Taxonomy IDs" = "staxids",
    "% Identity" = "pident",
    "% Identities" = "pidents",
    "# Reads" = "numreads",
    "Taxonomy ID" = "taxid",
    "Taxonomic Name" = "taxname",
    "Taxonomic Names" = "taxnames",
    "Species" = "species",
    "Lineage" = "lineage",
    "Lineage Ranks" = "lineage_ranks",
    "Reference BP covered" = "slen",
    "Positive Samples" = "positive_samples"
)

#' Convert code to natural names
rename_fields <- function(names) {
    for (i in seq_along(field_name_map)) {
        from_name = field_name_map[[i]]
        to_name = names(field_name_map)[[i]]
        match = pmatch(from_name, names)
        if (!is.na(match)) {
            names[[match]] <- to_name
        }
    }
    names
}


#' Concatenate fields
#'
#' Helper for summarize() joining strings
#'
#' @param strings List of strings to combine
#' @param sep Separator (paste collapse)
#' @param split Split each string using this RE before combining
#' @param sort Sort strings before combining
#'
concat <- function(strings, sep=";", split=NULL, sort=TRUE) {
    if (!is.null(split)) {
        strings <- unlist(strsplit(strings, split))
    }
    if (sort) {
        strings <- sort(strings)
    }
    paste(collapse=sep, unique(strings))
}


#' Strip last element of lineage (materialized path)
#'
#' @param lineage Lineage string(s)
#' @param n Number of elements to strip
#'
ancestor <- function(lineage, n=1) {
    for (i in seq(n)) {
        lineage <- sub(";[^;]*$", "", lineage)
    }
    lineage
}


#' Parse command line options
#'
#' Create option structure from command line parameters
#'
#' @param args Arguments to parse (default = commandArgs)
#'
parse_options <- function(args = commandArgs(trailingOnly = TRUE)) {
    option_list <- list(
        make_option(c("--virushostdb"),
                    metavar = "FILE_OR_URL",
                    default = "ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv",
                    help = "Virus Host DB. Can be local file or remote URL."
                    ),
        make_option(c("-f", "--filter"),
                    metavar = "REGEX",
                    default = "^Viruses;",
                    help = "Global lineage filter (default: '%default')"
                    ),
        make_option(c("--filter-host"),
                    metavar = "NAME",
                    default = "Homo sapiens",
                    help = "Select only viruses targeting this host (default: '%default')"
                    ),
        make_option(c("-o", "--out"),
                    metavar = "FILE",
                    default = "result_%s.csv",
                    help = "Output file pattern (default: %default)"
                    ),
        make_option(c("--excel-out"),
                    metavar = "FILE",
                    default = "result.xlsx",
                    help = "Output in Excel format (default: %default)"
                    ),
        make_option(c("-p", "--pattern"),
                    metavar = "REGEX",
                    default = "(.*)\\.virus\\.csv",
                    help = "Input file pattern (default: '%default')"
                    ),
        make_option(c("--min-slen"),
                    metavar = "N",
                    default = 200,
                    help = "Minimum covered subject length (default: %default)"
                    ),
        make_option(c("--min-reads"),
                    metavar = "N",
                    default = 3,
                    help = "Minimum read count (default: %default)"
                    ),
        make_option(c("--merge-samples"),
                    metavar = "REGEX",
                    help = "Rewrite sample name using this regex to merge units processed as samples"),
        make_option(c("--in-list"),
                    metavar = "FILE",
                    help = "File with list of input files (alternative specification)"
                    )
    )

    usage <- "usage: %prog [options] [samples]"

    description <- paste(c(
        "",
        "Gathers classifications from CSV and creates report"
    ))

    opt <- parse_args2(
          OptionParser(
              option_list = option_list,
            usage = usage,
            description = description
        ),
        args = args
    )

    if (length(opt$args) == 0 && is.null(opt$options$in_list)) {
        stop("Need at least one sample argument")
    }

    opt
}

#' Parse file arguments
#'
#' Creates dataframe with `path` filled with each file
#' name. Directories are expanded with option `pattern`. Verifies that
#' all files are readible.
#'
#' @param opt List with $arg parameter pointing to directory or list
#'     of files.
locate_files <- function(opt) {
    message("Locating files from arguments...")

    if (length(opt$args) == 1 && file_test("-d", opt$args[1])) {
        message("... single directory passed, globbing for ", opt$options$pattern)
        ## Just one directory argument, globbing for files
        samples <- tibble(
            path=list.files(
                path=opt$args,
                pattern=opt$options$pattern,
                full.names=TRUE
            ))
        if (nrow(samples) == 0) {
            stop("No samples in directory")
        }
    } else {
        message("... scanning ", length(opt$args), " files passed as arguments")
        ## Individually listed arguments
        samples <- tibble(path=opt$args)
    }
    if (!is.null(opt$options$in_list)) {
        message("... loading files from file ", opt$options$in_list)
        samples <- bind_rows(
            samples,
            read_table(opt$options$in_list,
                       col_names='path',
                       col_types=c(col_character()))
        )
    }
    message("... found ", nrow(samples), " sample data files")

    file_ok <- file_test("-f", samples$path) &
        file.access(samples$path, 4) == 0
    if (any(!file_ok)) {
        message("Unable to read files:")
        print(samples$path[which(!file_ok)])
        stop("Cannot read file(s). Aborting.")
    }
    message("... all files readable")
    message("... done")
    samples %<>%
        mutate(
            unit = sub(paste0(".*/", opt$options$pattern), "\\1", path)
        )
}

#' Load CSVs
#'
#'
load_files <- function(samples, opt) {
    message("Loading files...")
    cov_cols <- cols(
        sample = col_character(),
        words = col_character(),
        log_evalue = col_double(),
        slen = col_integer(),
        n_frag = col_integer(),
        sacc = col_character(),
        stitle = col_character(),
        staxids = col_character(),
        pident = col_double(),
        numreads = col_integer(),
        taxid = col_integer(),
        taxname = col_character(),
        species = col_character(),
        lineage = col_character(),
        lineage_ranks = col_character()
    )
    samples <- samples %>%
        mutate(
            data = map(path, read_csv, col_types=cov_cols)
        ) %>%
        unnest(cols=c(data)) %>%
        select(-path)
    n_samples = length(unique(samples$sample))
    message("... found ", n_samples, " unique sample names")
    if (!is.null(opt$options$merge_samples)) {
        message("... rewriting sample names using regex ", opt$options$merge_samples)
        samples <- samples %>%
            mutate(
                sample = gsub(paste0(opt$options$merge_samples, ".*"),
                              "\\1", sample)
            )
        n_samples = length(unique(samples$sample))
        message("... found ", n_samples, "unique sample names after rewrite")
    }
    message("... found a total of ", nrow(samples), " detected organisms")
    message("... done")
    samples
}

#' Basic hit filtering
filter_hits <- function(samples, opt) {
    message("Filtering accession level bins...")
    message("... minimum `slen`: ", opt$option$min_slen)
    message("... minimum read count: ", opt$option$min_reads)
    message("... lineage regex: ", opt$option$filter)

    samples <- samples %>%
        filter(
            grepl(opt$option$filter, lineage),
            slen >= opt$option$min_slen,
            numreads >= opt$option$min_reads
        )
    message("... remaining detections: ", nrow(samples))
    samples
}

#' Load and add virus to host annotation from VirusHostDB (genome.jp)
filter_vhdb_host <- function(samples, opt) {
    message("Filtering by virual host '", opt$options$filter_host, "' ...")
    samples <- samples %>%
        mutate(
            host = vhdb_get_host(lineage, taxid, url=opt$options$virushostdb)
        ) %>%
        filter(
            grepl(opt$options$filter_host, host)
        )
    message("... remaining detections: ", nrow(samples))
    samples
}

#' Filter known respiratory viruses
filter_respiratory_viruses <- function(samples) {
    message("Filtering with positive list...")
    re <- paste0("(", paste(collapse="|", respiratory_viruses), ")")
    samples <- samples %>%
        filter(
            grepl(re, taxname, ignore.case = TRUE)
        )
    message("... remaining detections: ", nrow(samples))
    samples
}

#' Merge species
merge_species <- function(samples) {
    message("Merging species...")
    samples <- samples %>%
        group_by(sample, species) %>%
        summarize(
            taxid    = concat(taxid),
            taxname  = concat(taxname),
            numreads = sum(numreads),
            min_log_evalue = min(log_evalue),
            pident   = round(sum(pident * slen) / sum(slen), 1),
            slen     = concat(slen),
            n_frag   = sum(n_frag),
            sacc     = concat(sacc),
            staxids  = concat(staxids, split=";"),
            lineage  = concat(lineage, "|"),
            .groups="drop"
        ) %>%
        ungroup()
    message("... remaining detections: ", nrow(samples))
    samples
}

merge_samples <- function(samples) {
    samples %>%
        group_by(sample) %>%
        arrange(desc(numreads)) %>%
        summarize(
            taxnames = paste(collapse=";", taxname),
            numreads = paste(collapse=";", numreads),
            min_log_evalues = paste(collapse=";", min_log_evalue),
            pidents = paste(collapse=";", pident),
            slen = paste(collapse=";", slen),
            n_frag = paste(collapse=";", n_frag),
            saccs = paste(collapse=";", sacc),
            .groups="drop"
        ) %>%
        ungroup()
}


if (!interactive()) {
    opt <- parse_options()

    results <- list()
    summary <- tibble(
        Description=character(),
        Value=numeric(),
        Tab=character()
    )

    ## Determine file locations from arguments
    files <- locate_files(opt)
    summary %<>%
        add_row(
            Value=nrow(files),
            Description="Units (files) processed",
            Tab=""
        )

    ## Load the files
    calls <- load_files(files, opt)
    summary %<>%
        add_row(
            Value=nrow(calls),
            Description="Total detections",
            Tab=""
        )

    ## Filter calls by lineage, length and read count
    filtered_calls <- filter_hits(calls, opt)
    results$Detections <- filtered_calls
    summary %<>%
        add_row(
            Value=nrow(filtered_calls),
            Description="Filtered detections",
            Tab="Detections"
        ) %>%
        add_row(
            Description=paste("... Minimum subject bases covered:", opt$options$min_slen)
        ) %>%
        add_row(
            Description=paste("... Minimum read count:", opt$options$min_reads)
        ) %>%
        add_row(
            Description=paste("... Lineage Regex: ", opt$options$filter)
        )

    ## Filter calls by host
    load_vhdb(url=opt$options$virushostdb)
    host_filtered_calls <- filter_vhdb_host(filtered_calls, opt)
    results$`Host Filtered` = host_filtered_calls
    summary %<>%
        add_row(
            Value=nrow(host_filtered_calls),
            Description="Host filtered detections",
            Tab="Host Filtered"
        ) %>%
        add_row(
            Description=paste("... Host Regex:", opt$options$filter_host)
        )

    ## Count found viruses
    viruses_found <- host_filtered_calls %>%
        group_by(taxname) %>%
        summarize(positive_samples=n(), .groups="drop") %>%
        arrange(taxname)
    results$`Viruses Found` = viruses_found
    message("Found ", nrow(viruses_found), " distinct viruses across all samples")
    summary %<>%
        add_row(
            Value=nrow(viruses_found),
            Description="Distinct Viruses Found",
            Tab="Viruses Found"
        )

    ## Filter respiratory viruses (positive list)
    resp_viruses <- filter_respiratory_viruses(host_filtered_calls)
    results$`Respiratory Viruses` <- resp_viruses
    summary %<>%
        add_row(
            Value=nrow(resp_viruses),
            Description="Respiratory Virus detections",
            Tab="Respiratory Viruses"
        )

    row_per_species <- merge_species(resp_viruses)
    results$`Species` <-  row_per_species
    summary %<>%
        add_row(
            Value=nrow(row_per_species),
            Description="Detections after merging multiple calls per species",
            Tab="Species"
        )

    row_per_sample  <- merge_samples(row_per_species)
    results$`Samples` <- row_per_sample
    summary %<>%
        add_row(
            Value=nrow(row_per_sample),
            Description="Positive Samples",
            Tab="Samples"
        )

    results$Summary <- summary

    for (i in seq_along(results)) {
        results[[i]] %<>% rename_with(rename_fields)
    }

    write_xlsx(rev(results), opt$options$excel)
    for (name in names(results)) {
        fn <- sprintf(opt$options$out, name)
        write_csv(results[[name]], fn)
    }

}
