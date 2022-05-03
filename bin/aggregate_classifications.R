#!/usr/bin/env Rscript
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

libdir <- file.path(
    dirname(dirname(
        sub("--file=", "", grep("--file",
                                commandArgs(trailingOnly = FALSE),
                                value = TRUE)
            )
    )),
    "R"
)
source(file.path(libdir, "virushostdb.R"))

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
    "Respiratory",
    "Bocavirus",
    "Adenovirus"
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
    "Unit" = "unit",
    "Sample" = "sample",
    "Sample(s)" = "samples",
    "Sample Count" = "n_samples",
    "Title Summary" = "words",
    "Log E-Value" = "log_evalue",
    "Log E-Value (min)" = "min_log_evalue",
    "Log E-Value(s) (min)" = "min_log_evalues",
    "Contig Count" = "n_frag",
    "Contig Count(s)" = "n_frags",
    "Reference Accession" = "sacc",
    "Reference Accession(s)" = "saccs",
    "Reference Title" = "stitle",
    "Reference Taxonomy IDs" = "staxids",
    "Reference Genome Size" = "genome_size",
    "Reference Size Source" = "genome_size_source",
    "% Identity" = "pident",
    "% Identities" = "pidents",
    "(old) Read Count" = "numreads",
    "(old) Read Count(s)" = "numreadss",
    "Read Count" = "numreads2",
    "Read Count(s)" = "numreadss2",
    "Taxonomy ID" = "taxid",
    "Taxonomy IDs" = "taxids",
    "Taxonomic Name" = "taxname",
    "Taxonomic Name(s)" = "taxnames",
    "Species Name" = "species",
    "Lineage" = "lineage",
    "Lineage(s)" = "lineage",
    "Lineage Ranks" = "lineage_ranks",
    "Scaffold BP" = "bp",
    "Scaffold BP(s)" = "bps",
    "% Genome" = "genome_coverage",
    "% Genome(s)" = "genome_coverages",
    "% Aligned" = "contig_coverage",
    "Positive Samples" = "positive_samples",
    "Host" = "host",
    "Known Respiratory Virus" = "respiratory"
)

#' Convert code to natural names
rename_fields <- function(names) {
    for (i in seq_along(field_name_map)) {
        from_name <- field_name_map[[i]]
        to_name <- names(field_name_map)[[i]]
        if (!nzchar(to_name)) {
            message("Can't rename ", from_name, " to empty string")
        } else {
            idx <- match(from_name, names)
            if (!is.na(idx)) {
                names[[idx]] <- to_name
            }
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
concat <- function(strings, sep=";", split=NULL, sort=FALSE, unique=FALSE) {
    if (!is.null(split)) {
        strings <- unlist(strsplit(strings, split))
    }
    if (sort) {
        strings <- sort(strings)
    }
    if (unique) {
        strings <- unique(strings)
    }
    paste(collapse=sep, strings)
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
        make_option(c("--min-bp"),
                    metavar = "N",
                    default = 200,
                    help = "Minimum sequence length (default: %default)"
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
                    ),
        make_option(c("--in-coverage-list"),
                    metavar = "FILE",
                    help = "File with list of coverage input files"
                    ),
        make_option(c("--in-scaffold-list"),
                    metavar = "FILE",
                    help = "File with list of scaffold input files"
                    ),
        make_option(c("--in-basecov-list"),
                    metavar = "FILE",
                    help = "File with list of base coverage input files"
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
load_calls <- function(samples, opt) {
    message("Loading calls...")
    call_cols <- cols(
        sample = col_character(),
        words = col_character(),
        log_evalue = col_double(),
        slen = col_integer(),
        n_frag = col_integer(),
        contig_coverage = col_double(),
        sacc = col_character(),
        stitle = col_character(),
        staxids = col_character(),
        pident = col_double(),
        numreads = col_integer(),
        taxid = col_integer(),
        genome_size = col_number(),
        taxname = col_character(),
        species = col_character(),
        lineage = col_character(),
        lineage_ranks = col_character()
    )
    samples <- samples %>%
        mutate(
            data = map(path, read_csv, col_types=call_cols)
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
    samples <- samples %>%
        select(
            unit, sample,
            words,
            numreads,
            log_evalue, pident, n_frag, contig_coverage,
            species, taxname, stitle,
            everything()
        )

    if (all(samples$unit == samples$sample)) {
        samples <- select(samples, -unit)
    }
    message("... done")
    samples
}

load_coverages <- function(calls, coverage_filelist_file) {
    message("Loading coverage files...")
    call_cols <- cols(
        "#rname" = col_character(),
        "startpos" = col_number(),
        "endpos" = col_number(),
        "numreads" = col_integer(),
        "covbases" = col_integer(),
        "coverage" = col_double(),
        "meandepth"  = col_double(),
        "meanbaseq" = col_double(),
        "meanmapq" = col_double()
    )
    coverages <- coverage_filelist_file %>%
        read_table(
            col_names='path',
            col_types=c(col_character())
        ) %>%
        mutate(
            data = map(path, read_tsv, col_types=call_cols)
        ) %>%
        unnest(cols=c(data)) %>%
        select(-path) %>%
        dplyr::rename(rname = "#rname") %>%
        mutate(rname = gsub("_pilon$", "", rname)) %>%
        separate(rname, c("sample", "sacc"), sep="\\.") %>%
        mutate(sacc = sub("(\\d)_(\\d{1,3})", "\\1:\\2", sacc)) %>%
        separate(sacc, c("sacc", "fragment"), sep=":", fill = "right") %>%
        group_by(sample, sacc) %>%
        summarize(numreads2 = sum(numreads), .groups="drop")

    left_join(calls, coverages, by=c("sample", "sacc")) %>%
        relocate(
            numreads2,
            .after = words
        )
}

load_scaffolds <- function(calls, scaffold_filelist_file) {
    message("Loading scaffold files...")
    scaffold_cols <- cols(
        "acc" = col_character(),
        "bin" = col_character(),
        "sstart" = col_integer(),
        "send" = col_integer(),
        "qacc" = col_character(),
        "qstart" = col_integer(),
        "qend" = col_integer(),
        "qlen" = col_integer(),
        "reversed" = col_logical(),
        "bp" = col_integer(),
        "scaffold" = col_character()
    )
    scaffolds <- scaffold_filelist_file %>%
        read_table(
            col_names="path",
            col_types=c(col_character())
        ) %>%
        mutate(
            data = map(path, read_csv, col_types=scaffold_cols)
        ) %>%
        unnest(cols=c(data)) %>%
        select(-path) %>%
        separate("acc", c("sample", "sacc"), sep="\\.", remove = TRUE) %>%
        mutate(sacc = sub("(\\d)_(\\d{1,3})", "\\1:\\2", sacc)) %>%
        separate(sacc, c("sacc", "fragment"), sep=":", fill = "right") %>%
        group_by(sample, sacc, sstart, send) %>%
        mutate(bp = mean(bp)/length(bp)) %>%
        group_by(sample, sacc) %>%
        summarize(
            bp = sum(bp)
        )

    left_join(calls, scaffolds, by=c("sample", "sacc")) %>%
        relocate(
            bp,
            .after = pident
        )
}

#' Basic hit filtering
filter_hits <- function(samples, opt) {
    message("Filtering accession level bins...")
    message("... minimum `bp`: ", opt$option$min_bp)
    message("... minimum read count: ", opt$option$min_reads)

    samples <- samples %>%
        filter(
            bp >= opt$option$min_bp,
            numreads >= opt$option$min_reads
        )
    message("... remaining detections: ", nrow(samples))
    samples
}

annotate_viruses <- function(samples, opt) {
    message("Annotating Viruses")
    re <- paste0("(", paste(collapse="|", respiratory_viruses), ")")

    samples <- samples %>%
        filter(
            grepl("^Viruses;", lineage)
        ) %>%
        mutate(
            host = vhdb_get_host(lineage, taxid, url=opt$options$virushostdb),
            respiratory = grepl(re, taxname, ignore.case = TRUE)
        )
    message("... found Viruses: ", nrow(samples))
    samples
}

#' Filter known respiratory viruses
filter_respiratory_viruses <- function(samples) {
    message("Filtering with positive list...")
    message("... remaining detections: ", nrow(samples))
    samples
}

#' Merge species
merge_species <- function(samples) {
    message("Merging species...")
    samples <- samples %>%
        group_by(sample, species) %>%
        summarize(
            # minimize taxononomic names list
            taxnames  = concat(taxname, ", ", sort=TRUE, unique=TRUE),
            # sum up the read count
            numreads2 = sum(numreads2),
            numreads = sum(numreads),
            # pick the best evalue
            min_log_evalue = min(log_evalue),
            # calculate weighted % identity
            pident   = round(sum(pident * slen) / sum(slen), 1),
            # collect each value for these:
            genome_coverages = concat(genome_coverage),
            bps       = concat(bp),
            n_frags   = concat(n_frag),
            saccs     = concat(sacc),
            # summarize these:
            staxids   = concat(staxids, ", ", split=";", sort=TRUE, unique=TRUE),
            taxids    = concat(taxid, ", ", sort=TRUE, unique=TRUE),
            # should be just one anyway:
            lineages  = concat(lineage, ", ", sort=TRUE, unique=TRUE),
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
            taxnames = concat(taxnames, " | "),
            numreadss2 = concat(numreads2, " | "),
            numreadss = concat(numreads, " | "),
            min_log_evalues = concat(min_log_evalue, " | "),
            pidents = concat(pident, " | "),
            genome_coverages = concat(genome_coverages, " | "),
            bps = concat(bps, " | "),
            n_frags = concat(n_frags, " | "),
            saccs = concat(saccs, " | "),
            .groups="drop"
        ) %>%
        ungroup()
}


if (!interactive()) {
    opt <- parse_options()

    results <- list()
    summary <- tibble(
        Description=character(),
        Count=numeric(),
        Tab=character()
    )

    ## Determine file locations from arguments
    files <- locate_files(opt)
    summary %<>%
        add_row(
            Count=nrow(files),
            Description="Units (files) processed",
            Tab=""
        )

    ## Load the files
    calls <- load_calls(files, opt)
    summary %<>%
        add_row(
            Count=nrow(calls),
            Description="Total detections",
            Tab=""
        )
    if (!is.null(opt$options$in_coverage_list)) {
        calls <- load_coverages(calls, opt$options$in_coverage_list)
    }
    if (!is.null(opt$options$in_scaffold_list)) {
        calls <- load_scaffolds(calls, opt$options$in_scaffold_list)
    }

    ## Compute Genome Coverage
    calls <- calls %>%
        mutate(
            genome_coverage = if_else(
                genome_size > 0,
                round(bp / genome_size * 100),
                NA_real_
            )
        ) %>%
        relocate(
            genome_coverage,
            .before=bp
        )

    ## Filter calls and read count
    filtered_calls <- filter_hits(calls, opt)
    results$Detections <- filtered_calls
    summary %<>%
        add_row(
            Count=nrow(filtered_calls),
            Description="Filtered detections",
            Tab="Detections"
        ) %>%
        add_row(
            Description=paste("... Minimum scaffold bp:", opt$options$min_bp)
        ) %>%
        add_row(
            Description=paste("... Minimum read count:", opt$options$min_reads)
        )

    ## Summary taxonomic stats
    species_found <- filtered_calls %>%
        group_by(species) %>%
        summarize(n_samples = length(unique(sample)), .groups="drop") %>%
        arrange(species)
    results$`Species Found` = species_found
    message("Found ", nrow(species_found), " distinct species across all samples")
    summary %<>%
        add_row(
            Count=nrow(species_found),
            Description="Per species positive sample counts",
            Tab="Species Found"
        )

    ## Filter & Annotate Viruses
    load_vhdb(url=opt$options$virushostdb)
    virus_calls <- annotate_viruses(filtered_calls, opt)
    results$`Virus Detections` <- virus_calls
    summary %<>%
        add_row(
            Count=nrow(virus_calls),
            Description="Annotated virus detections",
            Tab="Virus Detections"
        )

    host_respvir_calls <- virus_calls %>%
        filter(
            grepl(opt$options$filter_host, host),
            respiratory == TRUE
        ) %>%
        merge_species()
    results$`Respiratory Viruses` <- host_respvir_calls
    summary %<>%
        add_row(
            Count=nrow(host_respvir_calls),
            Description=paste0("Respiratory Viruses (", opt$options$filter_host, ")"),
            Tab="Respiratory Viruses"
        )

    positive_samples <- host_respvir_calls %>%
        merge_samples()
    results$`Positive Samples` <- positive_samples
    summary %<>%
        add_row(
            Count=nrow(positive_samples),
            Description="Positive Samples",
            Tab="Positive Samples"
        )

    results$`Summary` <- summary

    for (i in seq_along(results)) {
        results[[i]] %<>% rename_with(rename_fields)
    }

    write_xlsx(rev(results), opt$options$excel)
    message("Writing CSVs...")
    for (name in names(results)) {
        fname <- tolower(sub(" ", "_", name))
        fn <- sprintf(opt$options$out, fname)
        message("... ", fn)
        write_csv(results[[name]], fn)
    }

}
