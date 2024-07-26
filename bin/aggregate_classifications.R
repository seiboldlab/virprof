#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)

    library(magrittr)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(tibble)
    library(lubridate)
    library(fuzzyjoin)

    library(openxlsx)

    library(BiocGenerics)
    library(S4Vectors)
    library(stats4)
})

## Find out where we are. Sadly difficult to do.
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

#' Default list of terms marking known respiratory viruses
#' Normally, this will be overridden by a file supplied
#' from the pipeline configuration.
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
                       colWidths = "auto", maxWidth = 80,
                       firstCol = TRUE, firstRow = TRUE,
                       ...) {
    message("Writing XLSX...")
    message("... output: ", file)
    mw <- getOption("openxslx.maxWidth")
    options(openxlsx.maxWidth = maxWidth)
    write.xlsx(sheets, file, colWidths = colWidths,
               firstCol=firstCol, firstRow=firstRow,
               ...)
    options(openxlsx.maxWidth = mw)
    message("... done")
}

#' Translation table for R "code" names to "English"
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
    "% Identity(s)" = "pidents",
    "Read Count" = "numreads",
    "Read Count(s)" = "numreadss",
    "Taxonomy ID" = "taxid",
    "Taxonomy IDs" = "taxids",
    "Taxonomic Name" = "taxname",
    "Taxonomic Name(s)" = "taxnames",
    "Species Name" = "species",
    "Lineage" = "lineage",
    "Lineage(s)" = "lineages",
    "Lineage Ranks" = "lineage_ranks",
    "Scaffold BP" = "bp",
    "Scaffold BP(s)" = "bps",
    "Max. Scaffold BP" = "max_bp",
    "Max. Scaffold BP(s)" = "max_bps",
    "Subject Length" = "slen",
    "% Genome" = "genome_coverage",
    "% Genome(s)" = "genome_coverages",
    "% Aligned" = "contig_coverage",
    "Positive Samples" = "positive_samples",
    "Host" = "host",
    "Known Respiratory Virus" = "respiratory",
    "Sequencing Units" = "unit_names",
    "Unit Count" = "num_units",
    "Trimmed Reads" = "trimmed_read_count",
    "Mapped Reads" = "mapped_read_count",
    "% Homopolymer" = "pcthp",
    "% Homopolymer(s)" = "pcthps",
    "Entropy" = "entropy",
    "Entropy(s)" = "entropies"
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
        make_option(c("-o", "--out-csv"),
                    metavar = "FILE",
                    default = "result_%s.csv",
                    help = "Output file pattern (default: %default)"
                    ),
        make_option(c("--out-excel"),
                    metavar = "FILE",
                    default = "result.xlsx",
                    help = "Output in Excel format (default: %default)"
                    ),
        make_option(c("--out-rds"),
                    metavar = "FILE",
                    default = "result.rds",
                    help = "Output in RDS format (default: %default)"
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
        make_option(c("--min-aligned-bp"),
                    metavar = "N",
                    default = 150,
                    help = "Minimum aligned sequence length (default: %default)"
                    ),
        make_option(c("--min-reads"),
                    metavar = "N",
                    default = 3,
                    help = "Minimum read count (default: %default)"
                    ),
        make_option(c("--min-pident"),
                    metavar = "PERCENT",
                    default = 70,
                    help = "Minimum percent identity (default: %default)"
                    ),
        make_option(c("--max-pcthp"),
                    metavar = "PERCENT",
                    default = 12,
                    help = "Maximum percent homopolymer (default: %default)"
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
                    ),
        make_option(c("--in-fastaqc-list"),
                    metavar = "FILE",
                    help = "File with list of Fasta QC CSV input files"
                    ),
        make_option(c("--in-rnaseq-stats"),
                    metavar = "FILE",
                    help = "Stats summary file (.stats.rds) from rnaseq pipeline"
                    ),
        make_option(c("--in-white-list"),
                    metavar = "FILE",
                    help = "White list for selecting respiratory viruses"
                    ),
        make_option(c("--set-project"),
                    metavar = "STRING",
                    help = "Set project variable for export"
                    ),
        make_option(c("--set-label"),
                    metavar = "STRING",
                    help = "Set label variable for export"
                    ),
        make_option(c("--set-version"),
                    metavar = "STRING",
                    help = "Set version variable for export"
                    ),
        make_option(c("--set-pipeline"),
                    metavar = "STRING",
                    help = "Set pipeline variable for export"
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

#' Load Calls (bins, detections)
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
        pcthp = col_number(),
        entropy1 = col_number(),
        entropy2 = col_number(),
        entropy3 = col_number(),
        entropy4 = col_number(),
        entropy6 = col_number(),
        entropy10 = col_number(),
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
            species,
            taxname,
            words,
            numreads,
            pident,
            n_frag,
            contig_coverage,
            log_evalue,
            stitle,
            everything()
        ) %>%
        mutate(
            log_evalue = round(log_evalue)
        ) %>%
        ## skip the contig level values for the report
        ## we are using the scaffold ones instead
        select(-c(
            numreads,
            pcthp, entropy1, entropy2, entropy3, entropy4,
            entropy6, entropy10,
        ))

    if (all(samples$unit == samples$sample)) {
        samples <- select(samples, -unit)
    }
    message("... done")
    samples
}


#' Load coverage data and add read counts to calls
#'
#' This only grabs "numreads" from the samtools
#' coverage output.
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
        summarize(numreads = sum(numreads), .groups="drop")

    left_join(calls, coverages, by=c("sample", "sacc")) %>%
        relocate(numreads, .after = words)
}


#' Load the scaffold data and add to calls
#'
#' This only grabs the "bp" number
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
        summarize(bp = sum(bp), .groups="drop")

    left_join(calls, scaffolds, by=c("sample", "sacc")) %>%
        relocate(bp, .after = pident)
}

#' Load the FastaQC data and add to calls
#'
#' This only loads the percent homopolymer and average entropy. Values
#' are weighted by fragment length.
load_fastaqc <- function(calls, fastaqc_filelist_file) {
    message("Loading Fasta QC data...")
    fastaqc_cols <- cols(
        "acc" = col_character(),
        len_all = col_integer(),
        len_used = col_integer(),
        entropy1 = col_number(),
        entropy2 = col_number(),
        entropy3 = col_number(),
        entropy4 = col_number(),
        entropy6 = col_number(),
        entropy10 = col_number(),
        frac_hp = col_number()
    )
    fastaqc <- fastaqc_filelist_file %>%
        read_table(col_names = "path", col_types = cols(col_character())) %>%
        mutate(data = map(path, read_csv, col_types = fastaqc_cols)) %>%
        unnest(cols = c(data)) %>%
        select(-path) %>%
        mutate(acc = str_replace(acc, "_pilon$", "")) %>%
        separate(acc, c("sample", "sacc"), sep="\\.") %>%
        mutate(sacc = sub("(\\d)_(\\d{1,3})", "\\1:\\2", sacc)) %>%
        separate(sacc, c("sacc", "fragment"), sep=":", fill = "right") %>%
        select(-fragment) %>%
        group_by(sample, sacc) %>%
        summarize(
            pcthp = round(sum(frac_hp * len_used) / sum(len_used) * 100, 1),
            entropy = round(sum(
                (entropy1 + entropy2 + entropy3 + entropy4
                    + entropy6 + entropy10) * len_used) /
                sum(len_used) / 6, 2)
        )

    left_join(calls, fastaqc, by = c("sample", "sacc")) %>%
        relocate(pcthp, entropy, .after = bp)
}


#' Apply thresholds to collected hits
filter_hits <- function(samples, opt) {
    total <- nrow(samples)
    stats <- list()
    message("Filtering accession level bins...")
    message("... input detection count: ", total)

    message("... minimum `bp`: ", opt$options$min_bp)
    samples <- filter(samples, bp >= opt$options$min_bp)
    stats$min_bp <- total - nrow(samples)
    total <- nrow(samples)
    message("... remaining detections: ", total)

    message("... minimum aligned bp: ", opt$options$min_aligned_bp)
    samples <- filter(samples,
                      bp * contig_coverage / 100 >= opt$options$min_aligned_bp)
    stats$min_aligned_bp <- total - nrow(samples)
    total <- nrow(samples)
    message("... remaining detections: ", total)

    message("... minimum read count: ", opt$option$min_reads)
    samples <- filter(samples, numreads >= opt$option$min_reads)
    stats$min_reads <- total - nrow(samples)
    total <- nrow(samples)
    message("... remaining detections: ", total)

    message("... minimum % blast identity: ", opt$options$min_pident)
    samples <- filter(samples, pident >= opt$options$min_pident)
    stats$min_pident <- total - nrow(samples)
    total <- nrow(samples)
    message("... remaining detections: ", total)

    message("... maximum % homopolymers: ", opt$option$max_pcthp)
    samples <- filter(samples, pcthp <= opt$option$max_pcthp)
    stats$max_pcthp <- total - nrow(samples)
    total <- nrow(samples)
    message("... remaining detections: ", total)

    attr(samples, "stats") <- stats
    samples
}

annotate_viruses <- function(samples, opt) {
    message("Annotating Viruses")
    if (!is.null(opt$options$in_white_list)) {
        message("Loading whitelist from ", opt$options$in_white_list)
        respiratory_viruses <- readLines(opt$options$in_white_list)
    }
    re <- paste0("(", paste(collapse="|", respiratory_viruses), ")")

    samples <- samples %>%
        filter(
            grepl("^Viruses;", lineage)
        ) %>%
        mutate(
            host = vhdb_get_host(lineage, taxid, url=opt$options$virushostdb),
            respiratory = grepl(re, taxname, ignore.case = TRUE)
        ) %>%
        relocate(respiratory, .after = taxname)
    message("... found Viruses: ", nrow(samples))
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
            numreads = sum(numreads),
            # get largest scaffold
            max_bp    = max(bp),
            # calculate weighted % identity
            pident   = round(sum(pident * slen) / sum(slen), 1),
            # calculate weighted % homopolymer
            pcthp = round(sum(pcthp * bp) / sum(bp), 1),
            # calculate weighted entropy
            entropy = round(sum(entropy * bp) / sum(bp), 4),
            # collect each value for these:
            genome_coverages = concat(genome_coverage),
            bps       = concat(bp),
            n_frags   = concat(n_frag),
            # pick the best evalue
            min_log_evalue = min(log_evalue),
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
            numreadss = concat(numreads, " | "),
            max_bps = concat(max_bp, " | "),
            pidents = concat(pident, " | "),
            pcthps = concat(pcthp, " | "),
            entropies = concat(entropy, " | "),
            n_frags = concat(n_frags, " | "),
            genome_coverages = concat(genome_coverages, " | "),
            min_log_evalues = concat(min_log_evalue, " | "),
            saccs = concat(saccs, " | "),
            .groups="drop"
        ) %>%
        ungroup()
}


if (!interactive()) {
    opt <- parse_options()

    results <- list()
    results$`Run Info` <- list(
        virprof_version = opt$options$set_version,
        project = opt$options$set_project,
        label = opt$options$set_label,
        date = now(),
        pipeline = opt$options$set_pipeline,
        respiratory_viruses = respiratory_viruses,
        filter_host = opt$options$filter_host,
        min_bp = opt$options$min_bp,
        min_aligned_bp = opt$options$min_aligned_bp,
        min_reads = opt$options$min_reads,
        min_pident = opt$options$min_pident,
        max_pcthp = opt$options$max_pcthp
    )

    summary <- tibble(
        Description=character(),
        Count=numeric(),
        Tab=character()
    )

    ## Determine file locations from arguments
    files <- locate_files(opt)
    message("... units processed: ", nrow(files))
    summary %<>%
        add_row(
            Count=nrow(files),
            Description="Files processed",
            Tab=""
        )

    ## Load the files
    calls <- load_calls(files, opt)
    summary %<>%
        add_row(
            Count=nrow(calls),
            Description="Raw Detections",
            Tab=""
        )
    if (!is.null(opt$options$in_coverage_list)) {
        calls <- load_coverages(calls, opt$options$in_coverage_list)
    }
    if (!is.null(opt$options$in_scaffold_list)) {
        calls <- load_scaffolds(calls, opt$options$in_scaffold_list)
    }
    if (!is.null(opt$options$in_fastaqc_list)) {
        calls <- load_fastaqc(calls, opt$options$in_fastaqc_list)
    }
    message("Finished loading files")

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
    filter_stats <- attr(filtered_calls, "stats")
    attr(filtered_calls, "stats") <- NULL
    results$Detections <- filtered_calls
    summary %<>%
        add_row(
            Description="Filtered Detections",
            Count=nrow(filtered_calls),
            Tab="Detections"
        ) %>%
        add_row(
            Description=paste("... Minimum scaffold bp:", opt$options$min_bp),
            Count = filter_stats$min_bp
        ) %>%
        add_row(
            Description=paste("... Minimum aligned bp:", opt$options$min_aligned_bp),
            Count = filter_stats$min_aligned_bp
        ) %>%
        add_row(
            Description=paste("... Minimum read count:", opt$options$min_reads),
            Count = filter_stats$min_reads
        ) %>%
        add_row(
            Description=paste("... Minimum % blast identity:", opt$options$min_pident),
            Count = filter_stats$min_pident
        ) %>%
        add_row(
            Description=paste("... Maximum % homopolymers:", opt$options$max_pcthp),
            Count = filter_stats$max_pcthp
        )

    ## Summary taxonomic stats
    species_found <- filtered_calls %>%
        group_by(species) %>%
        summarize(
            n_samples = length(unique(sample)),
            numreads = sum(numreads),
            .groups="drop"
        ) %>%
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

    if (!is.null(opt$options$in_rnaseq_stats)) {
        message("Loading RNA-Seq stats...")
        message("... file=", opt$options$in_rnaseq_stats)
        stats <- readRDS(opt$options$in_rnaseq_stats)
        results$`Run Info`$rnaseq_run_info <- list(
            virprof_version = stats$virprof_version,
            pipeline = stats$pipeline,
            date = as.character(stats$date)
        )

        sample_sheet <- stats$sample_sheet %>%
            left_join(as_tibble(stats$coldata), by = c("unit", "sample")) %>%
            group_by(sample) %>%
            summarize(
                unit_names = paste(unit, collapse=", "),
                num_units = n(),
                trimmed_read_count = sum(num_processed),
                mapped_read_count = sum(num_mapped),
                .groups = "drop"
            ) %>%
            left_join(
                positive_samples %>%
                select(sample, taxnames, numreadss, max_bps,
                       pidents, genome_coverages, pcthps),
                by = "sample"
            )
        results$`Samples` <- sample_sheet
        summary %<>%
            add_row(
                Count=nrow(sample_sheet),
                Description="Analyzed Samples",
                Tab="Samples"
            )
    }

    results$`Summary` <- summary

    if (!is.null(opt$options$out_rds)) {
        message("Writing RDS...")
        saveRDS(results, opt$options$out_rds)
    }

    for (i in seq_along(results)) {
        if (!is.data.frame(results[[i]])) {
            results[[i]] %<>%
                enframe(name="property") %>%
                rowwise() %>%
                mutate(
                    value = paste(collapse=", ", value)
                )
        }
    }

    if (!is.null(opt$options$out_csv)) {
        message("Writing CSVs...")
        for (name in names(results)) {
            fname <- tolower(sub(" ", "_", name))
            fn <- sprintf(opt$options$out_csv, fname)
            message("... ", fn)
            write_csv(results[[name]], fn)
        }
    }

    if (!is.null(opt$options$out_excel)) {
        message("Writing Excel")
        for (i in seq_along(results)) {
            if (!is.data.frame(results[[i]])) {
                results[[i]] <- NULL
            } else {
                results[[i]] %<>% rename_with(rename_fields)
            }
        }

        results$Summary %<>% mutate(
            Tab = if_else(
                Tab=="", "",
                makeHyperlinkString(Tab, text = paste("[", Tab, "]")
            )
        ))
        class(results$Summary$Tab) <- c(class(results$Summary$Tab), "formula")
        results$Summary$"Link to Sheet" <- results$Summary$Tab
        results$Summary$Tab <- NULL
        write_xlsx(rev(results), opt$options$out_excel,
                   zoom=c(1.5, 1, 1, 1, 1, 1)*100)
    }
}
