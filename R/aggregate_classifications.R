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


#' List of known respiratory viruses
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


#' Blacklist of problematic VHDB entries
vhdb_blacklist <- tribble(
    ~vhdb_virus_name, ~vhdb_host_name,
    "IAS virus", "Homo sapiens"
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


#' Load VirusHostDB
load_vhdb <- function(url = NULL, reload = FALSE) {
    global_vhdb <- mget(".virushostdb", env=globalenv(),
                        ifnotfound=list(NULL))[[1]]
    if (!reload && !is.null(global_vhdb)) {
        return (global_vhdb)
    }
    if (is.null(url)) {
        url <- "ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv"
    }

    message("Loading VirusHostDB from '", url, "' ...")
    col_spec <- cols(
        `virus tax id` = col_integer(),    # NCBI taxid
        `virus name` = col_character(),
        `virus lineage` = col_character(),
        `refseq id` = col_character(),     # RefSeq Accession
        `KEGG GENOME` = col_character(),   # Dblink
        `KEGG DISEASE` = col_character(),
        DISEASE = col_character(),         # Disease Name
        `host tax id` = col_integer(),     # NCBI taxid
        `host name` = col_character(),
        `host lineage` = col_character(),
        pmid = col_character(),            # PubMed ID
        evidence = col_character(),        # source for host info
        ## if host tax id = 1, type of env sample
        `sample type` = col_character(),
        ## NCBI tax id for env sample (animal, plant, ...)
        `source organism` = col_integer()
    )
    vhdb <- read_tsv(url, col_types = col_spec)
    message("... processing")
    vhdb_filtered <- vhdb %>%
        rename_all(~ gsub(" ", "_", .)) %>%
        rename_all(tolower) %>%
        rename_all(~ paste0("vhdb_", .)) %>%
        mutate(
            vhdb_virus_lineage_name =
                paste(sep="; ",
                      vhdb_virus_lineage,
                      vhdb_virus_name)
        ) %>%
        anti_join(vhdb_blacklist,
                  by=c("vhdb_virus_name", "vhdb_host_name"))

    # Remove annoying spec attribute created by read_tsv
    attr(vhdb_filtered, "spec") <- NULL

    assign(".virushostdb", vhdb_filtered, env=globalenv())
    message("... done")
    vhdb_filtered
}


#' Compute list of hosts for each virus usign VHDB
vhdb_get_host <- function(lineages, taxids=NULL, url=NULL,
                          vhdb=load_vhdb(url=url)) {
    message("> Finding hosts using VirusHostDB for ", length(lineages), " hits")
    result <- tibble(
        lineage=lineages,
        host_names=""
    )

    ## If we have TaxIDs, check those first:
    if (!is.null(taxids)) {
        message(">  matching by NCBI tax ID")
        result$vhdb_virus_tax_id <- taxids
        result %<>%
            mutate(
                rowid=row_number()
            ) %>%
            left_join(select(vhdb, vhdb_virus_tax_id, vhdb_host_name),
                      by="vhdb_virus_tax_id") %>%
            group_by(rowid, lineage) %>%
            summarize(
                host_names = paste(collapse="; ", na.omit(vhdb_host_name)),
                .groups="drop"
            ) %>%
            ungroup() %>%
            select(lineage, host_names)
        message(">    found ", length(which(nzchar(result$host_names))))
    }


    matches <- tibble(lineage=character(), regex=character(), res=character())
    tosearch <- result %>%
        filter(!nzchar(host_names)) %>%
        group_by(lineage) %>%
        summarize(
            regex=paste0("^", gsub("([][$^()])", "\\\\\\1", unique(lineage))),
            .groups="drop"
        ) %>%
        ungroup()

    sum_linname <- function(hlin, hnam, vnam) {
        tibble(lineage=hlin, name=hnam) %>%
            filter(
                !is.na(lineage),
                !is.na(name)
            ) %>%
            mutate(
                domain=sub(";.*", "", lineage),
                kingdom=sub("^([^;]*;[^;]*);.*", "\\1", lineage),
                phylum=sub("^([^;]*;[^;]*;[^;]*);.*", "\\1", lineage)
            ) %>%
            group_by(domain) %>%
            mutate(n=n()/length(hlin)) %>%
            filter(n >= 0.1) %>%
            group_by(kingdom) %>%
            mutate(n=n()/length(hlin)) %>%
            filter(n >= 0.1) %>%
            group_by(phylum) %>%
            mutate(n=n()/length(hlin)) %>%
            filter(n >= 0.1) %>%
            pull(name) %>%
            unique() %>%
            sort() %>%
            paste(collapse=";")
    }

    for(i in seq(25)) {
        if (nrow(tosearch) == 0) {
            break;
        }
        if (i==1) {
            message("  matching by lineage prefix")
        }
        message(sprintf("    iteration %i: %i lineages left",
                        i, nrow(tosearch)))
        ## Match lineages
        found <- vhdb %>%
            regex_inner_join(
                tosearch,
                by=c(vhdb_virus_lineage_name="regex")
            ) %>%
            group_by(lineage, regex) %>%
            summarize(
                #res = paste(collapse=";", unique(vhdb_host_name))
                res = sum_linname(vhdb_host_lineage, vhdb_host_name, vhdb_virus_name),
                .groups="drop"
            ) %>%
            ungroup()
        matches <- rbind(matches, found)

        # Remove new matches from tosearch and drop one level lineage
        tosearch <- tosearch %>%
            anti_join(found, by="lineage") %>%  # remove found
            mutate(regex=sub(";[^;]*$", "", regex))  # drop lineage level

        # Re-use matches we already determined
        found <- matches %>%
            select(regex, res) %>%
            distinct() %>%
            inner_join(
                tosearch,
                by="regex"
            ) %>%
            select(lineage, regex, res)

        if (nrow(found) > 0) {
            message("      re-using ", nrow(found), " previous matches")
            matches <- rbind(matches, found)
            tosearch <- tosearch %>%
                anti_join(found, by="lineage")
        }
    }

    if (nrow(tosearch) > 0) {
        message(">  Failed to find ", nrow(tosearh), " lineages:")
        print(data)
    }

    result <- result %>%
        left_join(matches, by="lineage") %>%
        mutate(
            host_names = ifelse(nzchar(host_names), host_names, res)
        ) %>%
        pull(host_names)
    message("> ... done")
    result
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
    results$`Species Level` <-  row_per_species
    summary %<>%
        add_row(
            Value=nrow(row_per_species),
            Description="Detections after merging multiple calls per species",
            Tab="Species Level"
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
