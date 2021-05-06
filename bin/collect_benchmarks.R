#!/usr/bin/env Rscript
library(optparse)
library(tidyverse)

parse_options <- function(args = commandArgs(trailingOnly = TRUE)) {
    option_list <- list(
        make_option(c("--in-folder"),
                    metavar = "DIR",
                    help = "Load data from folder"),
        make_option(c("--in-rds"),
                    metavar = "FILE",
                    help = "Load data from RDS",
                    ),
        make_option(c("--out-pdf"),
                    metavar = "FILE",
                    help = "Generate output PDF"),
        make_option(c("--out-rds"),
                    metavar = "FILE",
                    help = "Write data to RDS"),
        make_option(c("--project"),
                    metavar = "NAME",
                    help = "Filter for data from this project")
    )
    usage <- "usage: %prog [options]"
    description <- paste(c(
        "",
        "Plots timings for pipeline runs"
    ))
    opt <- parse_args2(
        OptionParser(
            option_list = option_list,
            usage = usage,
            description = description
        ),
        args = args
    )
    opt
}

opt <- parse_options()


stack_to_subpipeline <- function(stacks) {
    res <- rep("N/A", length(stacks))
    res[grepl("qc_fastqc", stacks)] <- "QC"
    res[grepl("map_hisat2", stacks)] <- "Deplete Host Reads"
    res[grepl("dust_bbmap", stacks)] <- "Assemble"
    res[grepl("annotate_blast", stacks)] <- "Blast & Bin"
    res[grepl("scaffold_vp", stacks)] <- "Scaffold"
    res[grepl("summarize_vp", stacks)] <- "Report"
    res[grepl("index_bowtie2", stacks)] <- "Coverage"
    factor(res, levels=rev(c("QC", "Deplete Host Reads", "Assemble", "Coverage", "Blast & Bin", "Scaffold", "Report")))
}

thread_map = list(
    "UNKNOWN" = -1,
    "annotate_blast" = 24,
    "basecov_bedtools" = 1,
    "bin_vp" = 1,          
    "blastbin_vp" = 1,
    "blastfilter_vp" = 1,
    "coverage_samtools" = 1,
    "dust_bbmap" = 8,
    "extract_reads" = 4,
    "format_bbmap" = 4,    
    "index_bowtie2" = 8,
    "map_bowtie2" = 12,
    "map_hisat2" = 16,       
    "markdup_sambamba" = 8,
    "polish_pilon" = 1, 
    "scaffold_vp" = 1,
    "sort_bam" = 8,
    "spades" = 48,
    "summarize_vp" = 1,     
    "trim_bbmap" = 16
)

stage_to_threads <- function(stage) {
    unlist(thread_map[match(stage, names(thread_map), 1)])
}


if (!is_null(opt$options$in_folder)) {
    message("Searching for files in '", opt$options$in_folder, "/' ...")
    files <- list.files(path = opt$options$in_folder, pattern = "\\.txt$", full.names = TRUE, recursive=TRUE) %>%
        as_tibble() %>%
        separate(value, c("benchmarks", "stage", "stack", "filename"), sep="/", remove=FALSE) %>%
        mutate(filename = sub(".txt$", "", filename)) %>%
        separate(filename, c("sample", "part"), sep="\\.", fill="right") %>%
        separate(stack, c("project", "stack"), sep="\\.", extra="merge")

    if (!is_null(opt$options$project)) {
        message("Filtering ", nrow(files), " files for project '", opt$options$project, "' ...")
        files <- files %>%
            filter(project == opt$options$project)
    }

    if (nrow(files) == 0) {
        message("No files found!")
        quit()
    }

    message("Loading ", nrow(files), " files...")
    df <- files %>%
        mutate(data = map(value, read_tsv, na=c("-"), col_types="dtnnnnnnnn")) %>%
        unnest(cols=c(data)) %>%
        select(-benchmarks, -value)
}

if (!is_null(opt$options$in_rds)) {
    message("Loading '", opt$options$in_rds, "' ...")
    df2 <- readRDS(file = opt$options$in_rds)
    if ("df" %in% ls()) {
        message("Appending newly loaded files to data from RDS")
        df <- bind_rows(df2, df)
    } else {
        df <- df2
    }

    if (!is_null(opt$options$project)) {
        message("Filtering ", nrow(df), " row for project '", opt$options$project, "' ...")
        df <- df %>%
            filter(project == opt$options$project)
    }
}



if (!is_null(opt$options$out_rds)) {
    message("Saving data to '", opt$options$out_rds, "' ...")
    saveRDS(object = df, file = opt$options$out_rds)
}
    

if (!is_null(opt$options$out_pdf)) {
    message("Printing to ", opt$options$out_pdf)
    pdf(opt$options$out)
}

message("Preparing plot...")
df %>%
    group_by(stage, stack, sample) %>%
    summarize(s=sum(s), .groups="drop") %>%
    mutate(
        subpipeline = stack_to_subpipeline(stack),
        threads = stage_to_threads(stage),
        s = s * threads
    ) %>%
    mutate(
        sample = sub("[0-9]*$", "", sample),
        sample = sub("ALL", "Collect", sample)
    ) %>%
    ggplot(aes(x=sample, y=s, fill=subpipeline)) +
    geom_col() +
    theme_classic()


if (!is_null(opt$options$out)) {
    message("Closing file ", opt$options$out)
    dev.off()
}

message("DONE")
