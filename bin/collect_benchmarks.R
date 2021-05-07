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
    "basecov_bedtools" = 2,
    "bin_vp" = 2,
    "blastbin_vp" = 2,
    "blastfilter_vp" = 2,
    "coverage_samtools" = 2,
    "dust_bbmap" = 8,
    "extract_reads" = 4,
    "format_bbmap" = 4,
    "index_bowtie2" = 8,
    "map_bowtie2" = 12,
    "map_hisat2" = 16,
    "markdup_sambamba" = 8,
    "polish_pilon" = 2,
    "scaffold_vp" = 2,
    "sort_bam" = 8,
    "spades" = 48,
    "summarize_vp" = 2,
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
    pdf(opt$options$out, width=8, height=6)
}


message("Preparing plot...")

data <- df %>%
    group_by(stage, stack, sample) %>%
    summarize(s=sum(s), cpu_time=sum(cpu_time), .groups="drop") %>%
    mutate(
        subpipeline = stack_to_subpipeline(stack),
        subpipeline_stage = paste(subpipeline, stage),
        threads = stage_to_threads(stage),
        core_hours = s * threads / 3600 / 2
    ) %>%
    group_by(subpipeline) %>%
    mutate(
        stage_rank = subpipeline_stage %>% rank()
    ) %>%
    filter(
        sample != "ALL"
    ) %>% ungroup()

data %>%
    group_by(sample, subpipeline_stage) %>%
    summarise(core_hours = sum(core_hours), .groups="drop") %>%
    spread(subpipeline_stage, core_hours) %>%
    as.matrix() %>%
    dist() %>%
    hclust()

as.data.frame(data[c(4,5),])


g<-ggplot(data, aes(x=sample, y=core_hours, fill=subpipeline, alpha=stage_rank)) +
    scale_y_continuous(
        name="Core Hours",
        expand = expansion(mult = c(0, 0), add=c(.1,.1))
    ) +
    scale_x_discrete(name="Sample") +
    scale_fill_discrete(name=NULL) +
    scale_alpha_continuous(range = c(1, 1), guide=NULL) +
    geom_col(position="stack") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
    ) +
    ggtitle("Wall-Time by Sample and Step")
if (length(unique(df$sample)) > 20) {
    g <- g + theme(
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank()
             )
}
print(g)

g<-data %>%
    filter(threads > 2) %>%
    mutate(
        efficency = cpu_time / (s * threads),
        stage=factor(stage)
    ) %>%
    ggplot(aes(x=cpu_time, y=s*threads)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color="green") +
    geom_abline(intercept = 0, slope = 2, color="orange") +
    geom_abline(intercept = 0, slope = 4, color="red") +
    geom_label(
        data = data %>%
            filter(threads > 2) %>%
            group_by(stage) %>%
            summarize(threads=max(threads)),
        aes(label=threads),
        x=-Inf, y=Inf, vjust=1, hjust=0,
    ) +
    theme_classic() +
    facet_wrap(~stage, scales = "free") +
    ggtitle("CPU Utilization (wall-clock vs CPU time)") +
    labs(
        x = "CPU seconds used",
        y = "CPU seconda allocated"
    )

print(g)



message("Sample processing core hour stats")
data %>%
    group_by(sample) %>%
    summarize(core_hours = sum(core_hours), .groups="drop") %>%
    summary()


if (!is_null(opt$options$out)) {
    message("Closing file ", opt$options$out)
    dev.off()
}

message("DONE")
