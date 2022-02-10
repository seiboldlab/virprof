#!/usr/bin/env Rscript
library(optparse)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(grid)
library(ggnewscale)

libdir <- file.path(
    dirname(dirname(
        sub("--file=", "", grep("--file",
                                commandArgs(trailingOnly = FALSE),
                                value = TRUE)
            )
    )),
    "R"
)
source(file.path(libdir, "gene_plot.R"))

parse_options<- function(args = commandArgs(trailingOnly = TRUE)) {
    option_list <- list(
        make_option(c("--input"),
                    metavar = "FILE",
                    help = "Result calls file from blastbin [REQUIRED]"
                    ),
        make_option(c("--input-hits"),
                    metavar = "FILE",
                    help = "Result hits file from blastbin [REQUIRED]"
                    ),
        make_option(c("--input-features"),
                    metavar = "FILE",
                    help = "Result features file from blastbin"
                    ),
        make_option(c("--input-bam"),
                    metavar = "FILELIST",
                    help = "Comma separated list of sorted BAMs mapping reads to contigs (enables coverage plot)"
                    ),
        make_option(c("--scaffold-bam"),
                    metavar = "FILELIST",
                    help = "Comma separated list of sorted BAMs mapping reads to scaffolds (enables reference coverage plot"),
        make_option(c("--output"),
                    metavar = "FILE",
                    help = "Output PDF"
                    ),
        make_option(c("--output-rds"),
                    metavar = "FILE",
                    help = "Output RDS"
                    ),
        make_option(c("--page-width"),
                    metavar = "INCHES",
                    help = "Output PDF width (default: %default)",
                    default=22
                    ),
        make_option(c("--page-height"),
                    metavar = "INCHES",
                    help = "Ootput PDF height (default: %default)",
                    default=18
                    ),
        make_option(c("--plots-per-page"),
                    metavar = "N",
                    help = "Number of plots per page (default: %default)",
                    default=5
                    ),
        make_option(c("--max-plots"),
                    metavar = "N",
                    help = "Maximum number of plots to output (default: all)",
                    default=0),
        make_option(c("--min-slen"),
                    metavar = "LEN",
                    help = "Minimum slen for rendered hits (default: %default)",
                    default = 200),
        make_option(c("--min-reads"),
                    metavar = "N",
                    help = "Minimum read count for rendered hits (default: %default)",
                    default = 3),
        make_option(c("--warn-level"),
                    metavar = "0|1|2",
                    help = "Set level reporting warnings: 0=collected, 1=immediate, 2=exit (default: %default)",
                    default = 1),
        make_option(c("--cache-path"),
                    metavar = "PATH",
                    help = "Data fetched from Entrez will be cached in this location",
                    default = "/tmp/entrez_cache")
    )
    usage <- "usage: %prog [options]"
    description <- paste(c(
        "",
        "Plots alignment of contigs to selected reference sequences"
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



#' Renders per page plots
#'
#' Calls the function plot_item for each item, printing
#' plots_per_page plots at a time (per page).
plot_pages <- function(plots_per_page, items, plot_item) {
    if (length(items) == 0) {
        return()
    }
    npages <- as.integer((length(items) + plots_per_page - 1) / plots_per_page)

    for (page in 1:npages) {
        message("Printing page ", page, "/", npages, "...")
        first <- (page-1) * plots_per_page + 1
        last <- min(length(items), first + plots_per_page  - 1)
        plots <- vector("list", plots_per_page)
        for (i in seq_along(plots)) {
            plots[[i]] <- plot_spacer()
        }
        for (i in first:last) {
            plots[[i - first + 1]] <- plot_item(items[[i]])
        }
        print(wrap_plots(plots, ncol=1))
    }
}

print_pdf <- function(opt, vp) {
    message("Ordering calls ...")
    saccs <- vp@calls %>%
        # Get total read count per species
        group_by(species) %>%
        mutate(reads_per_species = sum(numreads)) %>%
        ungroup() %>%
        # Order by species total read count, then call read count
        arrange(desc(reads_per_species), desc(numreads)) %>%
        pull(sacc)

    if (opt$options$max_plots > 0) {
        message("Plotting only top ", opt$options$max_plots, " calls")
        saccs <- head(saccs, opt$options$max_plots)
    }
    message("Writing output to ", opt$options$output)
    pdf(
        file = opt$options$output,
        width = opt$options$page_width,
        height = opt$options$page_height
    )

    message("Generating plots...")
    plot_pages(opt$options$plots_per_page, saccs, function(acc) {
        plot(vp, accession=acc)
    })
    message("Finished")
    dev.off()
}

if (!interactive()) {
    if (T) {
        opt <- parse_options()
        options(warn=opt$options$warn_level)
    } else {
        opt <- parse_options(args=c(
                                 "--input", "out.virus.csv",
                                 "--input-hits", "out.hits.csv",
                                 "--input-features", "out.features.csv",
                                 "--input-bam", "out.bam",
                                 "--scaffold-bam", "out.scaffold.bam",
                                 "--output", "out.pdf",
                                 "--output-rds", "out.rds",
                                 "--max-plots", "5"))
    }
    vp <- VirProfFromCSV(opt$options$input)
    vp <- place_contigs(vp)
    if (!is.null(opt$options$input_bam)) {
        vp <- coverage_depth(vp, opt$options$input_bam)
    }
    if (!is.null(opt$options$scaffold_bam)) {
        vp <- coverage_depth(vp, opt$options$scaffold_bam, scaffold=TRUE)
    }
    if (!is.null(opt$options$output_rds)) {
        saveRDS(vp, opt$options$output_rds)
    }
    if (!is.null(opt$options$output)) {
        print_pdf(opt, vp)
    }
}
