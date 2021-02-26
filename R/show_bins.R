#!/usr/bin/env Rscript
library(here)
library(optparse)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(grid)
source(here("R", "gene_plot.R"))

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
        make_option(c("--output"),
                    metavar = "FILE",
                    help = "Output PDF (default: %default)",
                    default = "out.pdf"
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


#' Get depth from BAM file
coverage_depth <- function(fname) {
    run_bedtools <- function(strand=NULL, split=FALSE, fragment=FALSE) {
        message(
            "Running bedtools genomecov on ", fname,
            if (!is.null(strand)) paste(" counting", strand, "strand"),
            if (split) " not counting gaps",
            if (fragment) " counting whole fragment",
            "..."
        )
        proc <- pipe(paste(
            "bedtools", "genomecov",
            "-ibam", fname,
            "-d",
            if (is.null(strand)) "" else paste("-du -strand", strand),
            if (split) "-split" else "",
            if (fragment) "-fs" else ""
        ))
        read_tsv(proc,
                 col_types = "cii",
                 col_names = c("contig", "pos", "depth"),
                 skip = 1)
    }

    plus <- run_bedtools(strand="+", split=TRUE) %>% rename(plus=depth)
    minus <- run_bedtools(strand="-", split=TRUE) %>% rename(minus=depth)
    result <- full_join(plus, minus, by=c("contig", "pos"))
    result$total = result$plus + result$minus
    result
}

#' Linewrap lineage string
break_lineage <- function(lineage, maxlen=200, insert="\n") {
    lineage %>%
        gsub(" ", "_", .) %>%
        gsub(";_", " ", .) %>%
        strwrap(maxlen) %>%
        gsub(" ", ";_", .) %>%
        gsub("_", " ", .) %>%
        paste(collapse=paste0(";", insert))
}

#' The actual plot function
plot_ranges <- function(reference, hits, depths, feature_tables) {
    message("Plotting ", reference$sacc)
    cat(".")
    heights <- c(
        labels = 1.5,
        coverage = 4,
        contigs = 1,
        contig_hits = .5,
        alignments = 4,
        subject = 1,
        annotations = 1
    )
    ymax <- cumsum(rev(heights))
    ymin <- ymax - rev(heights)
    heights <- as.list(heights)
    ymax <- as.list(ymax)
    ymin <- as.list(ymin)
    fill <- list(
        contigs     = "darkblue",
        contig_hits = "darkgreen",
        alignments  = "green",
        subject     = "darkred",
        annotations = "orange"
    )

    label_tpl <- paste0(
        "Keywords: {reference$words}\n",
        "{reference$sacc}: \"{reference$stitle}\"\n",
        "{lineage}\n",
        "Read Count: {reference$numreads}, ",
        "log E-Value: {reference$log_evalue}, ",
        "Identity: {reference$pident}%, ",
        "Genome Coverage {reference$slen}bp ({reference$genome_coverage}%)"
    )
    ylabel_tpl <- paste0(
        "{reference$species}\n",
        "{reference$sacc}"
    )
    lineage <- break_lineage(reference$lineage)

    cat(".")
    ## Calculate alignment positions for hits
    hits <- hits %>%
        mutate(flip = sstart > send) %>%
        select(sacc, qacc, qstart, qend, sstart, send, cstart, cend, flip, contig, pident) %>%
        mutate(
            aleft  = cstart + qstart - 1,
            aright = cstart + qend - 1
        )


    cat(".")
    ## Get unique contigs from hits
    contigs <- hits %>%
        group_by(qacc, flip, cstart, cend, contig) %>%
        summarize(.groups="drop")

    cat(".")
    ## Make compressed X scale
    xtrans <-
        bind_rows(
            hits %>% mutate(start=sstart, stop=send) %>% select(start, stop),
            hits %>% mutate(start=cstart, stop=cend) %>% select(start, stop)
        ) %>%
        arrange(start, stop) %>%
        compress_axis()

    cat(".")
    ## Setup corners for trapezoid showing alignment mapping
    alignment_boxes <- bind_rows(
        hits %>% mutate(y = ymax$alignments, x = aleft),
        hits %>% mutate(y = ymax$alignments, x = aright),
        hits %>% mutate(y = ymin$alignments, x = send),
        hits %>% mutate(y = ymin$alignments, x = sstart)
    ) %>%
        mutate(
            alignment = paste(qacc, aleft, aright, send, sstart)
        )

    cat(".")
    p <- ggplot() +
        scale_y_continuous(limits = c(0, 16)) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle=90, hjust=1, size=6),
            axis.text.y = element_blank()
        ) +
        scale_x_continuous(trans = xtrans)  +
        ## Contigs
        geom_rect(
            data = contigs,
            aes(xmin = cstart, xmax = cend),
            ymin = ymin$contigs, ymax = ymax$contigs,
            fill = fill$contigs
        ) +
        ## Alignments on contigs
        geom_rect(
            data = hits,
            aes(xmin = aleft, xmax = aright),
            ymin = ymin$contig_hits, ymax = ymax$contig_hits,
            fill = fill$contig_hits
        ) +
        ## Contig labels
        geom_text_repel(
            data = contigs,
            aes(x = (cstart+cend)/2, label = contig),
            y = ymin$labels,
            size = 2,
            nudge_y = 3,
            direction  = "both",
            angle = 0,
            vjust = 0,
            segment.size = 0.5
        ) +
        ## Plot Subject sequence pieces
        geom_rect(
            data = hits,
            aes(xmin=sstart, xmax=send),
            ymin=ymin$subject, ymax=ymax$subject, fill=fill$subject
        ) +
        ## Plot Alignment connection
        geom_polygon(
            data = alignment_boxes,
            aes(x = x, y = y, group = alignment, fill = pident),
            color=fill$alignments, alpha=.5
        ) +
        scale_fill_continuous(limits = c(70, 100))


    cat(".")
    if(!is.null(depths)) {
        label_tpl <- paste0(label_tpl, ", Average Read Depth: {mean_depth}")

        depths <- depths %>%
            rename(qacc=contig) %>%
            inner_join(contigs, by="qacc")

        mean_depth <- round(mean(depths$total), 1)
        max_depth <- max(depths$total)
        norm_depth <- pmin(mean_depth*5, max_depth)

        depth_data <- depths %>%
            mutate(
                x = if_else(flip, cend - pos, cstart + pos),
                fplus = if_else(total == 0, 0.5, if_else(flip, minus, plus) / total),
                total = total / norm_depth,
                plus = plus / norm_depth,
                minus = minus / norm_depth,
                swap_if(flip, minus, plus),
            ) %>%
            gather(
                total, plus, minus, fplus, key="sense", value="y"
            ) %>%
            select(
                qacc, sense, x, y
            )

        ## Plot depths above contigs
        p <- p +
            geom_line(
                data = depth_data,
                aes(
                    x = x,
                    y = y * heights$coverage + ymin$coverage,
                    color=sense,
                    group=interaction(sense,qacc)
                )
            )
    }

    cat(".")
    if (!is.null(feature_tables)) {
        ## Annotate subject sequences
        subject_annotations <- annotate_subjects(
            hits$sstart, hits$send,
            reference$sacc,
            feature_tables
        )
        ## Plot annotations
        p <- p +
            geom_rect(
                data = subject_annotations,
                aes(xmin=start, xmax=stop),
                ymin=ymin$annotations, ymax=ymax$annotations, fill=fill$annotations
            ) +
            geom_text(
                data = subject_annotations,
                aes(x = (start+stop)/2, label = name),
                y = (ymin$annotations + ymax$annotations)/2,
                hjust = .5,
                vjust = .5,
                size=3
            )
    }

    cat(".")
    ## Add labels
    p <- p +
        xlab(str_glue(label_tpl)) +
        ylab(str_glue(ylabel_tpl))

    cat(".")
    message(" [DONE]")
    p
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
        plots <- vector("list", last - first + 1)
        for (i in first:last) {
            plots[[i - first + 1]] <- plot_item(items[[i]])
        }
        print(wrap_plots(plots, ncol=1))
    }
}

run <- function() {
    message("Loading calls...")
    calls <- read_csv(opt$options$input, col_types = cols())
    message("... ", nrow(calls), " calls found")
    calls <- filter(calls, slen >= opt$options$min_slen)
    message("... ", nrow(calls), " calls left after removing calls with <",
            opt$options$min_slen, " bp on subject")
    calls <- filter(calls, numreads >= opt$options$min_reads)
    message("... ", nrow(calls), " calls left after removing calls with <",
            opt$options$min_reads, " mapped reads")
    if (nrow(calls) == 0) {
        message("No calls left, exiting")
        return(NULL)
    }
    calls <- calls %>%
        mutate(
            genome_coverage = if_else(genome_size > 0,
                                      round(slen / genome_size * 100, 1),
                                      NA_real_)
        )

    message("Loading alignments...")
    alignments <- read_csv(opt$options$input_hits, col_types = cols())
    message("... ", nrow(alignments), " alignments found")
    alignments <- filter(alignments, sacc %in% calls$sacc)
    message("... ", nrow(alignments), " alignments matching filtered calls")
    alignments <- alignments %>%
        mutate(
            ## Abbreviate spades contig names
            contig = sub("^NODE_([0-9]+)_.*", "\\1", qacc),
            ## Mark reversed contigs
            reversed = sstart > send,
            ## Alwaus have sstart > send
            swap_if(reversed, qstart, qend),
            swap_if(reversed, sstart, send)
        ) %>%
        ## Check if we should flip contig to have majority alignments lined up
        group_by(qacc) %>%
        mutate(flip = as.logical(median(reversed))) %>%
        ungroup() %>%
        ## Flip contigs
        mutate(
            qstart = if_else(flip, qlen - qstart + 1, qstart),
            qend = if_else(flip, qlen - qend + 1, qend)
        ) %>%
        select(-reversed)

    if (is.null(opt$options$input_features)) {
        features_tables <- NULL
    } else {
        message("Loading feature tables...")
        feature_tables <- read_csv(opt$options$input_features, col_types = cols())
        message("... ", nrow(feature_tables), " feature annotations found")
    }

    message("Placing contigs ...")
    contigs <- place_contigs(alignments)  %>% select(sacc, qacc, cstart, cend)
    alignments <- alignments %>% left_join(contigs, by=c("sacc", "qacc"))

    message("Ordering calls ...")
    saccs <- calls %>%
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

    if (!is.null(opt$options$input_bam)) {
        message("Getting coverages...")
        depths <- opt$options$input_bam %>%
            strsplit(",") %>%
            unlist() %>%
            map_dfr(coverage_depth) %>%
            group_by(contig, pos) %>%
            summarize_all(sum)
    } else {
        depths <- NULL
    }

    message("Generating plots...")
    plot_pages(opt$options$plots_per_page, saccs, function(acc) {
        plot_ranges(
            calls %>% filter(sacc == acc),  # per call data
            alignments %>% filter(sacc == acc),  # per alignment data
            depths,  # coverage depth data
            feature_tables %>% filter(acc == acc)
        )
    })
    message("Finished")
}

if (!interactive()) {
    if (T) {
        opt <- parse_options()
        options(warn=opt$options$warn_level)
    } else {
        opt <- parse_options(args=c("--input", "out.csv", "--input-bam", "out.bam", "--output", "out.pdf"))

        opt <- parse_options(args=c("--input", "out.csv", "--input-bam", "out.bam", "--output", "out.pdf",
                                    "--max-plots", "5"))
    }

    message("Writing output to ", opt$options$output)
    pdf(
        file = opt$options$output,
        width = opt$options$page_width,
        height = opt$options$page_height
    )
    run()
    dev.off()
}

