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

    #### CONFIG  ###
    cat(".")
    heights <- c(
        coverage = 5,
        contigs = 1,
        alignments = 3,
        subject = 1,
        annotations = 1
    )
    ymax <- cumsum(rev(heights))
    ymin <- ymax - rev(heights)
    ymax <- ymax - ymin[["coverage"]]
    ymin <- ymin - ymin[["coverage"]]
    heights <- as.list(heights)
    ymax <- as.list(ymax)
    ymin <- as.list(ymin)
    box_colors <- c(
        "Contig" = "darkblue",
        "Reference" = "darkred",
        "gene" = "orange",
        "product" = "darkorange"
    )
    label_colors <- c(
        "Contig" = "white",
        "Reference" = "white",
        "gene" = "black",
        "product" = "black"
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
    display_slots <-
        bind_rows(
            hits %>% mutate(start=sstart, stop=send) %>% select(start, stop),
            hits %>% mutate(start=cstart, stop=cend) %>% select(start, stop)
        ) %>%
        arrange(start, stop) %>%
        compress_axis(merge_dist=200, spacing = 100)

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

    contig_boxes <- data.frame(
        xmin = contigs$cstart,
        xmax = contigs$cend,
        ymin = ymin$contigs,
        ymax = ymax$contigs - 0.1,
        type = "Contig",
        label = contigs$contig
    )

    subject_boxes<- data.frame(
        xmin = hits$sstart,
        xmax = hits$send,
        ymin = ymin$subject,
        ymax= ymax$subject,
        type = "Reference",
        label = NA
    )

    boxes <- rbind(
        contig_boxes,
        subject_boxes
    )

    polygons <- alignment_boxes

    cat(".")
    p <- ggplot() +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle=90, hjust=1, size=6)
        ) +
        scale_y_continuous() +
        scale_x_continuous(trans = display_slots$trans)  +
        coord_cartesian(
            ylim = c(NA, max(unlist(ymax))),
            xlim = c(min(hits$cstart, hits$sstart, hits$cend, hits$send),
                     max(hits$cstart, hits$sstart, hits$cend, hits$send)),
            )

    cat(".")
    if(!is.null(depths)) {
        label_tpl <- paste0(label_tpl, ", Average Read Depth: {mean_depth}")

        depths <- depths %>%
            rename(qacc=contig) %>%
            inner_join(contigs, by="qacc")

        mean_depth <- round(mean(depths$total), 1)
        sd_depth <- sd(depths$total)
        max_depth <- max(depths$total)
        cap_depth <- min(max_depth, mean_depth + 3 * sd_depth)

        depth_data <- depths %>%
            mutate(
                y = total / cap_depth * heights$coverage,
                x = if_else(flip, cend - pos, cstart + pos),
                swap_if(flip, minus, plus),
                fplus = if_else(total == 0, 0.5, plus / total),
            ) %>%
            select(
                qacc, x, y
            )

        ## Plot depths above contigs
        p <- p +
            stat_summary_bin(
                geom="bar",
                orientation="x",
                fun="median",
                bins=min(2000, nrow(depth_data)),
                data = depth_data,
                aes(x = x, y = y)
            )
    }

    cat(".")
    if (!is.null(feature_tables)) {
        ## Annotate subject sequences
        subject_annotations <- annotate_subjects(
            display_slots$ranges$start,
            display_slots$ranges$end,
            reference$sacc,
            feature_tables
        )
        ## Remove gene for now
        subject_annotations <- subject_annotations %>%
            filter(key == "product") %>%
            mutate(key == "Annotation")

        ## Merge annotations spanning display ranges
        subject_annotations <- subject_annotations %>%
            group_by(gstart, gstop, key, value) %>%
            summarize(
                start=min(start),
                stop=max(stop),
                .groups="drop"
            )
        ## Compute label location
        subject_annotations <- subject_annotations %>%
            mutate(
                trans_s = display_slots$trans$trans(start),
                trans_e = display_slots$trans$trans(stop),
            )
        ## Remove duplicates
        subject_annotations <- subject_annotations %>%
            group_by(start, stop, key, value) %>%
            summarize(.groups="drop")
        ## Stack
        subject_annotations <- subject_annotations %>%
            mutate(
                bin = IRanges::disjointBins(IRanges::IRanges(start, stop))
            )

        bins <- max(subject_annotations$bin)
        y <- ymax$annotations
        #height <- bins ##(ymax$annotations - y)/bins

        annotation_boxes <- data.frame(
            xmin = subject_annotations$start,
            xmax = subject_annotations$stop,
            ymin = y - subject_annotations$bin + 0.9,
            ymax = y - subject_annotations$bin,
            type = subject_annotations$key,
            label = subject_annotations$value
        )
        boxes <- rbind(boxes, annotation_boxes)

        ## Plot annotations
    }

    boxes <- boxes %>%
        mutate(
            xmid = display_slots$trans$inv((
                display_slots$trans$trans(xmin) +
                display_slots$trans$trans(xmax)
            ) / 2),
            ymid = (ymin+ymax)/2
        )

    p <- p +
        new_scale_fill() +
        scale_fill_manual(values = box_colors) +
        ## Boxes
        geom_rect(
            data = boxes,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=type)
        ) +
        new_scale_color() +
        scale_color_manual(values = label_colors) +
        geom_text_repel(
            data = boxes,
            aes(x=xmid, ymid, label=label, color=type),
            size = 3,
            vjust = .5,
            hjust = .5,
            min.segment.length = 0.1,
            point.size = NA,  # don't shift around center
            max.overlaps = 100,
        ) +
        new_scale_fill() +
        scale_fill_continuous(limits = c(70, 100)) +
        new_scale_color() +
        geom_polygon(
            data = polygons,
            aes(x = x, y = y, group = alignment, fill = pident),
            color="lightblue", alpha=.5
        )

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
        depths <- load_coverage(opt$options$input_bam)
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
