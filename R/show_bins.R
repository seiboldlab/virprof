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
                    help = "Result file from blastbin [REQUIRED]"
                    ),
        make_option(c("--input-bam"),
                    metavar = "FILE",
                    help = "Sorted BAM mapping reads to contigs (enabled coverage plot)"
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
                    default = "/tmp/entrez_cache"),
        make_option(c("--no-annotate"),
                    dest = "do_annotate", action = "store_false", default = TRUE,
                    help = "Disable annotation layer")
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

#' Parses CSV from python blastbin
#'
#' The CSV contains one row per bin, combining multiple contigs and
#' alignments. For plotting, we split this here into multiple data
#' rows again - one per alignment. The fields qaccs, qranges, sranges,
#' reversed, pidents and numreadss each contain semi-colon separated
#' per contig data. The qranges and sranges also contain
#' dash-separated start/end positions. Lastly, we undo the sorting of
#' start/stop positions and flip contigs for which the majority of
#' alignments are reversed.
#'
#' @param data Data.frame containing data from blastbin output
#' @return data.frame containing alignments
#'   qstart, qstop 1-indexed first and last position of alignment in query
#'   sstart, sstop 1-indexed first and last position of alignment in subject
#'   flip TRUE if the contig coordinates (qstart, qstop) have been flipped
#'   qaccs Query accession
#'   sacc  Subject accession
#'   pidents Percent Identity
#'   numreadss Count of reads mapped to contig
#'   qlen Length of query/contig
#'
parse_blastbins <- function(data) {
    data %>%
        # Split alignments from row-per-bin data
        separate_rows(qranges, sranges, reversed, qaccs, pidents, numreadss, sep=";") %>%
        mutate(
            pidents=as.numeric(pidents),
            numreadss=as.numeric(numreadss),
            reversed=as.logical(reversed=="T"),
            ## FIXME: The contig length should come from python
            qlen=as.integer(sub("NODE_.*_length_([0-9]*)_.*", "\\1", qaccs)),
            contig=as.integer(sub("NODE_(.*)_length_([0-9]*)_.*", "\\1", qaccs))
        ) %>%
        separate(qranges, c("qstart", "qstop"), convert = TRUE) %>%
        separate(sranges, c("sstart", "sstop"), convert = TRUE) %>%
        # Restore qstart/qstop order
        mutate(
            swap_if(reversed, qstart, qstop)
        ) %>%
        ## Flip entire contig if majority of alignments are reversed
        group_by(qaccs) %>%
        mutate(
            flip = as.logical(median(reversed))
        ) %>%
        ungroup() %>%
        mutate(
            qstart = if_else(flip, as.integer(qlen - qstart + 1), qstart),
            qstop  = if_else(flip, as.integer(qlen - qstop  + 1), qstop)
        ) %>%
        arrange(-log_evalue, sacc, sstart) %>%
        ungroup()
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
        "keywords: {reference$words}\n",
        "{reference$sacc}: \"{reference$stitle}\"\n",
        "{lineage}\n",
        "nreads={reference$numreads}, ",
        "log(evalue)={reference$log_evalue}, ",
        "id={reference$pident}%, ",
        "qlen={reference$qlen}, ",
        "slen={reference$slen}"
    )
    ylabel_tpl <- paste0(
        "{reference$species}\n",
        "{reference$sacc}"
    )
    lineage <- break_lineage(reference$lineage)

    ## Calculate alignment positions for hits
    hits <- hits %>%
        select(sacc, qaccs, qstart, qstop, sstart, sstop, cstart, cstop, flip, contig, pidents) %>%
        mutate(
            aleft  = cstart + qstart - 1,
            aright = cstart + qstop - 1
        )


    ## Get unique contigs from hits
    contigs <- hits %>%
        group_by(qaccs, flip, cstart, cstop, contig) %>%
        summarize(.groups="drop")

    ## Make compressed X scale
    xtrans <-
        bind_rows(
            hits %>% mutate(start=sstart, stop=sstop) %>% select(start, stop),
            hits %>% mutate(start=cstart, stop=cstop) %>% select(start, stop)
        ) %>%
        arrange(start, stop) %>%
        compress_axis()

    ## Setup corners for trapezoid showing alignment mapping
    alignment_boxes <- bind_rows(
        hits %>% mutate(y = ymax$alignments, x = aleft),
        hits %>% mutate(y = ymax$alignments, x = aright),
        hits %>% mutate(y = ymin$alignments, x = sstop),
        hits %>% mutate(y = ymin$alignments, x = sstart)
    ) %>%
        mutate(
            alignment = paste(qaccs, aleft, aright, sstop, sstart)
        )

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
            aes(xmin = cstart, xmax = cstop),
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
            aes(x = (cstart+cstop)/2, label = contig),
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
            aes(xmin=sstart, xmax=sstop),
            ymin=ymin$subject, ymax=ymax$subject, fill=fill$subject
        ) +
        ## Plot Alignment connection
        geom_polygon(
            data = alignment_boxes,
            aes(x = x, y = y, group = alignment, fill = pidents),
            color=fill$alignments, alpha=.5
        ) +
        scale_fill_continuous(limits = c(70, 100))


    if(!is.null(depths)) {
        label_tpl <- paste0(label_tpl, ", depth={mean_depth}")

        depths <- depths %>%
            rename(qaccs=contig) %>%
            inner_join(contigs, by="qaccs")

        mean_depth <- round(mean(depths$total), 1)
        max_depth <- max(depths$total)
        norm_depth <- pmin(mean_depth*5, max_depth)

        depth_data <- depths %>%
            mutate(
                x = if_else(flip, cstop - pos, cstart + pos),
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
                qaccs, sense, x, y
            )

        ## Plot depths above contigs
        p <- p +
            geom_line(
                data = depth_data,
                aes(
                    x = x,
                    y = y * heights$coverage + ymin$coverage,
                    color=sense,
                    group=interaction(sense,qaccs)
                )
            )
    }

    if (!is.null(feature_tables)) {
        ## Annotate subject sequences
        subject_annotations <- annotate_subjects(
            hits$sstart, hits$sstop,
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

    ## Add labels
    p <- p +
        xlab(str_glue(label_tpl)) +
        ylab(str_glue(ylabel_tpl))

    p
}


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
    message("Loading data files ...")
    data <- read_csv(opt$options$input, col_types = cols()) %>%
        filter(
            slen >= opt$options$min_slen,
            numreads >= opt$options$min_reads
        )

    message("Parsing input data ...")
    alignments <- parse_blastbins(data)
    message("Placing contigs ...")
    contigs <- place_contigs(alignments)
    ranges <- merge(contigs, alignments)

    message("Ordering results ...")
    saccs <- ranges %>%
        group_by(species) %>%
        mutate(
            sumreads = sum(numreads)
        ) %>%
        ungroup() %>%
        group_by(species, sumreads, sacc, numreads) %>%
        summarize(
            log_evalue=min(log_evalue),
            lineage=paste(unique(lineage), collapse=" ::: "),
            .groups="drop"
        ) %>%
        arrange(desc(sumreads), desc(numreads)) %>%
        filter(!is.na(sacc)) %>%
        pull(sacc)

    if (length(saccs) ==  0) {
        return(NULL)
    }

    if (opt$options$max_plots > 0) {
        saccs <- head(saccs, opt$options$max_plots)
    }

    if (!is.null(opt$options$input_bam)) {
        depths <- coverage_depth(opt$options$input_bam)
    } else {
        depths <- NULL

    }

    if (opt$options$do_annotate) {
        message("Loading ", length(saccs), " Feature Tables from Entrez...")
        feature_tables <- load_feature_table(saccs, opt$options$cache_path)
    } else {
        feature_tables <- NULL
    }

    message("Writing output to ", opt$options$output, " ...")
    pdf(
        file = opt$options$output,
        width = opt$options$page_width,
        height = opt$options$page_height
    )

    plot_pages(opt$options$plots_per_page, saccs, function(acc) {
        reference <- data %>% filter(sacc == acc)
        df <- ranges %>% filter(sacc == acc)
        plot_ranges(reference, df, depths, feature_tables)
    })

    dev.off()
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
    run()
}

