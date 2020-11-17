#!/usr/bin/env Rscript
library(optparse)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(grid)


parse_options<- function(args = commandArgs(trailingOnly = TRUE)) {
    option_list <- list(
        make_option(c("--input"),
                    metavar = "FILE",
                    help = "Result file from blastbin"
                    ),
        make_option(c("--input-bam"),
                    metavar = "FILE",
                    help = "Sorted BAM mapping reads to contigs"
                    ),
        make_option(c("--output"),
                    metavar = "FILE",
                    help = "Output PDF (default=%default)",
                    default = "out.pdf"
                    ),
        make_option(c("--page-width"),
                    metavar = "INCHES",
                    help = "Output PDF width (default=%default)",
                    default=22
                    ),
        make_option(c("--page-height"),
                    metavar = "INCHES",
                    help = "Ootput PDF height (default=%default)",
                    default=18
                    ),
        make_option(c("--plots-per-page"),
                    metavar = "N",
                    help = "Number of plots per page",
                    default=5
                    ),
        make_option(c("--max-plots"),
                    metavar = "N",
                    help = "Maximum number of plots to output",
                    default=0),
        make_option(c("--min-alen"),
                    metavar = "LEN",
                    help = "Minimum alen for rendered hits",
                    default = 200),
        make_option(c("--min-reads"),
                    metavar = "N",
                    help = "Minuimum read count for rendered hits",
                    default = 3)
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


compress_axis <- function(occ_ranges) {
    if (nrow(occ_ranges) == 0) {
        return("identity")
    }
    cur_start <- cur_stop <- 0
    for (j in 1:nrow(occ_ranges)) {
        rng <- occ_ranges[j,]
        if (rng$start > cur_stop + 50) {
            cur_start <- rng$start
        } else {
            rng$start <- cur_start
        }
        cur_stop <- max(cur_stop, rng$stop)
        rng$stop <- cur_stop
        occ_ranges[j,] <- rng
    }
    ranges <- occ_ranges %>%
        group_by(start) %>%
        summarize(stop=max(stop), .groups="drop") %>%
        ungroup() %>%
        mutate(
            newstop = cumsum(stop-start+50),
            newstart = newstop-(stop-start)
        )

    inv <- function(x) {
        x_na = is.na(x)
        seen = is.na(x)
        for (j in 1:nrow(ranges)) {
            matches <- between(x, ranges$newstart[[j]], ranges$newstop[[j]])
            matches <- matches & !seen
            seen <- matches | seen
            x[matches] <- x[matches] - ranges$newstop[[j]] + ranges$stop[[j]]
        }
        x[x_na] <- NA
        x
    }

    trans <- function(x) {
        x_na = is.na(x)
        seen = is.na(x)
        for (j in 1:nrow(ranges)) {
            matches <- between(x, ranges$start[[j]], ranges$stop[[j]])
            matches <- matches & !seen
            seen <- matches | seen
            x[matches] <- x[matches] - ranges$stop[[j]] + ranges$newstop[[j]]
        }
        x[x_na] <- NA
        x
    }

    breaks <- function(x) {
        x <- unique(sort(c(ranges$start, ranges$stop)))
        x
    }

    scales::trans_new("compress_axis", trans, inv, breaks)
}


#' Turn each alignment into a row of its own
split_alignments <- function(data) {
    data %>%
        separate_rows(qranges, sranges, reversed, qaccs, pidents, numreadss, sep=";") %>%
        mutate(
            pidents=as.numeric(pidents),
            numreadss=as.numeric(numreadss),
            reversed=as.logical(reversed=="T")
        ) %>%
        separate(qranges, c("qstart", "qstop"), convert = TRUE) %>%
        separate(sranges, c("sstart", "sstop"), convert = TRUE) %>%
        arrange(-log_evalue, sacc, sstart)
}


#' Get depth from BAM file
coverage_depth_ <- function(fname, strand=NULL, split=FALSE, fragment=FALSE) {
    message("Running bedtools genomecov on ", fname,
            if (!is.null(strand)) paste(" counting", strand, "strand"),
            if (split) " not counting gaps",
            if (fragment) " counting whole fragment")
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
coverage_depth <- function(fname) {
    plus <- coverage_depth_(fname, strand="+", split=TRUE) %>% rename(plus=depth)
    minus <- coverage_depth_(fname, strand="-", split=TRUE) %>% rename(minus=depth)
    result <- full_join(plus, minus, by=c("contig", "pos"))
    result$total = result$plus + result$minus
    result
}


#' Create offsets for qacc positions in plot
place_contigs <- function(alignments) {
    if (nrow(alignments) == 0) {
        emptyres <- data.frame(
            sacc = character(),
            qaccs = character(),
            cstart = numeric(),
            cstop = numeric()
            )
        return(emptyres)
    }
    contigs <- alignments %>%
        group_by(sacc, qaccs) %>%
        summarize(
            offset = min(sstart),
            cstart = offset + min(qstart) -1,
            cstop = cstart + as.integer(sub("NODE_.*_length_([0-9]*)_.*", "\\1", unique(qaccs))),
            .groups = "drop"
        ) %>%
        select(-offset) %>%
        arrange(sacc, cstart)
    last_stop <- -1
    last_sacc <- ""
    for (i in 1:nrow(contigs)) {
        contig <- contigs[i,]
        if (contig$sacc != last_sacc) {
            last_sacc <- contig$sacc
            last_stop <- -10
        }
        if (contig$cstart <= last_stop + 10) {
            offset <- last_stop - contig$cstart + 10
            #message("Shifting ", contig$sacc, ":", contig$qaccs,
            #        " by ", offset)
            contig$cstart <- contig$cstart + offset
            contig$cstop <- contig$cstop + offset
            contigs[i,] <- contig
        }
        last_stop <- contig$cstop
    }
    contigs
}


break_lineage <- function(lineage, maxlen=200, insert="\n") {
    lineage %>%
        gsub(" ", "_", .) %>%
        gsub(";_", " ", .) %>%
        strwrap(maxlen) %>%
        gsub(" ", ";_", .) %>%
        gsub("_", " ", .) %>%
        paste(collapse=paste0(";", insert))
}


plot_ranges <- function(reference, hits, depths) {
    heights <- c(
        labels = 1.5,
        coverage = 4,
        contigs = 1,
        contig_hits = .5,
        alignments = 4,
        subject = 1
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
        subject     = "darkred"
    )

    label_tpl <- paste0(
        "keywords: {reference$words}\n",
        "{reference$sacc}: \"{reference$stitle}\"\n",
        "{lineage}\n",
        "nreads={reference$numreads}, ",
        "log(evalue)={reference$log_evalue}, ",
        "id={reference$pident}%, ",
        "qlen={reference$qlen}, ",
        "alen={reference$alen}"
    )

    ylabel_tpl <- paste0(
        "{reference$species}\n",
        "{reference$sacc}"
    )

    hits <- hits %>%
        ## decide on orientation of hits
        group_by(qaccs) %>%
        mutate(
            crev = as.logical(median(reversed)),
            cleft  = max(cstart),
            cright = max(cstop)
        ) %>%
        ungroup() %>%
        mutate(
            aleft  = if_else(crev, cright - qstart, cleft + qstart),
            aright = if_else(crev, cright - qstop,  cleft + qstop)
        ) %>%
        select(qaccs, contig, pidents,
               crev, cleft, cright, sstart, sstop, aleft, aright)

    contigs <- hits %>%
        group_by(qaccs, crev, cleft, cright, contig) %>%
        summarize(.groups="drop")

    if (!is.null(depths)) {
        label_dpl <- paste0(label_tpl, ", depth={mean_depth}")

        depths <- depths %>%
            rename(qaccs=contig) %>%
            inner_join(contigs, by="qaccs")

        mean_depth <- round(mean(depths$total), 1)
        max_depth <- max(depths$total)
        norm_depth <- pmin(mean_depth*5, max_depth)


        depth_data <- depths %>%
            mutate(
                x = if_else(crev, cright - pos, cleft + pos),
                fplus = if_else(total == 0, 0.5, if_else(crev, minus, plus) / total),
                plusx = if_else(crev, minus, plus) / norm_depth,
                minusx = if_else(crev, plus, minus) / norm_depth,
                plus = plusx,
                minus = minusx,
                total = total / norm_depth
            ) %>%
            gather(
                total, plus, minus, fplus, key="sense", value="y"
            ) %>%
            select(
                qaccs, sense, x, y
            )
    }

    p <- ggplot() +
        scale_y_continuous(limits = c(0, 16)) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle=90, hjust=1, size=6),
            axis.text.y = element_blank()
        )

    ## Contigs
    p <- p +
        geom_rect(
            data = contigs,
            aes(xmin = cleft, xmax = cright),
            ymin = ymin$contigs, ymax = ymax$contigs,
            fill = fill$contigs
        )

    ## Alignments on contigs
    p <- p +
        geom_rect(
            data = hits,
            aes(xmin = aleft, xmax = aright),
            ymin = ymin$contig_hits, ymax = ymax$contig_hits,
            fill = fill$contig_hits
        )

    ## Contig labels
    p <- p +
        geom_text_repel(
            data = contigs,
            aes(x = (cleft+cright)/2, label = contig),
            y = ymin$labels,
            size = 2,
            nudge_y = 3,
            direction  = "both",
            angle = 0,
            vjust = 0,
            segment.size = 0.5,
            )

    if(!is.null(depths)) {
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

    ## Make compressed X scale
    occ_ranges <-
        bind_rows(
            hits %>% mutate(start=sstart, stop=sstop) %>% select(start, stop),
            hits %>% mutate(start=cleft, stop=cright) %>% select(start, stop)
        ) %>%
        arrange(start, stop)

    p <- p +
        scale_x_continuous(trans = compress_axis(occ_ranges))

    ## Plot Subject sequence pieces
    p <- p +
        geom_rect(data = hits,
                  aes(xmin=sstart, xmax=sstop),
                  ymin=ymin$subject, ymax=ymax$subject, fill=fill$subject
                  )
    ## Plot Alignment connection
    alignment_boxes <- bind_rows(
        hits %>% mutate(y = ymax$alignments, x = if_else(crev, aright, aleft)),
        hits %>% mutate(y = ymax$alignments, x = if_else(crev, aleft, aright)),
        hits %>% mutate(y = ymin$alignments, x = sstop),
        hits %>% mutate(y = ymin$alignments, x = sstart)
    ) %>%
        mutate(
            alignment = paste(qaccs, aleft, aright, sstop, sstart)
        )
    p <- p +
        geom_polygon(
            data = alignment_boxes,
            aes(x = x, y = y, group = alignment, fill = pidents),
            color=fill$alignments, alpha=.5
        ) +
        scale_fill_continuous(limits = c(70, 100))

    ## Add labels
    lineage <- break_lineage(reference$lineage)
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
    data <- read_csv(opt$options$input, col_types = cols()) %>%
        filter(
            alen >= opt$options$min_alen,
            numreads >= opt$options$min_reads
        )

    pdf(
        file = opt$options$output,
        width = opt$options$page_width,
        height = opt$options$page_height
    )

    alignments <- split_alignments(data)
    contigs <- place_contigs(alignments)
    ranges <- merge(contigs, alignments)

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

    if (opt$options$max_plots > 0) {
        saccs <- head(saccs, opt$options$max_plots)
    }

    if (length(saccs) > 0) {
        if (!is.null(opt$options$input_bam)) {
            depths <- coverage_depth(opt$options$input_bam)
        } else {
            depths <- NULL
        }

        plot_pages(opt$options$plots_per_page, saccs, function(acc) {
            reference <- data %>% filter(sacc == acc)
            df <- ranges %>% filter(sacc == acc) %>%
                mutate(
                    contig=sub("NODE_(.*)_length_([0-9]*)_.*", "\\1", qaccs)
                )
            plot_ranges(reference, df, depths)
        })
   }

    dev.off()
}

if (!interactive()) {
    if (T) {
        opt <- parse_options()
    } else {
        opt <- parse_options(args=c("--input", "out.csv", "--input-bam", "out.bam", "--output", "out.pdf"))

        opt <- parse_options(args=c("--input", "out.csv", "--input-bam", "out.bam", "--output", "out.pdf",
                                    "--max-plots", "5"))
    }
    run()
}

#opt <- parse_options(args=c("--input", "inspire.qc.ref_hg38.deplete.assemble_meta.coverage.ref_hg38.annotate_blastE2MegaBest.blastfilter_vpU200.ref_NT.annotate_blastE10Best.ref_NcbiTaxonomy.blastbin_vp/C5042X4_120617_L007.virus.csv", "--output", "out.pdf"))
