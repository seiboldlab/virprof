requireNamespace("IRanges", quietly = TRUE)
requireNamespace("rentrez", quietly = TRUE)
requireNamespace("rlang", quietly = TRUE)


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



#' Conditionally swap contents of columns
#'
#' For use with \code{mutate}
#'
#' @param cond If true, column contents will be swapped
#' @param x first column
#' @param y second column
#' @return [Tibble] Tibble with first and second column swapped where ``cond`` is TRUE.
#' @examples
#' data.frame(a=c(1,4,8), b=c(2,3,4)) %>%
#'   mutate(
#'     swapped = a<b,
#'     swap_if(swapped, a, b)
#'   )
#' @export
swap_if <- function(cond, x, y) {
    out_x <- if_else(cond, y, x)
    out_y <- if_else(!cond, y, x)
    setNames(
        tibble::tibble(out_x, out_y),
        c(substitute(x), substitute(y))
    )
}


#' Generates a "compressed" axis
#'
#' Takes a data frame of start/stop positions to include in the axis and translates
#' coordinates to omit ranges not between start/stop pairs.
#'
#' @param include_ranges Object with two columns containing start and stop positions
#' @param merge_dist Size of gap between ranges that should be merged
#' @param spacing Spacing between output rangess
#' @export
compress_axis <- function(include_ranges, merge_dist = 50, spacing = 50) {
    if (nrow(include_ranges) == 0) {
        stop("Cannot commpress axis with empty range")
    }
    if (spacing < 2) {
        stop("Spacing must be at least 2")
    }

    ranges <- IRanges::IRanges(include_ranges[[1]], include_ranges[[2]]) %>%
        IRanges::reduce(min.gapwidth = merge_dist) %>%
        as.data.frame() %>%
        mutate(
            toend = cumsum(end - start + spacing + 1) - spacing,
            tostart = toend - end + start
        ) %>%
        select(start, end, tostart, toend)

    convert <- function(x, from_start, from_end, to_start, to_end) {
        result <- vector("integer", length(x))
        result[is.na(x)] <- NA
        ## between start and end
        scale <- (to_start - to_end) / (from_start - from_end)
        for (j in seq_along(from_start)) {
            matches <- between(x, from_start[[j]], from_end[[j]])
            matches[is.na(matches)] <- FALSE
            result[matches] <- (x[matches] - from_start[[j]]) * scale[[j]] + to_start[[j]]
        }
        ## between end and start
        scale <- (
            as.double(to_start[-1] - to_end[-length(to_end)])
            /
            as.double(from_start[-1] - from_end[-length(to_end)])
        )
        for (j in seq_along(from_start[-1])) {
            matches <- between(x, from_end[[j]], from_start[[j+1]])
            matches[is.na(matches)] <- FALSE
            result[matches] <- (x[matches] - from_end[[j]]) * scale[[j]] + to_end[[j]]
        }
        ## outside first start and last end
        result[x<from_start[[1]] & !is.na(x)] <- x[x<from_start[[1]] & !is.na(x)] -
            from_start[[1]] + to_start[[1]]
        result[x>from_end[[length(from_end)]] & !is.na(x)] <- x[x>from_end[[length(from_end)]] & !is.na(x)] -
            from_end[[length(from_end)]] + to_end[[length(to_end)]]
        result
    }

    ## Translate regular coordinates to compressed coordinates
    trans <- function(x) {
        convert(x, ranges$start, ranges$end, ranges$tostart, ranges$toend)
    }

    ## Reverse translation
    inv <- function(x) {
        convert(x, ranges$tostart, ranges$toend, ranges$start, ranges$end)
    }

    ## Compute breaks (start/stop positions)
    breaks <- function(x) {
        unique(sort(c(ranges$start, ranges$end)))
    }

    list(
        "trans" = scales::trans_new("compress_axis", trans, inv, breaks),
        "ranges" = ranges
    )
}


#' Shifts overlapping intervals outwards until they have defined spacing
#'
#' @param contigs data.frame with cstart and cstop columns
#' @param spacing minimum distance between end and start of successive intervals
#' @param max_iter abort if no success after this many iterations
center_spread_ranges <- function(contigs, spacing, max_iter = 100) {
    ## Sort the intervals
    contigs <- arrange(contigs, cstart, cend)
    ## Make IRanges object from intervals with half spacing on either side
    ranges <- IRanges::IRanges(contigs$cstart, contigs$cend) + spacing / 2

    ## Do the following until nothing overlaps anymore
    iter <- 0
    while (!IRanges::isDisjoint(ranges) && iter < max_iter) {
        iter <- iter + 1
        ## Merge overlapping ranges
        mranges <- IRanges::reduce(ranges, with.revmap = TRUE)

        ## Extract list of ranges composing each overlapping range
        bins <- S4Vectors::mcols(mranges)$revmap
        ## Get number of ranges in each merged range
        bin_sizes <- lengths(bins)
        ## Calculate left shift required for each section of overlapping ranages
        ## Sum of widths of ranges for each set of overlapping:
        bin_widths <- sapply(bins, function(x) sum(IRanges::width(ranges)[x]))
        ## Minus actual width is what we need to shift
        bin_left_shift <- (bin_widths - IRanges::width(mranges) + 1) / 2
        ## Determine range start offset in bin
        range_start_offset <- unlist(sapply(bins, function(x) cumsum(IRanges::width(ranges[x])))) - IRanges::width(ranges)
        range_newstart <- range_start_offset + rep(IRanges::start(mranges), bin_sizes)
        range_right_shift <- range_newstart - IRanges::start(ranges)
        ranges <- IRanges::shift(ranges, range_right_shift - rep(bin_left_shift, bin_sizes))
    }
    if (iter >= max_iter) {
        print(paste("WARNING: aborted spacing out ranges after", iter, "iterations"))
    }
    ranges <- ranges - spacing / 2
    contigs$cstart <- IRanges::start(ranges)
    contigs$cend <- IRanges::end(ranges)
    contigs
}



#' Caching wrapper around \code{rentrez::entrez_fetch}
#'
#' @path Path to cache directory. Set to \code{""} to disable cache use.
cached_entrez_fetch <- function(db, id, rettype, path, retmode = "") {
    if (!nzchar(path)) {
        return(
            rentrez_entrez_fetch(
                db = db, id = id,
                rettype = rettype, retmode = retmode
            )
        )
    }
    dir.create(path, showWarnings = FALSE)
    result <- vector("list", length(id))
    for (i in seq_along(id)) {
        fname <- paste0(paste(sep="__", db, id[[i]], rettype, retmode), ".dat")
        if (file.exists(file.path(path, fname))) {
            result[[i]] <- read_file(file.path(path, fname))
        } else {
            result[[i]] <- rentrez::entrez_fetch(db=db, id=id[[i]], rettype=rettype, retmode=retmode)
            tmpname <- tempfile(fname, path, '.tmp')
            write_file(result[[i]], tmpname)
            file.rename(tmpname, file.path(path, fname))
        }
    }
    paste(collapse="\n", result)
}

#' Fetches and parses feature table for ``accs`` from Entrez
#'
#' @param accs List of accession numbers
#' @return Data frame with columns start, end, type, key and value
load_feature_table <- function(accs, cache_path) {
    feature_tables <-
        cached_entrez_fetch(db = "nucleotide", id=accs, rettype = 'ft', path = cache_path) %>%
        gsub("\n$", "", .) %>%
        str_split("(^|\n)>") %>%
        unlist() %>%
        str_split("\n", n = 2)
    feature_tables <- feature_tables[-1]
    tables <- vector("list", length(feature_tables))
    for (i in seq_along(feature_tables)) {
        acc <- sub("Feature gb\\|([^.]*)\\..*", "\\1", feature_tables[[i]][1])
        text <- feature_tables[[i]][2] %>%
            gsub("(^|[\n\t])[<>]", "\\1", .)  # remove <1, >123 from incomplete
        df <- read.csv(
            textConnection(text),
            col.names = c('start', 'end', 'type', 'key', 'value'),
            colClasses = c('numeric', 'numeric', 'character', 'character', 'character'),
            header = FALSE,
            sep = "\t"
        )
        if (nrow(df) > 0) {
            df$acc = acc
            tables[[i]] <- df
        }
    }
    df <- do.call("bind_rows", tables) %>% as_tibble()
    if (nrow(df) > 0) {
        df$type[df$type == ""] <- NA
        df <- df %>%
            fill(c(start, end, type), .direction = 'down') %>%
            filter(key!="")
    }
    df
}

#' Create annotated ranges for input ranges using feature table
#'
#' @param starts Start positions of ranges to annotate
#' @param ends End positions of ranges to annotate
#' @param acc Accession number of reference sequence
#' @param feature_table Feature table data
#' @param key Which feature key to use
#' @return data.frame with columns start, stop and name indicating the
#'     annotated regions. Columns gstart and gstop indicated the full
#'     range of the annotations.
#' @seealso load_feature_table
annotate_subjects <- function(starts, ends, acc_, feature_table) {
    ## Convert start/end lists into iranges object
    ranges <- IRanges::IRanges(starts, ends) %>%
        IRanges::reduce()

    ft <- feature_table %>%
        rowid_to_column() %>%
        filter(acc == acc_) %>%
        mutate(swap_if(start > end, start, end)) %>%
        {IRanges::IRanges(.$start, .$end, name = .$rowid)} %>%
        IRanges::subsetByOverlaps(ranges)

    overlaps <- IRanges::intersect(ranges, ft)
    hits <- IRanges::findOverlaps(overlaps, ft)
    overlaps_out <- overlaps[S4Vectors::from(hits)]
    ft_out <- ft[S4Vectors::to(hits)]
    index_out <- as.integer(names(ft_out))
    res <- data.frame(
        start = IRanges::start(overlaps_out),
        stop  = IRanges::end(overlaps_out),
        gstart = IRanges::start(ft_out),
        gstop = IRanges::end(ft_out),
        key = feature_table$key[index_out],
        value = feature_table$value[index_out]
    )
    res
}

#' Execute bedtools genomecov
#'
#' @param strand +/- to get stranded coverage (-du -strand +, du makes
#'     it so mates are counted towards same strand).
#' @param split If true, don't count gaps (-split)
#' @param fragment Count coverage for full fragment length (-fs)
run_bedtools <- function(fname, strand=NULL, split=FALSE, fragment=FALSE) {
    message(
        "Running bedtools genomecov on ", fname,
        if (!is.null(strand)) paste(" counting", strand, "strand"),
        if (split) " not counting gaps between spliced alignments",
        if (fragment) " counting whole fragment",
        "..."
    )
    proc <- pipe(paste(
        "bedtools", "genomecov",
        paste("-ibam", strsplit(fname, ",")[[1]]),
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




#' VP data class

setClass(
    "VirProf",
    slots = list(
        calls = "data.frame",
        alignments = "data.frame",
        depths = "data.frame",
        scaffold_depths = "data.frame",
        features = "data.frame"
    )
)


VirProfFromCSV <- function
(
    calls_fname,  # base filename  (.virus.csv)
    hits_fname = NULL,   # individual hits (.hits.csv),
    features_fname = NULL,  # feature table
    min_slen = 0,
    min_reads = 0
) {
    res <- new("VirProf")
    stopifnot(file.exists(calls_fname))
    if (is.null(hits_fname)) {
        hits_fname <- gsub("virus\\.csv$", "hits.csv", calls_fname)
    }
    stopifnot(file.exists(hits_fname))
    if (is.null(features_fname)) {
        fname <- gsub("virus\\.csv$", "features.csv", calls_fname)
        if (file.exists(fname)) {
            features_fname <- fname
        }
    } else {
        stopifnot(file.exists(features_fname))
    }
    message("Loading calls...")
    calls <- read_csv(calls_fname, col_types = cols())
    message("... ", nrow(calls), " calls found")
    calls <- filter(calls, slen >= opt$options$min_slen)
    message("... ", nrow(calls), " calls left after removing calls with <",
            opt$options$min_slen, " bp on subject")
    calls <- filter(calls, numreads >= opt$options$min_reads)
    message("... ", nrow(calls), " calls left after removing calls with <",
            opt$options$min_reads, " mapped reads")
    if (nrow(calls) == 0) {
        message("No calls left, exiting")
        return(res)
    }
    calls <- calls %>%
        mutate(
            genome_coverage = if_else(genome_size > 0,
                                      round(slen / genome_size * 100, 1),
                                      NA_real_)
        )
    res@calls <- calls
    message("Loading alignments...")
    alignments <- read_csv(hits_fname, col_types = cols())
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
    res@alignments <- alignments
    if (is.null(features_fname)) {
        message("Not loading feature tables")
    } else {
        message("Loading feature tables...")
        feature_tables <- read_csv(opt$options$input_features, col_types = cols())
        message("... ", nrow(feature_tables), " feature annotations found")
        res@features <- feature_tables
    }
    return(res)
}


setMethod("show", "VirProf", function(object) {
    cat("# VirProf Results\n")
    cat("# ", nrow(object@calls), " calls\n")
    cat("# ", nrow(object@alignments), " alignments\n")
    cat("# ", nrow(object@features), " features\n")
})




#' Create offsets for qacc positions in plot
#'
#' @param alignments data.frame with columns sacc, qacc, sstart, send, qstart, qend, qlen
#' @return alignments with added columns cstart, cend
setGeneric("place_contigs", function(object, ...) {
    standardGeneric("place_contigs")
})
setMethod(
    "place_contigs",
    signature("VirProf"),
    function(object, spacing = 10) {
        alignments <- object@alignments
        if (nrow(alignments) == 0) {
            empty_res <- data.frame(
                sacc = character(),
                qacc = character(),
                cstart = numeric(),
                cend = numeric()
            )
            return(empty_res)
        }
        ## Align each contig to have leftmost qstart be above sstart
        contigs <- alignments %>%
            group_by(sacc, qacc) %>%
            slice_min(qstart, n = 1, with_ties = FALSE) %>%
            mutate(
                cstart = (sstart+send)/2 - (qstart+qend)/2 + 1,
                cend = cstart + qlen
            )  %>%
            group_by(sacc) %>%
            group_modify(~ center_spread_ranges(., spacing)) %>%
            ungroup() %>%
            select(sacc, qacc, cstart, cend)
        object@alignments <- alignments %>% left_join(contigs, by=c("sacc", "qacc"))
        return(object)
    })


#' Get depth from BAM file
#'
#' @param fname Path to sorted BAM file
setGeneric("coverage_depth", function(object, ...) standardGeneric("coverage_depth"))
setMethod("coverage_depth", signature("VirProf"), function(object, fname, scaffold=FALSE) {
    plus <- run_bedtools(fname, strand="+", split=TRUE) %>% rename(plus=depth)
    minus <- run_bedtools(fname, strand="-", split=TRUE) %>% rename(minus=depth)
    result <- full_join(plus, minus, by=c("contig", "pos"))
    result$total = result$plus + result$minus
    if (scaffold) {
        object@scaffold_depths <- result
    } else {
        object@depths <- result
    }
    return(object)
})


##plot_ranges <- function(reference, hits, depths, featureshowMethods("plot")
#' The actual plot function
setMethod("plot", signature("VirProf"), function
(
    x,
    accession,
    box_colors = c(
        "Contig" = "darkblue",
        "Reference" = "darkred",
        "gene" = "orange",
        "product" = "darkorange"
    ),
    label_colors = c(
        "Contig" = "white",
        "Reference" = "white",
        "gene" = "black",
        "product" = "black"
    ),
    heights = c(
        scaffold_coverage=5,
        coverage = 5,
        contigs = 1,
        alignments = 3,
        subject = 1,
        annotations = 1
    ),
    label_tpl = paste0(
        "Keywords: {reference$words}\n",
        "{reference$sacc}: \"{reference$stitle}\"\n",
        "{lineage}\n",
        "Read Count: {reference$numreads}, ",
        "log E-Value: {reference$log_evalue}, ",
        "Identity: {reference$pident}%, ",
        "Genome Coverage {reference$slen}bp ({reference$genome_coverage}%)"
    ),
    ylabel_tpl = paste0(
        "{reference$species}\n",
        "{reference$sacc}"
    )
)
{
    message("Plotting ", accession)
    reference <- filter(x@calls, sacc == accession)
    hits <- filter(x@alignments, sacc == accession)
    depths <- x@depths
    feature_tables <- x@features

    #### CONFIG  ###
    cat(".")
    ymax <- cumsum(rev(heights))
    ymin <- ymax - rev(heights)
    ymax <- ymax - ymin[["coverage"]]
    ymin <- ymin - ymin[["coverage"]]
    heights <- as.list(heights)
    ymax <- as.list(ymax)
    ymin <- as.list(ymin)

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
                y = total / cap_depth * heights$coverage + ymin$coverage,
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
    if (!is.null(x@scaffold_depths)) {
        message("Showing scaffold depths")
        depths <- x@scaffold_depths %>%
            filter(sacc == accession)
#        %>%
#            mutate(contig=gsub("_pilon", "", contig)) %>%
#            separate(contig, into=c("sample", "sacc"), sep=".")
#            rename(qacc=contig) %>%
#            inner_join(contigs, by="qacc")

        sd_depth <- sd(depths$total)
        max_depth <- max(depths$total)
        cap_depth <- min(max_depth, mean_depth + 3 * sd_depth)

        depth_data <- depths %>%
            mutate(
                y = total / cap_depth * heights$scaffold_coverage + ymin$scaffold_coverage,
                x = pos
            ) %>%
            select(
                sacc, x, y
            )
        message("1")

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
        message("2")
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
})
