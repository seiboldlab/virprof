requireNamespace("IRanges", quietly = TRUE)


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
    ## Do nothing if the range list is empty
    if (nrow(include_ranges) == 0) {
        return("identity")
    }

    ## Create translation ranges
    ranges <- IRanges::IRanges(include_ranges[[1]], include_ranges[[2]]) %>%
        IRanges::reduce(min.gapwidth = merge_dist) %>%
        as.data.frame() %>%
        mutate(
            stop = end,
            newstop = cumsum(stop - start + spacing) - spacing + 1,
            newstart = newstop - stop + start
        ) %>%
        select(start, stop, newstart, newstop) %>%
        arrange(newstart, newstop)

    ## Extend outside ranges
    extra <- 10000000
    ranges$start[[1]] <- ranges$start[[1]] - extra
    ranges$newstart[[1]] <- ranges$newstart[[1]] - extra
    ranges$stop[[nrow(ranges)]] <- ranges$stop[[nrow(ranges)]] + extra
    ranges$newstop[[nrow(ranges)]] <- ranges$newstop[[nrow(ranges)]] + extra

    ## Translate regular coordinates to compressed coordinates
    trans <- function(x) {
        x_na = is.na(x)
        seen = is.na(x)
        for (j in 1:nrow(ranges)) {
            matches <- between(x, ranges$start[[j]], ranges$stop[[j]] + 1)
            matches <- matches & !seen
            seen <- matches | seen
            x[matches] <- x[matches] - ranges$stop[[j]] + ranges$newstop[[j]]
        }
        x[x_na] <- NA
        if (any(!seen)) {
            print(ranges)
            print(x[!seen])
            stop("Values out of display range!")
        }
        x
    }

    ## Reverse translation
    inv <- function(x) {
        x_na = is.na(x)
        seen = is.na(x)
        for (j in 1:nrow(ranges)) {
            matches <- between(x, ranges$newstart[[j]], ranges$newstop[[j]])
            matches <- matches & !seen
            seen <- matches | seen
            x[matches] <- x[matches] - ranges$newstop[[j]] + ranges$stop[[j]]
        }

        if (any(!seen)) {
            print(ranges)
            print(x)
            print(seen)
            warning("Values out of display range!")
        }
        x[x_na] <- NA
        x
    }

    ## Compute breaks (start/stop positions)
    breaks <- function(x) {
        x <- unique(sort(c(ranges$start, ranges$stop)))
        x
    }

    scales::trans_new("compress_axis", trans, inv, breaks)
}

#' Shifts overlapping intervals outwards until they have defined spacing
#'
#' @param contigs data.frame with cstart and cstop columns
#' @param spacing minimum distance between end and start of successive intervals
center_spread_ranges <- function(contigs, spacing) {
    ## Sort the intervals
    contigs <- arrange(contigs, cstart, cstop)
    ## Make IRanges object from intervals with half spacing on either side
    ranges <- IRanges::IRanges(contigs$cstart, contigs$cstop) + spacing / 2

    ## Do the following until nothing overlaps anymore
    max_iter <- 10
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
    contigs$cstop <- IRanges::end(ranges)
    contigs
}


#' Create offsets for qacc positions in plot
#'
#' @param alignments data.frame with columns sacc, qaccs, sstart, sstop, qstart, qstop, qlen
#' @return data.frame with columns sacc, qaccs, cstart, cstop
place_contigs <- function(alignments, merge_dist = 50, spacing = 10) {
    if (nrow(alignments) == 0) {
        emptyres <- data.frame(
            sacc = character(),
            qaccs = character(),
            cstart = numeric(),
            cstop = numeric()
            )
        return(emptyres)
    }

    ## Align each contig to have leftmost qstart be above sstart
    contigs <- alignments %>%
        group_by(sacc, qaccs) %>%
        slice_min(qstart, n = 1, with_ties = FALSE) %>%
        mutate(
            ##cstart = sstart - qstart + 1,
            cstart = (sstart+sstop)/2 - (qstart+qstop)/2 + 1,
            cstop = cstart + qlen
        )  %>%
        group_by(sacc) %>%
        group_modify(~ center_spread_ranges(., spacing))

    contigs
}
