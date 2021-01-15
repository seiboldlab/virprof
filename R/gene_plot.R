requireNamespace("IRanges", quietly = TRUE)
requireNamespace("rentrez", quietly = TRUE)
requireNamespace("rlang", quietly = TRUE)

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
#' @param max_iter abort if no success after this many iterations
center_spread_ranges <- function(contigs, spacing, max_iter = 100) {
    ## Sort the intervals
    contigs <- arrange(contigs, cstart, cstop)
    ## Make IRanges object from intervals with half spacing on either side
    ranges <- IRanges::IRanges(contigs$cstart, contigs$cstop) + spacing / 2

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
    contigs$cstop <- IRanges::end(ranges)
    contigs
}


#' Create offsets for qacc positions in plot
#'
#' @param alignments data.frame with columns sacc, qaccs, sstart, sstop, qstart, qstop, qlens
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
            cstart = (sstart+sstop)/2 - (qstart+qstop)/2 + 1,
            cstop = cstart + qlens
        )  %>%
        group_by(sacc) %>%
        group_modify(~ center_spread_ranges(., spacing)) %>%
        ungroup()

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
annotate_subjects <- function(starts, ends, acc, feature_table, key = "gene") {
    acc_ <- acc
    key_ <- key
    subj <- IRanges::IRanges(starts, ends) %>%
        IRanges::reduce()
    ft <- feature_table %>%
        filter(acc == acc_, key == key_) %>%
        mutate(swap_if(start > end, start, end)) %>%
        {IRanges::IRanges(.$start, .$end, name = .$value)} %>%
        IRanges::subsetByOverlaps(subj)
    overlaps <- IRanges::intersect(subj, ft)
    hits <- IRanges::findOverlaps(overlaps, ft)
    overlaps_out <- overlaps[S4Vectors::from(hits)]
    ft_out <- ft[S4Vectors::to(hits)]
    data.frame(
        start = IRanges::start(overlaps_out),
        stop  = IRanges::end(overlaps_out),
        gstart = IRanges::start(ft_out),
        gstop = IRanges::end(ft_out),
        name = names(ft_out)
    )
}
