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


#' Create offsets for qacc positions in plot
#'
#' @param alignments data.frame with columns sacc, qacc, sstart, send, qstart, qend, qlen
#' @return alignments with added columns cstart, cend
place_contigs <- function(alignments, spacing = 10) {
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

    alignments %>% left_join(contigs, by=c("sacc", "qacc"))
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


#' Get depth from BAM file
#'
#' @param fname Path to sorted BAM file
coverage_depth <- function(fname) {
    plus <- run_bedtools(fname, strand="+", split=TRUE) %>% rename(plus=depth)
    minus <- run_bedtools(fname, strand="-", split=TRUE) %>% rename(minus=depth)
    result <- full_join(plus, minus, by=c("contig", "pos"))
    result$total = result$plus + result$minus
    result
}


}

