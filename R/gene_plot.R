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
        select(start, stop, newstart, newstop)

    ## Translate regular coordinates to compressed coordinates
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
