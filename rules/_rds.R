#' Faster saveRDS/readRDS using external process for compression
#'
#' Uses `pigz` to handle compression outside of R main process. Will
#' be many times faster saving than regular saveRDS without much
#' increase in size. Reading regular gzip'ed RDS is only slightly
#' faster as deserialization takes most the CPU time.
#'
#' If `pv` is available, it will be used to show progress. For
#' saveRDS, it shows compressed bytes written, time elapsed and
#' current throughput. For readRDS, it shows compressed bytes read,
#' time elapased and a progress bar.
#'
#' @param object The object to store
#' @param file The file name (can't be a connection)
#' @param level 0, ..., 11, overrides compression level
#' @param threads Override auto-detected number of threads
#'
#' @export
saveRDS <- function(object, file = "", level = NULL, threads = NULL) {
    args <- list()
    args$pigz <- "pigz --stdout"
    if (!is.null(level)) {
        if (length(level) != 1) {
            rlang::abort("level must be length 1")
        }
        if (!is.numeric(level)) {
            rlang::abort("level must be numeric")
        }
        if (level < 0 || level > 11) {
            rlang::abort("level must be between 0 and 11 (inclusive)")
        }
        args$level <- paste0("-", as.integer(level))
    }
    if (!is.null(threads)) {
        if (length(threads) != 1) {
            rlang::abort("threads must be length 1")
        }
        if (!is.numeric(threads)) {
            rlang::abort("threads must be numeric")
        }
        if (!is.numeric(threads) || threads < 0) {
            rlang::abort("threads must be above 0")
        }
        args$threads <- paste("-p", threads)
    }
    if (interactive() && nzchar(Sys.which("pv"))) {
        args$pv <- "| pv --rate --timer --bytes "
    }
    if (length(file) != 1) {
        rlang::abort("file must be length 1")
    }
    if (!is.character(file)) {
        str(file)
        rlang::abort("file must be character")
    }
    if (!nzchar(file)) {
        rlang::abort("file must be non-empty")
    }
    args$pipe <- ">"
    args$file <- file
    if (!nzchar(Sys.which("pigz"))) {
        warning(
            "fast saveRDS: cannot find pigz in PATH, falling back to base::saveRDS"
        )
        return(base::saveRDS(object, file))
    }
    cmd <- paste(args, collapse = " ")
    con <- pipe(cmd)
    on.exit(close(con))
    base::saveRDS(object, con)
}

#' @rdname saveRDS
readRDS <- function(file = "", threads = NULL) {
    args <- list()
    if (length(file) != 1) {
        rlang::abort("file must be length 1")
    }
    if (!is.character(file)) {
        rlang::abort("file must be character")
    }
    if (!nzchar(file)) {
        rlang::abort("file must be non-empty")
    }
    if (interactive() && nzchar(Sys.which("pv"))) {
        args$pv <- "pv --progress --timer --bytes"
        args$file <- file
        args$pigz <- "| pigz --decompress --stdout"
    } else {
        args$pigz <- "pigz --decompress --stdout"
        args$file <- file
    }

    if (!is.null(threads)) {
        if (length(threads) != 1) {
            rlang::abort("threads must be length 1")
        }
        if (!is.numeric(threads)) {
            rlang::abort("threads must be numeric")
        }
        if (!is.numeric(threads) || threads < 0) {
            rlang::abort("threads must be above 0")
        }
        args$threads <- paste("-p", threads)
    }

    if (!nzchar(Sys.which("pigz"))) {
        warning(
            "fast readRDS: cannot find pigz in PATH, falling back to base::readRDS"
        )
        return(base::readRDS(object, file))
    }

    cmd <- paste(args, collapse = " ")
    con <- pipe(cmd)
    on.exit(close(con))
    force(base::readRDS(con))
}
