#' Convenience function for reading WindNinja output files
#'
#' @param path string: path to output directory
#' @param pattern string: an optional regular expression (e.g. \code{"_270_10_100m_"}) passed on to grep() to filter files
#' @param vars character: a character list of variables to return 
#' @param format string: output format to look for
#'
#' @return A raster stack
#'
#' @export
wn_read <- function(path, pattern, vars = c("vel", "ang"), format = "asc") {
    if (missing(pattern) || is.na(pattern) || !nzchar(pattern)) pattern <- NULL
    fls <- fs::dir_ls(path, regexp = pattern)
    if (!is.null(format) && nzchar(format)) fls <- fls[grepl(paste0(format, "$"), fls)]
    if (!is.null(vars)) fls <- fls[grepl(paste0("(", paste0(vars, collapse = "|"), ")\\."), fls)]
    out <- raster::stack(fls)
    names(out) <- basename(fls)
    out
}
