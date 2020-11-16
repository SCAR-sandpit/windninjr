#' Demonstration files bundled with the package
#'
#' @param what string: which file?
#' \itemize{
#'  \item "missoula_valley_elevation" - the Missoula Valley elevation geotiff file bundled with WindNinja itself
#' }
#'
#' @return The path to the file
#'
#' @export
wn_demo_file <- function(what = "missoula_valley_elevation") {
    system.file("extdata/missoula_valley.tif", package = "windninjr")
}
