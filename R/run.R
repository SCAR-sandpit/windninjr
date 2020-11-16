#' Run WindNinja
#'
#' @param config list: configuration object as returned by e.g. \code{\link{wn_config_domain_average}}
#' @param output_dir string: path to output directory
#' @param create_output_dir logical: if \code{TRUE} and the output directory does not exist, create it
#'
#' @return A list
#'
#' @examples
#' \dontrun{
#'   demfile <- wn_demo_file("missoula_valley_elevation")
#'   my_config <- wn_config_domain_average(elevation = demfile,
#'                                         input_speed = 10, input_direction = 270)
#'   res <- wn_run(my_config)
#' }
#'
#' @export
wn_run <- function(config, output_dir = tempfile(), create_output_dir = TRUE) {
    exe <- wn_find_exe()
    ## create the output dir if necessary
    if (!fs::dir_exists(output_dir) && isTRUE(create_output_dir)) fs::dir_create(output_dir, recurse = TRUE)
    if (!fs::dir_exists(output_dir)) stop("output directory ", output_dir, " does not exist, either create it or specify create_output_dir = TRUE")
    ## add output_path to config
    config$output_path <- output_dir
    ## write config to file
    config_file <- tempfile(tmpdir = output_dir, fileext = ".cfg")
    wn_write_config_file(config_file, config)
    ## output_path should make this unnecessary
    ##cwd <- getwd()
    ##on.exit(setwd(cwd))
    ##setwd(output_dir)
    out <- sys::exec_wait(exe, config_file)
    list(output_dir = output_dir, out = out)
}

# Find or set the path to the WindNinja command-line executable
#
# @param path string: optionally pass a path to look in
# @param error logical: only relevant if \code{path} is not specified (i.e. we are searching for the executable). If the executable is not found and \code{error} is \code{TRUE}, raise an error. If \code{FALSE}, do not raise an error but return NULL
#
# @return the path to the executable, or (if error is \code{FALSE}) NULL if it was not found
#
# @examples
# \dontrun{
#   wn_find_exe() ## search on the system path for the executable
#   wn_find_exe(path = "c:/my/custom/location") ## 
# }
#
# @export
wn_find_exe <- function(path, error = FALSE) {
    wn_opts <- getOption("windninjr")
    if (!missing(path) && !is.null(path) && nzchar(path)) {
        if (wn_exe_test(path)) {
            wn_opts$exe_path <- path
            options(windninjr = wn_opts)
            return(path)
        } else if (fs::is_dir(path)) {
            temp <- fs::dir_ls(path, regexp = "WindNinja_cli", recurse = TRUE, ignore.case = TRUE)
            if (length(temp) == 1) {
                return(wn_find_exe(path = temp, error = error))
            } else {
                stop("path '", path, "' does not appear to be the WindNinja command-line executable")
            }
        } else {
            stop("path '", path, "' does not appear to be the WindNinja command-line executable")
        }
    }
    if (!is.null(wn_opts)) {
        if (!is.null(wn_opts$exe_path)) {
            ## already set
            ## check that it actually exists
            if (file.exists(wn_opts$exe_path)) {
                ## ok
                return(wn_opts$exe_path)
            } else {
                ## it's disappeared somehow
                wn_opts$exe_path <- NULL
                options(windninjr = wn_opts)
                ## try again
                return(wn_find_exe(error = error))
            }
        }
    } else {
        wn_opts <- list()
    }
    if (wn_exe_test("WindNinja_cli")) {
        myexe <- unname(Sys.which("WindNinja_cli"))
    } else {
        my_os <- get_os()
        myexe <- NULL
        if (my_os == "windows") {
            dd <- dir("c:/Program Files/WindNinja", recursive = TRUE, pattern = "WindNinja_cli", full.names = TRUE)
            if (length(dd) == 1 && wn_exe_test(dd)) {
                myexe <- dd
            }
        }
        if (is.null(myexe)) {
            if (error) {
                stop("could not find the WindNinja command-line executable.\n You will need to install it yourself and either add it to the system path or pass the path to wn_find_exe.")
            }
            return(NULL)
        }
    }
    wn_opts$exe_path <- myexe
    options(windninjr = wn_opts)
    myexe
}

## internal: test a potential executable path
wn_exe_test <- function(path) {
    try({
        system(paste0(path, " --help"), intern = TRUE)
        return(TRUE)
    },silent = TRUE)
    FALSE
}

## adapted from http://conjugateprior.org/2015/06/identifying-the-os-from-r/
get_os <- function() {
    if (.Platform$OS.type=="windows") return("windows")
    sysinf <- Sys.info()
    if (!is.null(sysinf)){
        os <- sysinf["sysname"]
        if (tolower(os)=="darwin")
            os <- "osx"
    } else {
        os <- .Platform$OS.type
        if (grepl("^darwin", R.version$os,ignore.case=TRUE))
            os <- "osx"
        if (grepl("linux-gnu", R.version$os,ignore.case=TRUE))
            os <- "linux"
    }
    os <- tolower(os)
    if (!os %in% c("windows","linux","unix","osx"))
        stop("unknown operating system: ",os)
    os
}
