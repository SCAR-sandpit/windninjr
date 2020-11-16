## adapted from base::match.arg
match_arg <- function(arg, choices, several_ok = FALSE, case_insensitive = TRUE) {
    argname <- sub("^.*\\$", "", deparse(substitute(arg)))
    if (!(is.character(arg) && length(arg) == 1)) stop(argname, " must be a string (a character vector of length 1)")
    if (case_insensitive) {
        choices <- tolower(choices)
        arg <- tolower(arg)
    }
    if (!several_ok) {
        if (identical(arg, choices)) 
            return(arg[1L])
        if (length(arg) > 1L) 
            stop(argname, " must be of length 1")
    } else if (length(arg) == 0L) {
        stop(argname, " must be of length >= 1")
    }
    i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L)) stop(sprintf("%s should be one of %s", argname, paste(dQuote(choices), collapse = ", ")))
    i <- i[i > 0L]
    if (!several_ok && length(i) > 1) stop(sprintf("%s matches more than one of the possible choices %s", argname, paste(dQuote(choices), collapse = ", ")))
    choices[i]
}
