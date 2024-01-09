index <- function(info) {
  list(run = unlist(info$index["cases_inc"]), state = unlist(info$index))
}
