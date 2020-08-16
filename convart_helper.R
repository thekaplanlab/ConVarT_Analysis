
# Fast match operator

# !NOT RUN! #

`%fin%` <- function(x, table) {
  stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}
