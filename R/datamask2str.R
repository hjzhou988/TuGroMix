datamask2str <- function(var) {
  require(dplyr)
  require(rlang)
  rlang::as_name(rlang::quo(expr = {{var}}))
}
