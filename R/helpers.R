check_data <- function(W, A, Y) {
  if(any(is.na(W)) || any(is.na(A)) || any(is.na(Y))) {
    stop("NA values found in either W, A or Y. Please make sure there are no missing values.")
  }
  if(!all(A %in% c(0,1))) {
    stop("`A` must have values in 0,1.")
  }
  if(any(Y < 0)) {
    stop("Y must be nonnegative for relative risk minimization to be possible using this package.")
  }
}



quotemeta <- function(string) {
  str_replace_all(string, "(\\W)", "\\\\\\1")
}
