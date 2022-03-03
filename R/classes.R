

#' function to check if object is of class stanfit
#' @param x any R object
#' @return TRUE if object is of class stanfit. FALSE otherwise
#'
#' @examples
#'  x <- c(1, 5, 4)
#'  is.stanfit(x)
#' @export is.stanfit
#'
is.stanfit <- function(x) {
  methods::is(x, 'stanfit')
}

#' function to check if object is of class stanfit optim
#' @param x any R object
#' @return TRUE if object is of class stanfitOptim. FALSE otherwise
#'
#' @examples
#'  x <- c(1, 5, 4)
#'  is.stanfit(x)
#' @export is.stanfitOptim
#'
is.stanfitOptim <- function(x) {
  methods::is(x, 'stanfitOptim')
}


#' function to check if object is of class BMADR and its subcalsses
#' @param x any R object
#' @return TRUE if object is of class BMADR. FALSE otherwise
#'
#' @examples
#'  x <- c(1, 5, 4)
#'  is.BMADR(x)
#' @export is.BMADR
#'
is.BMADR <- function(x) {
  inherits(x, c("BMADR", "LP", "BS"))
}


#' function to check if object is of class BMADR and its subcalsses
#' @param x any R object
#' @return TRUE if object is of class stanfitOptim. FALSE otherwise
#'
#' @examples
#'  x <- c(1, 5, 4)
#'  is.BMADR2(x)
#' @export is.BMADR2
#'
is.BMADR2 <- function(x) {
  inherits(x, c("BMADR", "LP", "BS"), TRUE)
}


#' function to check if object is of class BMADR and its subcalsses
#' @param x any R object
#' @return TRUE if object is of class BMADRQ. FALSE otherwise
#'
#' @examples
#'  x <- c(1, 5, 4)
#'  is.BMADRQ(x)
#' @export is.BMADRQ

is.BMADRQ <- function(x) {
  inherits(x, c("BMADRQ", "LP", "BS"))
}


#' function to check if object is of class BMADR and its subcalsses
#' @param x any R object
#' @return 1, 2 or 3 if object is of class BMABMDRQ using Laplace or BMABMDRQ using MCMC.
#'
#' @examples
#'  x <- c(1, 5, 4)
#'  is.BMADRQ2(x)
#' @export is.BMADRQ2

is.BMADRQ2 <- function(x) {
  inherits(x, c("BMADRQ", "LP", "BS"), TRUE)
}


