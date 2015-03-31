#' Create a b-spline design matrix.
#'
#' @param x The points at which the b-spline is evaluated.
#' @param xl The lower bound of the domain.
#' @param xr The upper bound of the domain.
#' @param ndx The granularity of the basis.
#' @param bdeg The degree of the B-spline polynomial components.
#'
#' @author pschulam
#'
#' @export
#'
bspline <- function(x, xl, xr, ndx, bdeg) {
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg*dx, xr + bdeg*dx, by = dx)
  B <- splines::spline.des(knots, x, bdeg + 1, 0 * x)$design
  B
}

bspline_penalty <- function(p, ord) {
  D <- diag(p)
  for (k in 1:ord) D <- diff(D)
  t(D) %*% D
}