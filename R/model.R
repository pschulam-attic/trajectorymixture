logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' Fit a mixture of mean trajectories.
#'
#' @param trajectories Data frame containing trajectory ids in the first column, observation times in the second, and observation values in the third.
#' @param K The number of components in the mixture.
#' @param bfn The basis function mapping times to a design matrix.
#' @param maxiter The maximum number of EM iterations.
#' @param eps The smallest relative change in log likelihood needed to stop.
#'
#' @author pschulam
#'
#' @export
#'
trajectorymixture <- function(trajectories, K, bfn,
                              maxiter = 100, eps = 1e-5) {

  idcol <- 1
  timecol <- 2
  valuecol <- 3

  id <- unique(trajectories[[idcol]])
  N <- length(id)

  g <- trajectories[[idcol]]
  x <- trajectories[[timecol]]
  y <- trajectories[[valuecol]]

  B <- bfn(x)
  X <- cbind(1, B)
  P <- matrix(0, ncol(X), ncol(X))
  P[-1, -1] <- bspline_penalty(ncol(B), 1)
  P <- P + diag(1e-4, ncol(X))

  z <- matrix(runif(N * K), N, K)
  z <- sweep(z, 1, rowSums(z), "/")

  probs <- rep(0, K)
  mean_coef <- matrix(0, ncol(X), K)
  variance <- 1.0
  logl <- -Inf

  for (itr in 1:maxiter) {
    rss <- 0.0
    ss1 <- array(0, c(ncol(X), ncol(X), K))
    ss2 <- array(0, c(ncol(X), K))

    probs[] <- colSums(z) / sum(z)

    for (i in id) {
      ix <- g == i
      X_i <- X[ix, ]
      y_i <- y[ix]
      z_i <- z[id == i, ]

      for (j in 1:K) {
        yhat <- as.vector(X_i %*% mean_coef[, j])
        rss <- rss + z_i[j] * sum((y_i - yhat) ^ 2)
        ss1[, , j] <- ss1[, , j] + z_i[j] * t(X_i) %*% X_i
        ss2[, j] <- ss2[, j] + z_i[j] * t(X_i) %*% y_i
      }
    }

    variance[] <- rss / nrow(trajectories)

    for (j in 1:K) {
      mean_coef[, j] <- solve(ss1[, , j], ss2[, j])
    }

    logl_old <- logl
    logl <- 0.0

    for (i in id) {
      ix <- g == i
      X_i <- X[ix, ]
      y_i <- y[ix]
      ll <- log(probs)

      for (j in 1:K) {
        m <- as.vector(X_i %*% weights[, j])
        V <- variance * diag(length(y_i))
        ll[j] <- ll[j] + mvtnorm::dmvnorm(y_i, m, V, log = TRUE)
      }

      ls <- logsumexp(ll)
      z[id == i, ] <- exp(ll - ls)
      logl <- logl + ls
    }

    rel_diff <- (logl - logl_old) / abs(logl)
    convergence <- rel_diff < eps
    if (convergence) break
  }

  list(probs, mean_coef, variance)
}
