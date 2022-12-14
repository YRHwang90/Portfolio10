rmvn <- function(n,sigma) {
  Sh <- with(svd(sigma), v%*%diag(sqrt(d))%*%t(u))
  matrix(stats::rnorm(ncol(sigma)*n),ncol=ncol(sigma))%*%Sh
}
