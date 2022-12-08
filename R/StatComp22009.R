#' @title A illustrate dataset
#' @name Eckerle4
#' @description A dataset used to illustrate the brownian distance covariance of \code{dcovR} and \code{dcovC}.
#' @examples
#' \dontrun{
#' data(Eckerle4)
#' attach(Eckerle4)
#' plot(Eckerle4)
#' }
NULL

#' @title Computing brownian distance covariance using R.
#' @description Computing brownian distance covariance of the input samples using Rcpp.
#' @param x a group of samples
#' @param y a group of samples
#' @return the brownian distance covariance of the two groups of samples
#' @examples
#' \dontrun{
#' data(Eckerle4)
#' attach(Eckerle4)
#' dcov <- dcovR(Eckerle4$wavelength, Eckerle4$transmitance)
#' }
#' @export
dcovR<-function(x,y){
  a <- as.matrix(dist(x))
  b <- as.matrix(dist(y))
  aa <- a
  bb <- b
  n <- ncol(a)
  a1 <- rowMeans(a)
  a2 <- colMeans(a)
  b1 <- rowMeans(b)
  b2 <- colMeans(b)
  for (i in 1:n){
    for (j in 1:n){
      aa[i,j] <- a1[i]+a2[j]
      bb[i,j] <- b1[i]+b2[j]
    }
  }
  A <- a + mean(a) - aa
  B <- b + mean(b) - bb
  return(sqrt(mean(A*B)))
}

#' @title Computing brownian distance correlation using R
#' @description Computing brownian distance correlation of the input samples using Rcpp
#' @param x a group of samples
#' @param y a group of samples
#' @return the brownian distance correlation of the two groups of samples
#' @examples
#' \dontrun{
#' data(Eckerle4)
#' attach(Eckerle4)
#' dcor <- dcor(Eckerle4$wavelength, Eckerle4$transmitance)
#' }
#' @export
dcor<-function(x, y){
  a <- dcovR(x, y)^2
  b <- dcovR(x, x) * dcovR(y, y)
  return(sqrt(a / b))
}

#' @title Computing the p-value using R
#' @description Computing the p-value of the brownian distance covariance, pearson covariance and spearman covariance using R
#' @param x a group of samples
#' @param y a group of samples
#' @param cov_type the tpye of covariance, cov() or dcov()
#' @param repl the times of replication
#' @param Method the tpye of covariance in cov(), pearson or spearman
#' @return the p-value
#' @examples
#' \dontrun{
#' data(Eckerle4)
#' attach(Eckerle4)
#' pvalue <- p_value(Eckerle4$wavelength, Eckerle4$transmitance, cov, 1000, 'pearson')
#' }
#' @export
p_value <- function(x, y, cov_type, repl, Method=NULL){
  a <- matrix(0, ncol = 2, nrow = length(x))
  a[, 1] <- x;a[, 2] <- y
  f <- function(data, i){
    m <- data[, 1];n <- data[i, 2]
    k <- nrow(data)
    ifelse(is.null(Method), return(k * cov_type(m, n)^2), return(k * cov_type(m, n, method = Method)^2))
  }
  c <- boot(a, statistic = f, R = repl)
  t <- c(c$t0, c$t)
  p <- length(t[t >= t[1]]) / length(t)
  return(p)
}


