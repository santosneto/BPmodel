#'@name ZABP
#'
#'@aliases dZARBS
#'@aliases pZARBS
#'@aliases qZARBS
#'@aliases rZARBS
#'@aliases plotZARBS
#'@aliases meanZARBS
#'
#'@title Zero-Adjusted Beta-Prime (ZABP) distribution for fitting a GAMLSS
#'
#'@description The functions dZABP, pZABP, qZABP and rZABP define
#'the density, distribution function, quantile function and random generation for
#'the ZABP.
#'
#'@usage dZABP(x, mu = 1.0, sigma = 1.0, nu = 0.1, log = FALSE)
#'
#' @param x,q vector of quantiles.
#' @param mu vector of scale parameter values.
#' @param sigma vector of shape parameter values.
#' @param nu vector of mixture parameter values.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A., Barros, M. (2016) A methodology for stochastic inventory models based on a zero-adjusted Birnbaum-Saunders distribution.
#'\emph{Applied Stochastic Models in Business and Industry.}, 32(1), 74--89. doi:\email{10.1002/asmb.2124}.
#'
#'Santos-Neto, M., Cysneiros, F.J.A., Leiva, V., Barros, M. (2016) Reparameterized Birnbaum-Saunders
#'regression models with varying precision. \emph{Electronic Journal of Statistics}, 10, 2825--2855. doi: \email{10.1214/16-EJS1187}.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@export

dZABP <- function(x,
                  mu = 1.0,
                  sigma = 1.0,
                  nu = 0.1,
                  log = FALSE){

  if (any(mu < 0) || any(is.na(mu)))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0) || any(is.na(sigma)))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0) |  any(nu >= 1) )  stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(x < 0))  stop(paste("x must be positive", "\n", ""))


  aux <- function(x, mu, sigma, nu) return(ifelse(x==0, log(nu), dBP(x, mu, sigma, log = TRUE)))

  log.lik <- mapply(aux, x = x, mu = mu, sigma = sigma, nu = nu)

  if(log == FALSE) fy <- exp(log.lik) else fy <- log.lik

}

#'@rdname ZABP
#'
#'@export
pZABP <- function(q,
                  mu = 1.0,
                  sigma = 1.0,
                  nu = 0.1,
                  lower.tail = TRUE,
                  log.p = FALSE){

  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0) |  any(nu >= 1) )  stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(q < 0))  stop(paste("y must be positive", "\n", ""))

  cdf0 <- mapply(pBP, q, mu, sigma)

  auxcdf <- function(q, nu) return(ifelse((q==0), nu, nu + (1-nu)*cdf0))

  cdf <- mapply(auxcdf, q = q, nu = nu)

  if(lower.tail == TRUE) cdf  <- cdf else  cdf <- 1-cdf
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf)
  cdf

}


#'@rdname ZABP
#'
#'@export
qZABP <- function(p,
                  mu = 1.0,
                  sigma = 1.0,
                  nu = 0.1,
                  lower.tail = TRUE,
                  log.p = FALSE){

  if (any(mu <= 0)) stop(paste("mu must be positive ", "\n", ""))
  if (any(sigma < 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)|any(nu >= 1)) stop(paste("nu must be beetwen 0 and 1 ", "\n", ""))
  if (any(p < 0) | any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))

  if (log.p == TRUE) p <- exp(p) else p <- p
  if (lower.tail == TRUE) p <- p else p <- 1 - p

  p0 <- (p-nu)/(1-nu)

  aux0 <- function(p0, mu, sigma) return(ifelse(p0 <=0, 0,  qBP(p0, mu, sigma)))

  q <- mapply(aux0, p0 = p0, mu = mu, sigma = sigma)

  q
}


#'@rdname ZABP
#'
#'
#'@importFrom stats runif
#'@export
rZABP <- function(n,
                  mu = 1.0,
                  sigma = 1.0,
                  nu = 0.1){

  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)|any(nu >= 1)) stop(paste("nu must be beetwen 0 and 1 ", "\n", ""))
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))



  r <- mapply(qZABP, p = runif(n), mu = mu, sigma = sigma, nu = nu)

  r
}





