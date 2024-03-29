% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zabp.R
\name{ZABP}
\alias{ZABP}
\alias{dZABP}
\alias{dZARBS}
\alias{pZARBS}
\alias{qZARBS}
\alias{rZARBS}
\alias{plotZARBS}
\alias{meanZARBS}
\alias{pZABP}
\alias{qZABP}
\alias{rZABP}
\title{Zero-Adjusted Beta-Prime (ZABP) distribution for fitting a GAMLSS}
\usage{
dZABP(x, mu = 1.0, sigma = 1.0, nu = 0.1, log = FALSE)

pZABP(q, mu = 1, sigma = 1, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

qZABP(p, mu = 1, sigma = 1, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

rZABP(n, mu = 1, sigma = 1, nu = 0.1)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{mu}{vector of scale parameter values.}

\item{sigma}{vector of shape parameter values.}

\item{nu}{vector of mixture parameter values.}

\item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
}
\description{
The functions dZABP, pZABP, qZABP and rZABP define
the density, distribution function, quantile function and random generation for
the ZABP.
}
\references{
Leiva, V., Santos-Neto, M., Cysneiros, F.J.A., Barros, M. (2016) A methodology for stochastic inventory models based on a zero-adjusted Birnbaum-Saunders distribution.
\emph{Applied Stochastic Models in Business and Industry.}, 32(1), 74--89. doi:\email{10.1002/asmb.2124}.

Santos-Neto, M., Cysneiros, F.J.A., Leiva, V., Barros, M. (2016) Reparameterized Birnbaum-Saunders
regression models with varying precision. \emph{Electronic Journal of Statistics}, 10, 2825--2855. doi: \email{10.1214/16-EJS1187}.
}
\author{
Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
}
