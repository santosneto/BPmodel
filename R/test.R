#'@name test
#'
#'@aliases grad_test_bp
#'@aliases grad_test_ga
#'@aliases grad_test_ig
#'@aliases grad_test_wei
#'@aliases grad_test_rbs
#'@aliases wald_test_bp
#'@aliases wald_test_ga
#'@aliases wald_test_ig
#'@aliases wald_test_wei
#'@aliases score_test_bp
#'@aliases score_test_ga
#'@aliases score_test_ig
#'@aliases score_test_wei
#'@aliases score_test_rbs
#'
#'@title Precision test
#'
#'@description Tests the null hypothesis of precision fixed in RBS models against the alternative of precision variable.
#'
#'@usage grad_test_bp(modelh0,modelh1)
#'
#' @param modelh0 model under null hypothesis.
#' @param modelh1 model under alternative hypothesis.
#'
#' @return A list with class "htest" containing the following components:
#' @return \code{statistic}	the value of the test statistic.
#' @return \code{parameter}	the degrees of freedom for the test statistic.
#' @return \code{p.value}	the p-value for the test.
#' @return \code{method}	a character string indicating what type of likelihood ratio test was performed.
#' @return \code{data.name} a character string giving the name(s) of the data
#'
#'@author
#'Manoel Santos-Neto \email{manoelferreira@uaest.ufcg.edu.br}
#'
#'
#'@importFrom stats pchisq
#'
#'@export

grad_test_bp <- function(modelh0, modelh1)
{
  UalphaH0.gam <- function(modelh0, modelh1) {
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dBP(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }

    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }

    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  METHOD <- "Gradient test"
  DNAME <- deparse(substitute(modelh0))
  DNAME <- paste(DNAME, "vs", deparse(substitute(modelh1)))
  Ua <- UalphaH0.gam(modelh0, modelh1)[-1]
  alphah1 <- modelh1$sigma.coefficients[-1]
  G <- t(Ua) %*% alphah1
  q <- modelh1$sigma.df
  gl <- q - 1
  PVAL <- stats::pchisq(G, q - 1, lower.tail = F)
  names(gl) <- "df"
  names(G) <- "G"
  RVAL <- list(statistic = G, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}


#'@rdname test
#'
#'@importFrom stats pchisq
#'
#'@export

grad_test_ga <- function(modelh0, modelh1)
{
  UalphaH0.gam <- function(modelh0, modelh1) {
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dGA(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }

    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }

    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  METHOD <- "Gradient test"
  DNAME <- deparse(substitute(modelh0))
  DNAME <- paste(DNAME, "vs", deparse(substitute(modelh1)))
  Ua <- UalphaH0.gam(modelh0, modelh1)[-1]
  alphah1 <- modelh1$sigma.coefficients[-1]
  G <- t(Ua) %*% alphah1
  q <- modelh1$sigma.df
  gl <- q - 1
  PVAL <- stats::pchisq(G, q - 1, lower.tail = F)
  names(gl) <- "df"
  names(G) <- "G"
  RVAL <- list(statistic = G, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}


#'@rdname test
#'
#'@importFrom stats pchisq
#'
#'@export
#'

grad_test_ig <- function(modelh0, modelh1)
{
  UalphaH0.gam <- function(modelh0, modelh1) {
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dIG(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }


    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }


    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  METHOD <- "Gradient test"
  DNAME <- deparse(substitute(modelh0))
  DNAME <- paste(DNAME, "vs", deparse(substitute(modelh1)))
  Ua <- UalphaH0.gam(modelh0, modelh1)[-1]
  alphah1 <- modelh1$sigma.coefficients[-1]
  G <- t(Ua) %*% alphah1
  q <- modelh1$sigma.df
  gl <- q - 1
  PVAL <- stats::pchisq(G, q - 1, lower.tail = F)
  names(gl) <- "df"
  names(G) <- "G"
  RVAL <- list(statistic = G, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}


#'@rdname test
#'
#'@importFrom stats pchisq
#'@export

grad_test_wei <- function(modelh0, modelh1)
{
  UalphaH0.gam <- function(modelh0, modelh1) {
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dWEI3(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }

    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }

    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  METHOD <- "Gradient test"
  DNAME <- deparse(substitute(modelh0))
  DNAME <- paste(DNAME, "vs", deparse(substitute(modelh1)))
  Ua <- UalphaH0.gam(modelh0, modelh1)[-1]
  alphah1 <- modelh1$sigma.coefficients[-1]
  G <- t(Ua) %*% alphah1
  q <- modelh1$sigma.df
  gl <- q - 1
  PVAL <- stats::pchisq(G, q - 1, lower.tail = F)
  names(gl) <- "df"
  names(G) <- "G"
  RVAL <- list(statistic = G, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

#'@rdname test
#'
#'@importFrom stats pchisq
#'@importFrom RBS dRBS

#'@export

grad_test_rbs <- function(modelh0, modelh1)
{
  UalphaH0.gam <- function(modelh0, modelh1) {
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dRBS(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }


    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }

    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  METHOD <- "Gradient test"
  DNAME <- deparse(substitute(modelh0))
  DNAME <- paste(DNAME, "vs", deparse(substitute(modelh1)))
  Ua <- UalphaH0.gam(modelh0, modelh1)[-1]
  alphah1 <- modelh1$sigma.coefficients[-1]
  G <- t(Ua) %*% alphah1
  q <- modelh1$sigma.df
  gl <- q - 1
  PVAL <- stats::pchisq(G, q - 1, lower.tail = F)
  names(gl) <- "df"
  names(G) <- "G"
  RVAL <- list(statistic = G, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

#'@rdname test
#'
#'@importFrom stats pchisq
#'@importFrom Deriv Deriv
#'
#'@export

score_test_bp <- function(modelh0, modelh1)
{
  UalphaH0.gam <- function(modelh0, modelh1){
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dBP(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }


    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }


    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  vcovH0 <- function(modelh0, modelh1){
    mu <- modelh0$mu.fv
    sigma <-  modelh0$sigma.fv
    x <- modelh1$mu.x
    z <- modelh1$sigma.x
    linkstr <- modelh0$mu.link
    linkobj <- make.link(linkstr)
    mu.eta <- linkobj$mu.eta
    eta <- linkobj$linkfun
    etaH0 <- eta(modelh0$mu.fv)
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    y <- modelh0$y

    dai <- function(link)
    {
      switch(link, log = {mu.eta.2 <- function(eta) rep.int(1, length(eta))},
             identity = {mu.eta.2 <- function(eta) rep.int(0, length(eta))},
             sqrt = {mu.eta.2 <- function(eta) 1/eta} )
    }
    mu.eta.2 <- dai(linkstr)
    sigma.eta.2 <- dai(phi_linkstr)
    d_ai <- mu.eta.2(etaH0)
    d_bi <- sigma.eta.2(tauH0)
    ai <- mu.eta(etaH0)
    bi <- phi_mu.eta(tauH0)

    ll <- function(y,mu,sigma){
      a <- mu * (1 + sigma)
      b <- 2 + sigma
      fy <- (a - 1) * log(y) - (a + b) * log(1 + y) - lbeta(a, b)

      fy
    }


    dm <- Deriv::Deriv(ll,'mu')
    ds <- Deriv::Deriv(ll,'sigma')
    dmm <- Deriv::Deriv(Deriv(ll,'mu'),'mu')
    dss <- Deriv::Deriv(Deriv(ll,'sigma'),'sigma')
    dms <- Deriv::Deriv(Deriv(ll,'mu'),'sigma')

    d_mu <- dm(y,mu,sigma)
    d_sigma <- ds(y,mu,sigma)
    v <- dmm(y,mu,sigma)
    u <- dss(y,mu,sigma)
    s <- dms(y,mu,sigma)

    ci <- v*(ai^2) + d_mu*ai*d_ai
    mi <- s*ai*bi
    wi <- u*(bi^2) + d_sigma*bi*d_bi

    kbb <- crossprod(ci * x, x)
    kaa <- crossprod(wi * z, z)
    kba <- crossprod(mi * x, z)
    hess <- cbind(rbind(kbb, t(kba)), rbind(kba, kaa))
    vcov <- solve(-hess)
    return(vcov)
  }

  METHOD <- "Rao score test"
  DNAME <- deparse(substitute(modelh0))
  DNAME <- paste(DNAME, "vs", deparse(substitute(modelh1)))
  p <- modelh0$mu.df
  p1 <- p + 1
  varalpha <- vcovH0(modelh0, modelh1)[-(1:p1), -(1:p1)]
  Ua. <- UalphaH0.gam(modelh0, modelh1)[-1]
  SC <- t(Ua.) %*% varalpha %*% Ua.
  q <- modelh1$sigma.df
  gl <- q - 1
  PVAL <- pchisq(SC, q - 1, lower.tail = F)
  names(gl) <- "df"
  names(SC) <- "SC"
  RVAL <- list(statistic = SC, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

#'@rdname test
#'
#'@importFrom stats pchisq
#'@export

score_test_rbs <- function(modelh0, modelh1)
{


  UalphaH0.gam <- function(modelh0, modelh1){
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dRBS(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }


    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }


    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  vcovH0 <- function(modelh0, modelh1) {
    mu <- modelh0$mu.fv
    sigma <-  modelh0$sigma.fv
    x <- modelh1$mu.x
    z <- modelh1$sigma.x
    linkstr <- modelh0$mu.link
    linkobj <- make.link(linkstr)
    mu.eta <- linkobj$mu.eta
    eta <- linkobj$linkfun
    etaH0 <- eta(modelh0$mu.fv)
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    y <- modelh0$y

    dai <- function(link)
    {
      switch(link, log = {mu.eta.2 <- function(eta) rep.int(1, length(eta))},
             identity = {mu.eta.2 <- function(eta) rep.int(0, length(eta))},
             sqrt = {mu.eta.2 <- function(eta) 1/eta} )
    }
    mu.eta.2 <- dai(linkstr)
    sigma.eta.2 <- dai(phi_linkstr)
    d_ai <- mu.eta.2(etaH0)
    d_bi <- sigma.eta.2(tauH0)
    ai <- mu.eta(etaH0)
    bi <- phi_mu.eta(tauH0)

    ll <- function(y,mu,sigma){
      fy <-  0.5*sigma - 0.5*log((sigma + 1)) - 0.5 * log(mu) - 1.5 * log(y) + log((sigma * y) + y + (sigma * mu)) - (y * (sigma + 1))/(4 * mu) - (sigma*sigma*mu)/(4*y * (sigma + 1)) - 0.5 * log(16*pi)

      fy
    }

    dm <- Deriv::Deriv(ll,'mu')
    ds <- Deriv::Deriv(ll,'sigma')
    dmm <- Deriv::Deriv(Deriv(ll,'mu'),'mu')
    dss <- Deriv::Deriv(Deriv(ll,'sigma'),'sigma')
    dms <- Deriv::Deriv(Deriv(ll,'mu'),'sigma')

    d_mu <- dm(y,mu,sigma)
    d_sigma <- ds(y,mu,sigma)
    v <- dmm(y,mu,sigma)
    u <- dss(y,mu,sigma)
    s <- dms(y,mu,sigma)

    ci <- v*(ai^2) + d_mu*ai*d_ai
    mi <- s*ai*bi
    wi <- u*(bi^2) + d_sigma*bi*d_bi

    kbb <- crossprod(ci * x, x)
    kaa <- crossprod(wi * z, z)
    kba <- crossprod(mi * x, z)
    hess <- cbind(rbind(kbb, t(kba)), rbind(kba, kaa))
    vcov <- solve(-hess)
    return(vcov)
  }

  METHOD <- "Rao score test"
  DNAME <- deparse(substitute(modelh0))
  DNAME <- paste(DNAME, "vs", deparse(substitute(modelh1)))
  p <- modelh0$mu.df
  p1 <- p + 1
  varalpha <- round(vcovH0(modelh0, modelh1),5)[-(1:p1), -(1:p1)]
  Ua. <- UalphaH0.gam(modelh0, modelh1)[-1]
  SC <- t(Ua.) %*% varalpha %*% Ua.
  q <- modelh1$sigma.df
  gl <- q - 1
  PVAL <- pchisq(SC, q - 1, lower.tail = F)
  names(gl) <- "df"
  names(SC) <- "SC"
  RVAL <- list(statistic = SC, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}



#'@rdname test
#'
#'@importFrom stats pchisq
#'@export


score_test_ga <- function(modelh0, modelh1)
{


  UalphaH0.gam <- function(modelh0, modelh1){
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dGA(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }


    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }


    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  vcovH0 <- function(modelh0, modelh1) {
    mu <- modelh0$mu.fv
    sigma <-  modelh0$sigma.fv
    x <- modelh1$mu.x
    z <- modelh1$sigma.x
    linkstr <- modelh0$mu.link
    linkobj <- make.link(linkstr)
    mu.eta <- linkobj$mu.eta
    eta <- linkobj$linkfun
    etaH0 <- eta(modelh0$mu.fv)
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    y <- modelh0$y

    dai <- function(link)
    {
      switch(link, log = {mu.eta.2 <- function(eta) rep.int(1, length(eta))},
             identity = {mu.eta.2 <- function(eta) rep.int(0, length(eta))},
             sqrt = {mu.eta.2 <- function(eta) 1/eta} )
    }
    mu.eta.2 <- dai(linkstr)
    sigma.eta.2 <- dai(phi_linkstr)
    d_ai <- mu.eta.2(etaH0)
    d_bi <- sigma.eta.2(tauH0)
    ai <- mu.eta(etaH0)
    bi <- phi_mu.eta(tauH0)

    ll <- function(y,mu,sigma){
      fy <-  (1/sigma^2) * log(y/(mu * sigma^2)) - y/(mu*sigma^2) - log(y) - lgamma(1/sigma^2)

      fy
    }

    dm <- Deriv::Deriv(ll,'mu')
    ds <- Deriv::Deriv(ll,'sigma')
    dmm <- Deriv::Deriv(Deriv(ll,'mu'),'mu')
    dss <- Deriv::Deriv(Deriv(ll,'sigma'),'sigma')
    dms <- Deriv::Deriv(Deriv(ll,'mu'),'sigma')

    d_mu <- dm(y,mu,sigma)
    d_sigma <- ds(y,mu,sigma)
    v <- dmm(y,mu,sigma)
    u <- dss(y,mu,sigma)
    s <- dms(y,mu,sigma)

    ci <- v*(ai^2) + d_mu*ai*d_ai
    mi <- s*ai*bi
    wi <- u*(bi^2) + d_sigma*bi*d_bi

    kbb <- crossprod(ci * x, x)
    kaa <- crossprod(wi * z, z)
    kba <- crossprod(mi * x, z)
    hess <- cbind(rbind(kbb, t(kba)), rbind(kba, kaa))
    vcov <- solve(-hess)
    return(vcov)
  }

  METHOD <- "Rao score test"
  DNAME <- deparse(substitute(modelh0))
  DNAME <- paste(DNAME, "vs", deparse(substitute(modelh1)))
  p <- modelh0$mu.df
  p1 <- p + 1
  varalpha <- vcovH0(modelh0, modelh1)[-(1:p1), -(1:p1)]
  Ua. = UalphaH0.gam(modelh0, modelh1)[-1]
  SC = t(Ua.) %*% varalpha %*% Ua.
  q = modelh1$sigma.df
  gl = q - 1
  PVAL = pchisq(SC, q - 1, lower.tail = F)
  names(gl) = "df"
  names(SC) = "SC"
  RVAL <- list(statistic = SC, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

#'@rdname test
#'
#'@importFrom stats pchisq
#'@export

score_test_ig <- function(modelh0, modelh1)
{


  UalphaH0.gam <- function(modelh0, modelh1){
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dIG(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }


    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }


    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  vcovH0 <- function(modelh0, modelh1) {
    mu <- modelh0$mu.fv
    sigma <-  modelh0$sigma.fv
    x <- modelh1$mu.x
    z <- modelh1$sigma.x
    linkstr <- modelh0$mu.link
    linkobj <- make.link(linkstr)
    mu.eta <- linkobj$mu.eta
    eta <- linkobj$linkfun
    etaH0 <- eta(modelh0$mu.fv)
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    y <- modelh0$y

    dai <- function(link)
    {
      switch(link, log = {mu.eta.2 <- function(eta) rep.int(1, length(eta))},
             identity = {mu.eta.2 <- function(eta) rep.int(0, length(eta))},
             sqrt = {mu.eta.2 <- function(eta) 1/eta} )
    }
    mu.eta.2 <- dai(linkstr)
    sigma.eta.2 <- dai(phi_linkstr)
    d_ai <- mu.eta.2(etaH0)
    d_bi <- sigma.eta.2(tauH0)
    ai <- mu.eta(etaH0)
    bi <- phi_mu.eta(tauH0)

    ll <- function(y,mu,sigma){
      fy <-  (-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y) - ((y - mu)^2)/(2 * sigma^2 * (mu^2) * y))

      fy
    }

    dm <- Deriv::Deriv(ll,'mu')
    ds <- Deriv::Deriv(ll,'sigma')
    dmm <- Deriv::Deriv(Deriv(ll,'mu'),'mu')
    dss <- Deriv::Deriv(Deriv(ll,'sigma'),'sigma')
    dms <- Deriv::Deriv(Deriv(ll,'mu'),'sigma')

    d_mu <- dm(y,mu,sigma)
    d_sigma <- ds(y,mu,sigma)
    v <- dmm(y,mu,sigma)
    u <- dss(y,mu,sigma)
    s <- dms(y,mu,sigma)

    ci <- v*(ai^2) + d_mu*ai*d_ai
    mi <- s*ai*bi
    wi <- u*(bi^2) + d_sigma*bi*d_bi

    kbb <- crossprod(ci * x, x)
    kaa <- crossprod(wi * z, z)
    kba <- crossprod(mi * x, z)
    hess <- cbind(rbind(kbb, t(kba)), rbind(kba, kaa))
    vcov <- solve(-hess)
    return(vcov)
  }

  METHOD = "Rao score test"
  DNAME = deparse(substitute(modelh0))
  DNAME = paste(DNAME, "vs", deparse(substitute(modelh1)))
  p = modelh0$mu.df
  p1 = p + 1
  varalpha = vcovH0(modelh0, modelh1)[-(1:p1), -(1:p1)]
  Ua. = UalphaH0.gam(modelh0, modelh1)[-1]
  SC = t(Ua.) %*% varalpha %*% Ua.
  q = modelh1$sigma.df
  gl = q - 1
  PVAL = pchisq(SC, q - 1, lower.tail = F)
  names(gl) = "df"
  names(SC) = "SC"
  RVAL <- list(statistic = SC, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}


#'@rdname test
#'
#'@importFrom stats pchisq
#'@export

score_test_wei <- function(modelh0, modelh1)
{


  UalphaH0.gam <- function(modelh0, modelh1){
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dWEI3(x = y,mu = theta[1],sigma = theta[2],log = TRUE)
      return(L)
    }


    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }


    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  vcovH0 <- function(modelh0, modelh1) {
    mu <- modelh0$mu.fv
    sigma <-  modelh0$sigma.fv
    x <- modelh1$mu.x
    z <- modelh1$sigma.x
    linkstr <- modelh0$mu.link
    linkobj <- make.link(linkstr)
    mu.eta <- linkobj$mu.eta
    eta <- linkobj$linkfun
    etaH0 <- eta(modelh0$mu.fv)
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    y <- modelh0$y

    dai <- function(link)
    {
      switch(link, log = {mu.eta.2 <- function(eta) rep.int(1, length(eta))},
             identity = {mu.eta.2 <- function(eta) rep.int(0, length(eta))},
             sqrt = {mu.eta.2 <- function(eta) 1/eta} )
    }
    mu.eta.2 <- dai(linkstr)
    sigma.eta.2 <- dai(phi_linkstr)
    d_ai <- mu.eta.2(etaH0)
    d_bi <- sigma.eta.2(tauH0)
    ai <- mu.eta(etaH0)
    bi <- phi_mu.eta(tauH0)

    ll <- function(y,mu,sigma){
      mu2 <- mu/gamma((1/sigma) + 1)
      fy <-  log(sigma) + (sigma-1)*log(y) - sigma*log(mu2) - (y/mu2)^sigma

      fy
    }

    dm <- Deriv::Deriv(ll,'mu')
    ds <- Deriv::Deriv(ll,'sigma')
    dmm <- Deriv::Deriv(Deriv(ll,'mu'),'mu')
    dss <- Deriv::Deriv(Deriv(ll,'sigma'),'sigma')
    dms <- Deriv::Deriv(Deriv(ll,'mu'),'sigma')

    d_mu <- dm(y,mu,sigma)
    d_sigma <- ds(y,mu,sigma)
    v <- dmm(y,mu,sigma)
    u <- dss(y,mu,sigma)
    s <- dms(y,mu,sigma)

    ci <- v*(ai^2) + d_mu*ai*d_ai
    mi <- s*ai*bi
    wi <- u*(bi^2) + d_sigma*bi*d_bi

    kbb <- crossprod(ci * x, x)
    kaa <- crossprod(wi * z, z)
    kba <- crossprod(mi * x, z)
    hess <- cbind(rbind(kbb, t(kba)), rbind(kba, kaa))
    vcov <- solve(-hess)
    return(vcov)
  }

  METHOD = "Rao score test"
  DNAME = deparse(substitute(modelh0))
  DNAME = paste(DNAME, "vs", deparse(substitute(modelh1)))
  p = modelh0$mu.df
  p1 = p + 1
  varalpha = vcovH0(modelh0, modelh1)[-(1:p1), -(1:p1)]
  Ua. = UalphaH0.gam(modelh0, modelh1)[-1]
  SC = t(Ua.) %*% varalpha %*% Ua.
  q = modelh1$sigma.df
  gl = q - 1
  PVAL = pchisq(SC, q - 1, lower.tail = F)
  names(gl) = "df"
  names(SC) = "SC"
  RVAL <- list(statistic = SC, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

#'@rdname test
#'
#'@importFrom stats pchisq
#'@export

wald_test <- function(modelh1)
{
  METHOD = "Wald test"
  DNAME = deparse(substitute(modelh1))
  p = modelh1$mu.df
  p1 = p + 1
  vcov = vcov(modelh1)
  varalpha = vcov[-(1:p1), -(1:p1)]
  alphah1 = modelh1$sigma.coefficients[-1]
  W = t(alphah1) %*% solve(varalpha) %*% alphah1
  q = modelh1$sigma.df
  gl = q - 1
  PVAL = pchisq(W, q - 1, lower.tail = F)
  names(gl) = "df"
  names(W) = "W"
  RVAL <- list(statistic = W, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}















