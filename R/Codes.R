#'@name BP
#'
#'@aliases  BP
#'@aliases  dBP
#'@aliases  pBP
#'@aliases  qBP
#'@aliases  rBP
#'@aliases  hBP
#'@aliases  plotBP
#'@aliases  meanBP
#'
#'@title Reparameterized Beta Prime (BP) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{BP()} defines the BP distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dBP}, \code{pBP}, \code{qBP},
#'\code{rBP} and \code{hBP} define the density, distribution function, quantile function, random
#'genetation and hazard function for the \code{BP} parameterization of the RWL distribution.
#'
#'@usage BP(mu.link = "log", sigma.link = "log")
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log.d, log.p logical; if TRUE, density d are given as log(d).
#' @param log.p, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param log.q, log.p logical; if TRUE, quantiles q are given as log(q).
#' @param log.h, log.p logical; if TRUE, hazard h are given as log(h).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from.
#' @param to up to where to plot the distribution.
#' @param obj a fitted RWL object.
#' @param ... other graphical parameters for plotting.
#' @param title title of the plot.
#'
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a reparameterized weighted Lindley distribution in the \code{gamlss()} function.
#'
#'@note For the function RWL(), mu is the mean and sigma is the precision parameter of the reparameterized weighted Lindley distribution.
#'
#'@author
#'Manoel Santos-Neto \email{santosnetoce at protemail.com}
#'
#'@importFrom gamlss.dist checklink
#'
#'@export

BP <- function (mu.link = "log", sigma.link = "log")
{
  mstats <- checklink("mu.link", "Beta Prime", substitute(mu.link), c("log", "identity", "sqrt"))
  dstats <- checklink("sigma.link", "Beta Prime",substitute(sigma.link), c("log", "identity", "sqrt"))
  structure(list(family = c("BP", "Beta Prime"),
                 parameters = list(mu = TRUE,sigma = TRUE), nopar = 2, type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,


                 dldm = function(y, mu, sigma){
                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2
                   Phi <-  (1+sigma)
                   yast <- log(y) - log(1+y)
                   muast <- digamma(a) - digamma(b)
                   dldm <- Phi*(yast - muast)

                   dldm
                 },
                 d2ldm2 = function(mu, sigma){
                   Phi2 <- (1+sigma)^2
                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2
                   d2dldm2 <- -Phi2*(trigamma(a) - trigamma(b))

                   d2dldm2
                 },
                 dldd = function(y, mu, sigma){
                   Phi <-  (1+sigma)
                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2
                   ystar <- mu*log(y) - (1+mu)*log(1+y)
                   mustar <- mu*digamma(a) - (1+mu)*digamma(b) + digamma(Phi+1)

                   dldd <- ystar - mustar

                   dldd
                 },
                 d2ldd2 = function(mu,sigma){
                   Phi <-  (1+sigma)
                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2

                   d2ldd2 <- -(mu^2)*trigamma(a) + ((1+mu)^2)*trigamma(b) - trigamma(Phi+1)

                   d2ldd2

                 },
                 d2ldmdd = function(mu,sigma){

                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2
                   Phi <-  (1+sigma)
                   gammaast <- Phi*(trigamma(b) + mu*(trigamma(b)-trigamma(a)))

                   d2ldmdd <- gammaast

                   d2ldmdd

                 },
                 G.dev.incr = function(y, mu, sigma,...){-2*dBP(y, mu, sigma, log = TRUE)},
                 rqres = expression(rqres(pfun = "pBP", type = "Continuous", y = y, mu = mu, sigma = sigma)),
                 mu.initial = expression({mu <- mean(y)}),
                 sigma.initial = expression({sigma <-  mean(y)*(1+mean(y))/var(y) }),
                 mu.valid = function(mu) all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0)),
            class = c("gamlss.family","family"))
}

#'@rdname BP
#'@export

dBP <- function(x,mu=1,sigma=1,log=FALSE)
{
  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(x <= 0))
    stop(paste("x must be positive", "\n", ""))

  a <- mu*(1+sigma)
  b <- 2 + sigma

  fy <- extraDistr::dbetapr(x, shape1 = a, shape2 = b, scale = 1, log = log)
  fy

}

#'@rdname BP
#'
#'@export

pBP <-  function(q,mu=1,sigma=1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(q < 0))
    stop(paste("q must be positive", "\n", ""))

  a <- mu*(1+sigma)
  b <- 2 + sigma

  cdf <- extraDistr::pbetapr(q, shape1 = a, shape2 = b, scale=1, lower.tail = lower.tail,
                 log.p = log.p)
  cdf
}

#'@rdname BP
#'
#'@export

rBP <- function(n,mu=1,sigma=1)
{
  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(n <= 0))
    stop(paste("n must be a positive integer", "\n", ""))

  n <- ceiling(n)

  a <- mu*(1+sigma)
  b <- 2 + sigma

  r <- extraDistr::rbetapr(n,shape1=a,shape2=b,scale=1)

  r
}

#'@rdname BP
#'
#'@export

qBP <- function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(n <= 0))
    if (any(p <= 0) | any(p >= 1))
      stop(paste("p must be between 0 and 1", "\n", ""))

  a <- mu*(1+sigma)
  b <- 2 + sigma

  q <- extraDistr::qbetapr(p, shape1 = a, shape2 = b,scale=1, lower.tail = lower.tail, log.p = log.p)

  q
}


#'@name diag.RBS
#'
#'@aliases diag.RBS
#'@aliases residuals.pearson
#'
#'@title Diagnostic Analysis - Local Influnce
#'
#'@description Diagnostics for the BP model
#'
#'@param model Object of class \code{gamlss} holding the fitted model.
#'@param scheme Default is "case.weight". But, can be "response", "location" or "precision".
#'@param mu.link  Defines the mu.link, with "identity" link as the default for the mu parameter.
#'@param sigma.link Defines the sigma.link, with "identity" link as the default for the sigma parameter

#'
#'@return Local influence measures.
#'
#' @author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@importFrom pracma hessian ones
#'@importFrom stats make.link sd
#'@export

diag.BP <- function(model,mu.link = "log",sigma.link = "log",scheme="case.weight")
{

  x <- model$mu.x
  z  <- model$sigma.x
  y <- model$y
  p<-ncol(x)
  q<-ncol(z)

  linkstr <- mu.link
  linkobj <- make.link(linkstr)
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta  <- linkobj$mu.eta

  sigma_linkstr <- sigma.link
  sigma_linkobj <- make.link(sigma_linkstr)
  sigma_linkfun <- sigma_linkobj$linkfun
  sigma_linkinv <- sigma_linkobj$linkinv
  sigma_mu.eta  <- sigma_linkobj$mu.eta


  B=function(Delta,I,M)
  {
    B=(t(Delta)%*%(I-M)%*%Delta)
    return(B)
  }

  loglik <- function(vP)
  {
    betab = vP[1:p]
    alpha = vP[-(1:p)]
    eta   = as.vector(x%*%betab)
    tau   = as.vector(z%*%alpha)
    mu    = linkinv(eta)
    sigma = sigma_linkinv(tau)

    a <- mu*(1+sigma)
    b <- 2 + sigma

    fy <- (a-1)*log(y) - (a+b)*log(1+y) - lbeta(a,b)

    return(sum(fy))
  }

  muest <- model$mu.coefficients
  sigmaest <- model$sigma.coefficients
  x0<- c(muest,sigmaest)
  h0 <- hessian(loglik,x0)

  Ldelta= h0[(p+1):(p+q),(p+1):(p+q)]
  Lbeta=h0[1:p,1:p]
  b11=cbind(matrix(0, p, p), matrix(0, p, q))
  b12=cbind(matrix(0, q, p), solve(Ldelta))
  B1= rbind(b11, b12)  #parameter beta
  b211 =cbind(solve(Lbeta), matrix(0, p, q))
  b212= cbind(matrix(0, q, p), matrix(0, q, q))
  B2=rbind(b211,b212)  # parameter delta

  b311 =cbind(matrix(0, p, p), matrix(0, p, q))
  b312= cbind(matrix(0, q, p), matrix(0, q, q))
  B3=rbind(b311,b312)  # parameter theta

  if(scheme=="case.weight")
  {
    ############################Case Weight####################################

    mu <- model$mu.fv
    sigma <- model$sigma.fv
    eta <- linkfun(mu)
    ai <- mu.eta(eta)
    a <- mu*(1+sigma)
    b <- mu*(1+sigma)+sigma+2
    Phi <-  (1+sigma)
    yast <- log(y) - log(1+y)
    muast <- digamma(a) - digamma(b)
    dldm <- Phi*(yast - muast)
    Deltamu <- crossprod(x,diag(ai*dldm))


    tau <- sigma_linkfun(sigma)
    bi <- sigma_mu.eta(tau)
    ystar <- mu*log(y) - (1+mu)*log(1+y)
    mustar <- mu*digamma(a) - (1+mu)*digamma(b) + digamma(Phi+1)
    dldd <- ystar - mustar

    Deltasigma <- crossprod(z,diag(bi*dldd))

    Delta <- rbind(Deltamu,Deltasigma)

    ##################theta#########################
    BT<-B(Delta,solve(h0),B3)
    autovmaxthetaPC<- eigen(BT,symmetric=TRUE)$val[1]
    vetorpcthetaPC<- eigen(BT,symmetric=TRUE)$vec[,1]
    dmaxG.theta<-abs(vetorpcthetaPC)
    vCithetaPC<-2*abs(diag(BT))
    Cb0<-vCithetaPC
    Cb.theta<-Cb0/sum(Cb0)
    ######################betas########################
    BM<-B(Delta,solve(h0),B1)
    autovmaxbetaPC<-eigen(BM,symmetric=TRUE)$val[1]
    vetorpcbetaPC<-eigen(BM,symmetric=TRUE)$vec[,1]
    dmaxG.beta<-abs(vetorpcbetaPC)
    vCibetaPC<-2*abs(diag(BM))
    Cb1<-vCibetaPC
    Cb.beta<-Cb1/sum(Cb1)
    ####################alphas#########################
    BD<-B(Delta,solve(h0),B2)
    autovmaxdeltaPC<-eigen(BD,symmetric=TRUE)$val[1]
    vetordeltaPC<-eigen(BD,symmetric=TRUE)$vec[,1]
    dmaxG.alpha=abs(vetordeltaPC)
    vCideltaPC=2*abs(diag(BD))
    Cb2=vCideltaPC
    Cb.alpha=Cb2/sum(Cb2)

    result <- list(dmax.beta = dmaxG.beta,
                   dmax.alpha = dmaxG.alpha,
                   dmax.theta = dmaxG.theta,
                   Ci.beta = Cb.beta,
                   Ci.alpha = Cb.alpha,
                   Ci.theta = Cb.theta)
    return(result)
  }

  if(scheme=="response")
  {
    ############################Response####################################
    mu <- model$mu.fv
    sigma <- model$sigma.fv
    eta <- linkfun(mu)
    ai <- mu.eta(eta)
    tau <- sigma_linkfun(sigma)
    bi <- sigma_mu.eta(tau)
    sy<- sqrt((mu*(1+mu))/sigma)
    Phi <-  (1+sigma)

    dymu <- Phi*(1/(y*(1+y)))
    Deltamu <- crossprod(x,diag(ai*dymu*sy))
    p<-ncol(x)
    q<-ncol(z)
    dysigma <- mu*(1/(y*(1+y))) - 1/(1+y)
    Deltasigma <- crossprod(z,diag(bi*dysigma*sy))
    Delta <- rbind(Deltamu,Deltasigma)

    ###############thetas###########################
    BT<-B(Delta,solve(h0),B3)
    autovmaxthetaPC<- eigen(BT,symmetric=TRUE)$val[1]
    vetorthetaRP<- eigen(BT,symmetric=TRUE)$vec[,1]
    dmaxG.theta<-abs(vetorthetaRP)
    vCithetaRP<-2*abs(diag(BT))
    Cb0<-vCithetaRP
    Cb.theta<-Cb0/sum(Cb0)

    #################betas##########################
    BM=B(Delta,solve(h0),B1)
    autovmaxbetaRP <- eigen(BM,symmetric=TRUE)$val[1]
    vetorbetaRP <- eigen(BM,symmetric=TRUE)$vec[,1]
    dmaxG.beta <- abs(vetorbetaRP)
    vCibetaRP <- 2*abs(diag(BM))
    Cb1 <- vCibetaRP
    Cb.beta <- Cb1/sum(Cb1)
    ####################alpha#######################
    BD=B(Delta,solve(h0),B2)
    autovmaxdeltaRP <- eigen(BD,symmetric=TRUE)$val[1]
    vetordeltaRP <- eigen(BD,symmetric=TRUE)$vec[,1]
    dmaxG.alpha <- abs(vetordeltaRP)
    vCideltaRP <- 2*abs(diag(BD))
    Cb2 <- vCideltaRP
    Cb.alpha <- Cb2/sum(Cb2)


    result <- list(dmax.beta = dmaxG.beta,
                   dmax.alpha = dmaxG.alpha,
                   dmax.theta = dmaxG.theta,
                   Ci.beta = Cb.beta,
                   Ci.alpha = Cb.alpha,
                   Ci.theta = Cb.theta)
    return(result)
  }


}

#'@rdname diag.RBS
#'
#'@export
residuals.pearson <- function(model)
{

  mu <- fit$mu.fv
  phi <- fit$sigma.fv
  y <- model$y

  resP <- (sqrt(phi)*(y-mu))/sqrt(mu*(1+mu))


  resP
}


#'@name envelope
#'
#'@aliases envelope.BP
#'@aliases envelope.GA
#'@aliases envelope.IG
#'@aliases envelope.RBS
#'@aliases envelope.WEI3
#'
#'@title Envelopes
#'
#'@description A normal plot with simulated envelope of the residual is produced.
#'
#' @param model object of class \code{gamlss} holding the fitted model.
#' @param k number of replications for envelope construction. Default is 19.
#' @param res type of residuals to be extracted. Default is deviance.
#' @param precision If \code{precision = "fixed"} a model with precision fixed is used;
#' else a model with precision variable is used.
#' @param dist The function RBS() defines the RBS distribution.
#' @param xlabel a label for the x axis.
#' @param color the color of the envelope.
#' @param ylabel a label for the y axis.
#' @param font the font used in x and y axis.
#'
#'
#'@return A simulated envelope of the class RBS, ZARBS and ZAGA.
#'
#' @author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@references
#'Atkinson, A. C. (1985) Plots, transformations and regression : an introduction to graphical methods of diagnostic regression analysis. Oxford Science Publications, Oxford.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'Santos-Neto, M., Cysneiros, F.J.A., Leiva, V., Barros, M. (2016) Reparameterized Birnbaum-Saunders
#'regression models with varying precision. \emph{Electronic Journal of Statistics}, 10, 2825--2855. doi: \email{10.1214/16-EJS1187}.
#'
#'@importFrom graphics par points polygon
#'@importFrom stats qqnorm
#'@importFrom gamlss gamlss gamlss.control
#'@import ggplot2
#'@export

envelope.BP <- function(model,k=100,color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{

  n=model$N
  td  = model$residuals
  re <- matrix(0,n,k)

  for(i in 1:k)
  {
    y1 <-rnorm(n)
    re[,i] <- sort(y1)
  }
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)

  for(l in 1:n)
  {
    eo = sort(re[l,])
    e10[l] = eo[ceiling(k*0.01)]
    e20[l] = eo[ceiling(k*(1 - 0.01))]
    e11[l] = eo[ceiling(k*0.05)]
    e21[l] = eo[ceiling(k*(1 - 0.05))]
    e12[l] = eo[ceiling(k*0.1)]
    e22[l] = eo[ceiling(k*(1 - 0.1))]
  }

  a <- qqnorm(e10, plot.it = FALSE)$x
  r <- qqnorm(td, plot.it = FALSE)$x
  xb = apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x

  df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(text=element_text(size=30,family=font))
}

#'@rdname envelope
#'
#'@export
envelope.GA <- function(model,k=100,color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{

  n=model$N
  td  = model$residuals
  re <- matrix(0,n,k)

  for(i in 1:k)
  {
    y1 <- rnorm(n)
    re[, i] = sort(y1)
  }
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)

  for(l in 1:n)
  {
    eo = sort(re[l,])
    e10[l] = eo[ceiling(k*0.01)]
    e20[l] = eo[ceiling(k*(1 - 0.01))]
    e11[l] = eo[ceiling(k*0.05)]
    e21[l] = eo[ceiling(k*(1 - 0.05))]
    e12[l] = eo[ceiling(k*0.1)]
    e22[l] = eo[ceiling(k*(1 - 0.1))]
  }

  a <- qqnorm(e10, plot.it = FALSE)$x
  r <- qqnorm(td, plot.it = FALSE)$x
  xb = apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x

  df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel) +theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())  +theme(text=element_text(size=30,family=font))
}

#'@rdname envelope
#'
#'@export
envelope.IG <- function(model,k=100,color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{

  n=model$N
  td  = model$residuals
  re <- matrix(0,n,k)

  for(i in 1:k)
  {
    y1 <- rnorm(n)
    re[,i] <- sort(y1)
  }
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)

  for(l in 1:n)
  {
    eo = sort(re[l,])
    e10[l] = eo[ceiling(k*0.01)]
    e20[l] = eo[ceiling(k*(1 - 0.01))]
    e11[l] = eo[ceiling(k*0.05)]
    e21[l] = eo[ceiling(k*(1 - 0.05))]
    e12[l] = eo[ceiling(k*0.1)]
    e22[l] = eo[ceiling(k*(1 - 0.1))]
  }

  a <- qqnorm(e10, plot.it = FALSE)$x
  r <- qqnorm(td, plot.it = FALSE)$x
  xb = apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x

  df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel) +theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank()) +theme(text=element_text(size=30,family=font))
}

#'@rdname envelope
#'
#'@export
envelope.RBS <- function(model,k=100,color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{
  n <- model$N
  td  <-  model$residuals
  re <- matrix(0,n,k)

  for(i in 1:k)
  {
    y1 <- rnorm(n)
    re[,i] <- sort(y1)
  }
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)

  for(l in 1:n)
  {
    eo = sort(re[l,])
    e10[l] = eo[ceiling(k*0.01)]
    e20[l] = eo[ceiling(k*(1 - 0.01))]
    e11[l] = eo[ceiling(k*0.05)]
    e21[l] = eo[ceiling(k*(1 - 0.05))]
    e12[l] = eo[ceiling(k*0.1)]
    e22[l] = eo[ceiling(k*(1 - 0.1))]
  }

  a <- qqnorm(e10, plot.it = FALSE)$x
  r <- qqnorm(td, plot.it = FALSE)$x
  xb = apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x

  df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel) +theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(text=element_text(size=30,family=font))
}

#'@rdname envelope
#'
#'@export
envelope.WEI3 <- function(model,k=100,link=c("log","identity"),color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{

  n=model$N
  td  = model$residuals
  re <- matrix(0,n,k)

  for(i in 1:k)
  {
    y1 <- rnorm(n)
    re[,i] <- sort(y1)
  }
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)

  for(l in 1:n)
  {
    eo = sort(re[l,])
    e10[l] = eo[ceiling(k*0.01)]
    e20[l] = eo[ceiling(k*(1 - 0.01))]
    e11[l] = eo[ceiling(k*0.05)]
    e21[l] = eo[ceiling(k*(1 - 0.05))]
    e12[l] = eo[ceiling(k*0.1)]
    e22[l] = eo[ceiling(k*(1 - 0.1))]
  }

  a <- qqnorm(e10, plot.it = FALSE)$x
  r <- qqnorm(td, plot.it = FALSE)$x
  xb = apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x

  df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel) +theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank()) +theme(text=element_text(size=30,family=font))
}
