\name{shash}
\alias{shash}
\title{Sinh-arcsinh model family}
\usage{
shash(link = list("identity", "logeb", "identity", "slogit"), 
      b = 0.01, a1 = -10, a2 = 1)
}
\arguments{
\item{link}{vector of four characters indicating the link function for location, scale, skewness and kurtosis.}
\item{b}{positive parameter of the logeb link function, see Details.}
\item{a1}{parameter of the slogit link function, see Details.}
\item{a2}{parameter of the slogit link function, see Details.}
}
\value{
An object inheriting from class general.family.
}
\description{
The \code{shash} family implements the four-parameter sinh-arcsinh (shash) distribution of 
             Jones and Pewsey (2009). The location, scale, skewness and kurtosis of the density can depend 
             on additive smooth predictors. Useable only with gam, the linear predictors are specified 
             via a list of formulae.
}
\details{
The density function of the shash family is 
         \deqn{p(y|\mu,\sigma,\epsilon,\delta)=C(z) exp{-S(z)^2/2} / \sigma{2\pi(1+z^2)}^1/2, }
         where \eqn{ C(z)={1+S(z)^2}^1/2 }, \eqn{ S(z)=sinh{\delta sinh^(-1)(z)-\epsilon} } and
         \eqn{z=(y-\mu)/(\sigma\delta)}. Here \eqn{\mu} and \eqn{\sigma > 0} control, respectively, location and 
         scale, \eqn{\epsilon} determines skewness, while \eqn{\delta > 0} controls tailweight.
         \code{shash} can model skewness to either side, depending on the sign of \eqn{\epsilon}. 
         Also, shash can have tails that are lighter (\eqn{\delta>1}) or heavier (\eqn{0<\delta<1}) that a normal.
         For fitting purposes, here we are using \eqn{\tau = log(\sigma)} and \eqn{\phi = log(\delta)}.
         
The link function used for \eqn{\tau} is logeb with is \eqn{\eta = log{exp(\tau)-b}} so that the inverse link is
\eqn{\tau = log(\sigma) = log{exp(\eta)+b}}. The point is that we are don't allow \eqn{\sigma} to become smaller
than a small constant b. The link function used for \eqn{\phi} is \eqn{\eta = logit(c)} with \eqn{c = (\phi-a1)/(a2-a1)}, 
so that the inverse link is \eqn{\phi = log(\delta) = (a2-a1) * logistic(\eta) + a1}. The purpose is keeping 
\eqn{\delta} between \eqn{exp(a1)} and \eqn{exp(a2)} as the shape of the density becomes insensitive to \eqn{delta} if
this parameter takes extreme values.
}
\examples{

#########
# Shash dataset
#########
##  Simulate some data form shash
set.seed(847)
n <- 1000
x <- seq(-4, 4, length.out = n)

X <- cbind(1, x, x^2)
beta <- c(4, 1, 1)
mu <- X \%*\% beta 

sigma =  .5+0.4*(x+4)*.5            # Scale
eps = 2*sin(x)                      # Skewness
del = 1 + 0.2*cos(3*x)              # Kurtosis

dat <-  mu + (del * sigma) * sinh((1/del) * asinh(qnorm(runif(n))) + (eps/del))
dataf <- data.frame(cbind(dat, x))
names(dataf) <- c("y", "x")
plot(x, dat, xlab = "x", ylab = "y")

## Fit model
fit <- gam(list(y ~ s(x), # <- model for location 
                  ~ s(x),   # <- model for log-scale
                  ~ s(x),   # <- model for skewness
                  ~ s(x, k = 20)), # <- model for log-kurtosis
           data = dataf, 
           family = shash, # <- new family 
           optimizer = "efs") # Here we are using a new optimizer 

## Plotting truth and estimates for each parameters of the density 
muE <- fit$fitted[ , 1]
sigE <- exp(fit$fitted[ , 2])
epsE <- fit$fitted[ , 3]
delE <- exp(fit$fitted[ , 4])

par(mfrow = c(2, 2))
plot(x, muE, type = 'l', ylab = expression(mu(x)), lwd = 2)
lines(x, mu, col = 2, lty = 2, lwd = 2)
legend("top", c("estimated", "truth"), col = 1:2, lty = 1:2, lwd = 2)

plot(x, sigE, type = 'l', ylab = expression(sigma(x)), lwd = 2)
lines(x, sigma, col = 2, lty = 2, lwd = 2)

plot(x, epsE, type = 'l', ylab = expression(epsilon(x)), lwd = 2)
lines(x, eps, col = 2, lty = 2, lwd = 2)

plot(x, delE, type = 'l', ylab = expression(delta(x)), lwd = 2)
lines(x, del, col = 2, lty = 2, lwd = 2)

## Plotting true and estimated conditional density
par(mfrow = c(1, 1))
plot(x, dat, pch = '.', col = "grey", ylab = "y", ylim = c(-35, 70))
for(qq in c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)){
  est <- fit$family$qf(p=qq, mu = fit$fitted)
  true <- mu + (del * sigma) * sinh((1/del) * asinh(qnorm(qq)) + (eps/del))
  lines(x, est, type = 'l', col = 1, lwd = 2)
  lines(x, true, type = 'l', col = 2, lwd = 2, lty = 2)
}
legend("topleft", c("estimated", "truth"), col = 1:2, lty = 1:2, lwd = 2)
}
\author{
Matteo Fasiolo <matteo.fasiolo@gmail.com> and Simon N. Wood.
}
\references{
Jones, M. and A. Pewsey (2009). Sinh-arcsinh distributions. Biometrika 96 (4), 761–780.

            Wood, Simon N., Pya, N. and Safken, B. (2017). Smoothing parameter and model selection for 
            general smooth models. Journal of the American Statistical Association.
}

