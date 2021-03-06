% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SurvEval.R
\name{SurvEval}
\alias{SurvEval}
\title{Evaluate whether a true survival function is contained
in the credible set.}
\usage{
SurvEval(time.grid, true.surv, post.mean, radius)
}
\arguments{
\item{time.grid}{The time grid on which to evaluate the survival
function.}

\item{true.surv}{The true survival function.}

\item{post.mean}{The posterior mean of the survival function, given
as a function.}

\item{radius}{The radius of the credible set for the survival function}
}
\value{
\item{covered}{Indicator whether the true survival function is
completely covered by the credible set on the times contained in
\code{time.grid}. 0 = not completely covered, 1 = completely covered.}
}
\description{
This function is intended to evaluate the Bayesian procedure in a
simulation study. To that end, this function can be used to check
whether the true (user-defined) survival function is
contained in the credible set generated by the function
\link{BayesSurv}.
}
\examples{
#Demonstration on a simulated data set
library(simsurv)
library(ggplot2)
hazard.true <- function(t,x, betas, ...){1.2*(5*(t+0.05)^3 - 10*(t+0.05)^2 + 5*(t+0.05) ) + 0.7}
cumhaz.true <- Vectorize( function(t){integrate(hazard.true, 0, t)$value} )
surv.true <- function(t){exp(-cumhaz.true(t))}

sim.df <- data.frame(id = 1:1000)
df <- simsurv(x = sim.df, maxt = 1, hazard = hazard.true)

bs <- BayesSurv(df, "eventtime", "status")
surv.pm <- approxfun(bs$surv.eval.grid, bs$surv.post.mean)
SurvEval(bs$surv.eval.grid, surv.true, surv.pm, bs$surv.radius)


}
\references{
Castillo and Van der Pas (2020). Multiscale Bayesian survival
  analysis. <arXiv:2005.02889>.
}
\seealso{
\link{BayesSurv}, which computes the posterior mean of the
survival function as well as the radius for its credible set.
}
