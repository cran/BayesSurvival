#' Draw samples from the posterior for the hazard, using the piecewise
#' exponential (histogram) prior with independent Gamma heights
#'
#' The sampler is described in the Supplement to Castillo and Van der Pas
#' (2020). Most users of the package will not work with this function directly,
#' but instead use the main function \link{BayesSurv}, in which this particular
#' function is incorporated.
#'
#' The samples returned by this function are draws from the posterior for the
#' hazard function. To obtain draws from the posterior for the cumulative
#' hazard, one can use numerical integration. One way to achieve this is by
#' first finding the values of the cumulative hazard at the end of each interval,
#' e.g. by \code{t(apply(samples*time.max/K, 1, cumsum))}, where \code{samples}
#' is the output from the present function and \code{time.max} and \code{K} are
#' as described for \link{BayesSurv}, and then using \code{approxfun()} to linearly
#' interpolate in between. To obtain posterior samples from the survival, one
#' could then use \link{SurvivalFromCumhaz}.
#'
#' @seealso \link{BayesSurv}, which computes the posterior mean and credible
#'   bands for the cumulative hazard and survival functions, as well as the
#'   posterior mean for the hazard. Within \link{BayesSurv}, the present
#'   function is called.
#'
#' @references Castillo and Van der Pas (2020). Multiscale Bayesian survival
#'   analysis. <arXiv:2005.02889>.
#'
#' @param failures A vector of length \eqn{K} (the total number of intervals),
#'   containing for each interval the number of individuals who had an event
#'   during that interval.
#' @param exposures A vector of length \eqn{K} (the total number of intervals),
#'   containing for each interval the total amount of time all individuals
#'   together were under follow-up during that interval.
#' @param N The number of draws to take.
#' @param alpha.indep The shape parameter for the Gamma prior on the histogram
#'   height for each interval.
#' @param beta.indep The rate parameter for the Gamma prior on the histogram
#'   height for each interval.
#'
#' @return \item{samples}{A \eqn{N} by \eqn{K} (the number of draws by the
#' number of intervals) matrix, with each row containing a draw from the
#' posterior for the hazard, based on a histogram prior with independent
#' Gamma heights.}
#'
#' @export


SamplePosteriorIndepGamma <- function(failures, exposures, N = 1000, alpha.indep = 1.5, beta.indep = 1){

  K <- length(failures)

  samples <- matrix(0, nrow = N, ncol = K)

  for(k in 1:K){
    samples[ , k] <- rgamma(N, shape = (failures[k] + alpha.indep), rate = (exposures[k] + beta.indep))
  }

  return(samples)

}
