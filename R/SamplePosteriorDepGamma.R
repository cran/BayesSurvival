#' Draw samples from the posterior for the hazard, using the piecewise
#' exponential (histogram) prior with dependent Gamma heights
#'
#' The sampler is described in the Supplement to Castillo and Van der Pas
#' (2020) and uses MCMC within Gibbs, with a Gamma proposal with shape
#' parameter equal to the number of events in each interval plus some epsilon
#' (to prevent proposals equal to zero if there are no events in an interval)
#' and rate parameter equal to the parameter alpha (set by the user) divided by
#' histogram height on the previous interval, plus the total amount of time all
#' individuals were exposed during this interval. Most users of the package will
#' not work with this function directly, but instead use the main function
#' \link{BayesSurv}, in which this particular function is incorporated.
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
#' @param alpha.dep The main parameter \eqn{\alpha} for the dependent Gamma
#'   prior, as described in the documentation for \link{BayesSurv}. It is
#'   recommended to take \code{alpha.dep} smaller than \code{alpha0.dep}.
#' @param alpha0.dep The shape parameter for the Gamma prior on the histogram
#'   height for the first interval. It is recommended to take \code{alpha.dep}
#'   smaller than \code{alpha0.dep}.
#' @param beta0.dep The rate parameter for the Gamma prior on the histogram
#'   height for the first interval.
#'
#' @return \item{samples}{A \eqn{N} by \eqn{K} (the number of draws by the
#' number of intervals) matrix, with each row containing a draw from the
#' posterior for the hazard, based on a histogram prior with dependent
#' Gamma heights.}
#'
#' @export

SamplePosteriorDepGamma <- function(failures, exposures, N = 1000, alpha.dep = 1, alpha0.dep = 1.5, beta0.dep = 1){

  K <- length(failures)

  epsilon <- 0.001 #to prevent an initial 'current' estimate of zero
  samples <- matrix(0, nrow = N, ncol = K)

  samples[1, ] <- failures/(exposures + 1) + epsilon #initialization

  if(K == 1){
    samples[ , 1] <- rgamma(N, shape = (failures[1] + alpha0.dep), rate = (exposures[1] + beta0.dep))
  } #end if(K == 1)


  if(K == 2){
    for(i in 2:N){
      samples[i, 1] <- MCMCDepGammaFirst(samples[(i-1), 1], samples[(i-1), 2], failures[1], exposures[1], alpha.dep, alpha0.dep, beta0.dep)

      samples[i, K] <- rgamma(1, failures[K] + alpha.dep, alpha.dep/samples[i, (K-1)] + exposures[K]) #final interval, direct sampling possible
    }
  } #end if(K == 2)

 if(K > 2){
   for(i in 2:N){
     samples[i, 1] <- MCMCDepGammaFirst(samples[(i-1), 1], samples[(i-1), 2], failures[1], exposures[1], alpha.dep, alpha0.dep, beta0.dep)

     for(k in 2:(K-1)){
       samples[i, k] <- MCMCDepGammaIntermediate(samples[(i-1), k], samples[i, (k-1)], samples[(i-1), (k+1)], failures[k], exposures[k], alpha.dep)
     }

     samples[i, K] <- rgamma(1, failures[K] + alpha.dep, alpha.dep/samples[i, (K-1)] + exposures[K]) #final interval, direct sampling possible
   }
 } #end if(K > 2)

  return(samples)

}#end sample.posterior.dep.gamma
