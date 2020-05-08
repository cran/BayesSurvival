
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesSurvival

The goal of BayesSurvival is to perform unadjusted Bayesian survival
analysis for right censored time-to-event data. The main function
(BayesSurv) computes the posterior mean and a credible band for the
survival function and for the cumulative hazard, as well as the
posterior mean for the hazard, starting from a piecewise exponential
(histogram) prior with Gamma distributed heights that are either
independent, or have a Markovian dependence structure. A function
(PlotBayesSurv) is provided to easily create plots of the posterior
means of the hazard, cumulative hazard and survival function, with a
credible band accompanying the latter two. The priors and samplers are
described in more detail in the preprint ‘Multiscale Bayesian survival
analysis’ by Castillo and Van der Pas (2020+). In that paper it is also
shown that the credible bands for the survival function and the
cumulative hazard can be considered confidence bands (under mild
conditions) and thus offer reliable uncertainty quantification.

## Installation

You can install the released version of BayesSurvival from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("BayesSurvival")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BayesSurvival)
library(simsurv)
hazard.true <- function(t,x, betas, ...){1.2*(5*(t+0.05)^3 - 10*(t+0.05)^2 + 5*(t+0.05) ) + 0.7}
sim.df <- data.frame(id = 1:1000)
df <- simsurv(x = sim.df, maxt = 1, hazard = hazard.true)

bs <- BayesSurv(df, "eventtime", "status")
```
