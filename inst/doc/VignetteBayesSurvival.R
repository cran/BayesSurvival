## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning = FALSE---------------------------------------------------
library(BayesSurvival)
library(survival)
library(ggplot2)

## -----------------------------------------------------------------------------
cancer$status[cancer$status == 1] <- 0 #censored
cancer$status[cancer$status == 2] <- 1 #died

table(cancer$status)

## -----------------------------------------------------------------------------
res <- BayesSurv(df = cancer, #our data frame
          time = "time", #name of column with survival/censoring times
          event = "status", #name of column with status indicator
          prior = "Dependent", #use dependent Gamma prior
          )

## ---- warning = FALSE, fig.width = 5------------------------------------------
PlotBayesSurv(bayes.surv.object = res,
              object = "survival")

## ---- warning = FALSE, fig.width = 5------------------------------------------
PlotBayesSurv(bayes.surv.object = res,
              object = "survival",
              color = "orange",
              plot.title = "Survival function")

## ---- warning = FALSE, fig.width = 5------------------------------------------
gg <- PlotBayesSurv(bayes.surv.object = res,
              object = "survival")

km <- survfit( Surv(time, status) ~ 1, data = cancer ) #Kaplan-Meier

df.km <- data.frame(t = km$time, km = km$surv)

gg <- gg + geom_line(data = df.km, aes(x = t, y = km), colour = "black", size = 1, lty = 6) 
gg <- gg + labs(title = "With Kaplan-Meier + CI's")
gg

## ---- warning = FALSE, fig.width = 5------------------------------------------
PlotBayesSurv(bayes.surv.object = res,
              object = "cumhaz",
              plot.title = "Cumulative hazard")

## ---- warning = FALSE, fig.width = 5------------------------------------------
PlotBayesSurv(bayes.surv.object = res,
              object = "hazard",
              plot.title = "Hazard")

## -----------------------------------------------------------------------------
res <- BayesSurv(df = cancer, #our data frame
          time = "time", #name of column with survival/censoring times
          event = "status", #name of column with status indicator
          prior = "Independent", #use independent Gamma prior
          )

## ---- warning = FALSE, fig.width = 5------------------------------------------
PlotBayesSurv(bayes.surv.object = res,
              object = "survival",
              plot.title = "Survival")

PlotBayesSurv(bayes.surv.object = res,
              object = "cumhaz",
              plot.title = "Cumulative hazard")

PlotBayesSurv(bayes.surv.object = res,
              object = "hazard",
              plot.title = "Hazard")

