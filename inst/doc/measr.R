## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install-rstan, eval = FALSE----------------------------------------------
#  install.packages("rstan")

## ----rstan-test, eval = FALSE-------------------------------------------------
#  library(rstan)
#  
#  example(stan_model, package = "rstan", run.dontrun = TRUE)

## ----stan-repo, eval = FALSE--------------------------------------------------
#  install.packages("cmdstanr",
#                   repos = c("https://mc-stan.org/r-packages/",
#                             getOption("repos")))

## ----stan-dev, eval = FALSE---------------------------------------------------
#  # install.packages("remotes")
#  remotes::install_github("stan-dev/cmdstanr")

## ----check-toolchain, eval = FALSE--------------------------------------------
#  library(cmdstanr)
#  
#  check_cmdstan_toolchain()

## ----install-cmdstan, eval = FALSE--------------------------------------------
#  install_cmdstan(cores = 2)

## ----install-measr, eval = FALSE----------------------------------------------
#  install.packages("measr")

## ----measr-dev, eval = FALSE--------------------------------------------------
#  # install.packages("remotes")
#  remotes::install_github("wjakethompson/measr")

## ----load-pkg-----------------------------------------------------------------
library(measr)

## ----data---------------------------------------------------------------------
ecpe_data

ecpe_qmatrix

## ----est-hide, include = FALSE------------------------------------------------
ecpe_lcdm <- measr_dcm(data = ecpe_data, qmatrix = ecpe_qmatrix,
                       resp_id = "resp_id", item_id = "item_id",
                       method = "optim", type = "lcdm",
                       file = "fits/ecpe-optim-lcdm")

## ----est-show, eval = FALSE---------------------------------------------------
#  ecpe_lcdm <- measr_dcm(data = ecpe_data, qmatrix = ecpe_qmatrix,
#                         resp_id = "resp_id", item_id = "item_id",
#                         method = "optim", type = "lcdm")

## ----resp-prob, message = FALSE, warning = FALSE, error = FALSE---------------
ecpe_lcdm <- add_respondent_estimates(ecpe_lcdm)
measr_extract(ecpe_lcdm, "attribute_prob")

## -----------------------------------------------------------------------------
ecpe_lcdm <- add_reliability(ecpe_lcdm)
measr_extract(ecpe_lcdm, "classification_reliability")

