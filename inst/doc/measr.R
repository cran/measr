## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install-rstan, eval = FALSE----------------------------------------------
# install.packages("rstan")

## ----rstan-test, eval = FALSE-------------------------------------------------
# library(rstan)
# 
# example(stan_model, package = "rstan", run.dontrun = TRUE)

## ----stan-repo, eval = FALSE--------------------------------------------------
# install.packages(
#   "cmdstanr",
#   repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
# )

## ----stan-dev, eval = FALSE---------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("stan-dev/cmdstanr")

## ----check-toolchain, eval = FALSE--------------------------------------------
# library(cmdstanr)
# 
# check_cmdstan_toolchain()

## ----install-cmdstan, eval = FALSE--------------------------------------------
# install_cmdstan(cores = 2)

## ----install-measr, eval = FALSE----------------------------------------------
# install.packages("measr")

## ----measr-dev, eval = FALSE--------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("r-dcm/measr")

## ----load-pkg-----------------------------------------------------------------
library(measr)

## ----data---------------------------------------------------------------------
library(dcmdata)

ecpe_data

ecpe_qmatrix

## ----est-hide, include = FALSE------------------------------------------------
ecpe_spec <- dcm_specify(
  ecpe_qmatrix,
  identifier = "item_id",
  measurement_model = lcdm(),
  structural_model = unconstrained()
)

ecpe_lcdm <- dcm_estimate(
  ecpe_spec,
  data = ecpe_data,
  identifier = "resp_id",
  method = "optim",
  backend = "rstan",
  file = "fits/ecpe-optim-lcdm"
)

## ----est-show, eval = FALSE---------------------------------------------------
# ecpe_spec <- dcm_specify(
#   ecpe_qmatrix,
#   identifier = "item_id",
#   measurement_model = lcdm(),
#   structural_model = unconstrained()
# )
# 
# ecpe_lcdm <- dcm_estimate(
#   ecpe_spec,
#   data = ecpe_data,
#   identifier = "resp_id",
#   method = "optim",
#   backend = "rstan"
# )

## ----resp-prob, message = FALSE, warning = FALSE, error = FALSE---------------
ecpe_lcdm <- add_respondent_estimates(ecpe_lcdm)
measr_extract(ecpe_lcdm, "attribute_prob")

## -----------------------------------------------------------------------------
ecpe_lcdm <- add_reliability(ecpe_lcdm)
measr_extract(ecpe_lcdm, "classification_reliability")

