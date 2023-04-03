## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

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

