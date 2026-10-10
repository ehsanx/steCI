#' steCI: Confidence Intervals for Sequential Trial Emulation with Additive Hazards
#'
#' Estimation of causal marginal risk differences and constant hazard
#' differences via sequential trial emulation with weighted Aalen additive
#' hazards models, and a family of confidence-interval procedures
#' (model-based, robust, cluster-robust, nonparametric bootstrap, jackknife,
#' and a linearised estimating-function (LEF) / multiplier bootstrap).
#'
#' Reimplements, with citation, the estimator of Keogh et al. (2023,
#' Statistics in Medicine) and the CI-comparison framework of Limozin, Seaman
#' & Su (2025, Statistical Methods in Medical Research).
#' This implementation requires absorbing observed treatment (once treated,
#' always treated). Artificial-censoring weights are estimated in the control
#' arm to account for treatment initiation; treated-arm weights are one.
#'
#' @keywords internal
#' @importFrom stats aggregate ave coef delete.response glm model.matrix
#'   predict quantile qnorm reshape rnorm rpois runif sd stepfun terms
#' @importFrom dplyr %>% group_by mutate ungroup filter first lag lead
#'   row_number bind_rows
#' @importFrom rlang .data
#' @importFrom tidyr replace_na
#' @importFrom timereg aalen const
#' @importFrom survival Surv survfit
#' @importFrom utils globalVariables
"_PACKAGE"
