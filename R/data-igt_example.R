#' IGT example data
#'
#' This is an example data set to be used with the `igt_*_regression()` functions. It is an
#' augmented version of the example IGT data set found in the `hBayesDM` package. Each subject in
#' this data set performs 100 trials. As such, the values found in the columns containing the
#' augmented covariates will be identical for all trials belonging to the same subject.
#'
#' @format A data frame with 400 observations and 14 variables:
#' \describe{
#' \item{subjID}{Subject ID}
#' \item{trial}{Trial number}
#' \item{choice}{The number of the deck chosen}
#' \item{gain,loss}{The gain(s) and/or loss(es) associated with the choice}
#' \item{age}{Age of the subject}
#' \item{x1,x2,x3,x4}{Continuous covariates}
#' \item{sex}{Binary indicator of the subject's sex; 0 = female, 1 = male}
#' \item{cond1,cond2,cond3}{One-hot encoded categorical covariate indicating presence of one of
#'   three conditions}
#' }
#'
#' @seealso The
#'   \href{https://github.com/adamoshen/hBayesDMregression/data-raw/igt_example.R}{script} used to
#'   create the igt_example data set.
#'
"igt_example"
