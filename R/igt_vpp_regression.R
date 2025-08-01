#' @templateVar MODEL_FUNCTION igt_vpp_regression
#' @templateVar CONTRIBUTOR NA
#' @templateVar TASK_NAME Iowa Gambling Task
#' @templateVar TASK_CODE igt
#' @templateVar TASK_CITE (Ahn et al., 2008)
#' @templateVar MODEL_NAME Value-Plus-Perseverance
#' @templateVar MODEL_CODE vpp
#' @templateVar MODEL_CITE (Worthy et al., 2013)
#' @templateVar MODEL_TYPE Hierarchical
#' @templateVar DATA_COLUMNS "subjID", "choice", "gain", "loss"
#' @templateVar PARAMETERS \code{A} (learning rate), \code{alpha} (outcome sensitivity), \code{cons} (response consistency), \code{lambda} (loss aversion), \code{epP} (gain impact), \code{epN} (loss impact), \code{K} (decay rate), \code{w} (RL weight)
#' @templateVar REGRESSORS NA
#' @templateVar POSTPREDS "y_pred"
#' @templateVar LENGTH_DATA_COLUMNS 4
#' @templateVar DETAILS_DATA_1 \item{subjID}{A unique identifier for each subject in the data-set.}
#' @templateVar DETAILS_DATA_2 \item{choice}{Integer indicating which deck was chosen on that trial (where A==1, B==2, C==3, and D==4).}
#' @templateVar DETAILS_DATA_3 \item{gain}{Floating point value representing the amount of currency won on that trial (e.g. 50, 100).}
#' @templateVar DETAILS_DATA_4 \item{loss}{Floating point value representing the amount of currency lost on that trial (e.g. 0, -50).}
#' @templateVar LENGTH_ADDITIONAL_ARGS 1
#' @templateVar ADDITIONAL_ARGS_1 \item{payscale}{Raw payoffs within data are divided by this number. Used for scaling data. Defaults to 100.}
#'
#' @template model-documentation
#'
#' @export
#' @include hBayesDMregression_model.R
#' @include preprocess_funcs.R
#'
#' @references
#' Ahn, W. Y., Busemeyer, J. R., & Wagenmakers, E. J. (2008). Comparison of decision learning models using the generalization criterion method. Cognitive Science, 32(8), 1376-1402. https://doi.org/10.1080/03640210802352992
#'
#' Worthy, D. A., & Todd Maddox, W. (2013). A comparison model of reinforcement-learning and win-stay-lose-shift decision-making processes: A tribute to W.K. Estes. Journal of Mathematical Psychology, 59, 41-49. https://doi.org/10.1016/j.jmp.2013.10.001
#'
#' @examples
#' \dontrun{
#' # Run the model with the `igt_example` data set
#' data("igt_example")
#'
#' output <- igt_vpp_regression(
#'   data = igt_example,
#'   exclude_cols = "trial",
#'   regression_pars = c("lambda", "w"),
#'   niter = 2000,
#'   nwarmup = 1000,
#'   nchain = 4,
#'   ncore = 1
#' )
#'
#' # Preview posterior samples of the regression coefficients; a matrix with dimension
#' # (niter - nwarmup) x (# of covariates)
#' head(output$parVals$beta)
#'
#' # Preview posterior samples of the sigmas corresponding to the regression coefficients; same
#' # dimension as previous matrix
#' head(output$parVals$sigma_beta)
#'
#' # View Stan code used to fit the model
#' cat(output$model_code)
#'
#' # For visual diagnostics, see the "Getting started" vignette:
#' # vignette("getting-started", package="hBayesDMregression")
#' }

igt_vpp_regression <- hBayesDMregression_model(
  task_name       = "igt",
  model_name      = "vpp_regression",
  model_type      = "",
  data_columns    = c("subjID", "choice", "gain", "loss"),
  parameters      = list(
    "A" = c(0, 0.5, 1),
    "alpha" = c(0, 0.5, 2),
    "cons" = c(0, 1, 5),
    "lambda" = c(0, 1, 10),
    "epP" = c(-Inf, 0, Inf),
    "epN" = c(-Inf, 0, Inf),
    "K" = c(0, 0.5, 1),
    "w" = c(0, 0.5, 1)
  ),
  regressors      = NULL,
  postpreds       = c("y_pred"),
  preprocess_func = igt_regression_preprocess_func)
