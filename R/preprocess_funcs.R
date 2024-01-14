#' @noRd

igt_regression_preprocess_func <- function(raw_data, general_info, regression_pars, payscale = 100) {
  # Currently class(raw_data) == "data.table"

  # Use general_info of raw_data
  subjs   <- general_info$subjs
  n_subj  <- general_info$n_subj
  t_subjs <- general_info$t_subjs
  t_max   <- general_info$t_max

  # Initialize data arrays
  Ydata    <- array(-1, c(n_subj, t_max))
  RLmatrix <- array( 0, c(n_subj, t_max))

  # Write from raw_data to the data arrays
  for (i in 1:n_subj) {
    subj <- subjs[i]
    t <- t_subjs[i]
    DT_subj <- raw_data[raw_data$subjid == subj]

    Ydata[i, 1:t]    <- DT_subj$choice
    RLmatrix[i, 1:t] <- DT_subj$gain - abs(DT_subj$loss)
  }

  # Check that there is at least one covariate supplied
  if (ncol(raw_data) == 4) {
    stop("At least one covariate must be supplied.")
  }

  npars <- length(regression_pars)

  # Remove insensitive data columns to obtain covariate matrix
  insensitive_data_columns <- c("choice", "gain", "loss")
  covariate_matrix <- raw_data[, -insensitive_data_columns, with=FALSE]

  # Get distinct rows by subject id
  covariate_matrix <- unique(covariate_matrix, by="subjid")
  covariate_matrix <- covariate_matrix[, -"subjid"]

  # To avoid NOTEs by R CMD check
  .SD <- NULL

  # Determine which covariates are factors by assuming that any covariates that have two or fewer
  # distinct values over the entire data set are factors
  are_factors <- covariate_matrix[, lapply(.SD, as.factor)]
  are_factors <- are_factors[, lapply(.SD, nlevels)]
  are_factors <- are_factors[, lapply(.SD, function(x) {x <= 2})]
  are_factors <- stats::setNames(as.logical(are_factors), names(are_factors))

  # For non-factor covariates, we can calculate their standard deviations as usual
  # For factor covariates, we set their standard deviations to 1
  if (all(!are_factors)) {
    covariate_sds <- covariate_matrix[, lapply(.SD, stats::sd)]
    covariate_sds <- stats::setNames(as.numeric(covariate_sds), names(covariate_sds))
  } else {
    non_factor_sds <- covariate_matrix[, lapply(.SD, stats::sd), .SDcols=!are_factors]
    non_factor_sds <- as.numeric(non_factor_sds)

    covariate_sds <- numeric(ncol(covariate_matrix))
    covariate_sds[are_factors] <- 1
    covariate_sds[covariate_sds != 1] <- non_factor_sds
    covariate_sds <- stats::setNames(covariate_sds, colnames(covariate_matrix))
  }

  # Before converting to matrix, convert all factors to numeric, as per details in
  # `data.table:::as.matrix.data.table`, otherise a character matrix will be returned
  true_factor_columns <- sapply(covariate_matrix, is.factor)
  true_factor_columns <- names(true_factor_columns)[true_factor_columns]
  if (length(true_factor_columns > 1)) {
    covariate_matrix <- covariate_matrix[, lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols=true_factor_columns]
  }

  covariate_matrix <- as.matrix(covariate_matrix)
  covariate_precisions <- 1 / covariate_sds
  covariate_precisions <- matrix(
    covariate_precisions,
    nrow=npars, ncol=length(covariate_precisions), byrow=TRUE
  )

  # Wrap into a list for Stan
  data_list <- list(
    N        = n_subj,
    T        = t_max,
    Tsubj    = t_subjs,
    choice   = Ydata,
    outcome  = RLmatrix / payscale,
    sign_out = sign(RLmatrix),
    ncov = ncol(covariate_matrix),
    covariate_precisions = covariate_precisions,
    covariate_matrix = covariate_matrix
  )

  # Returned data_list will directly be passed to Stan
  return(data_list)
}
