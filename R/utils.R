reorder_regression_pars <- function(model_name, regression_pars) {
  model_pars_order <- switch(
    model_name,
    "orl_regression" = c("Arew", "Apun", "K", "betaF", "betaP"),
    "pvl_decay_regression" = c("A", "alpha", "cons", "lambda"),
    "vpp_regression" = c("A", "alpha", "cons", "lambda", "K", "w", "epP", "epN")
  )

  # Re-order user supplied `regression_pars` to match the order used in Stan code
  regression_pars <- regression_pars[order(match(regression_pars, model_pars_order))]

  regression_pars
}

create_names <- function(regression_pars, covariate_names, sigma) {
  names_out <- as.character(t(outer(regression_pars, covariate_names, paste, sep="_")))

  if (sigma) {
    names_out <- paste0("sigma_", names_out)
  }

  names_out
}

names_key <- function(regression_pars, covariate_names, sigma) {
  index_matrix <- outer(1:length(regression_pars), 1:length(covariate_names), paste, sep=",")
  index_matrix <- as.character(t(index_matrix))

  new_names <- as.list(create_names(regression_pars, covariate_names, sigma))

  if (sigma) {
    names(new_names) <- paste0("sigma_beta[", index_matrix, "]")
  } else {
    names(new_names) <- paste0("beta[", index_matrix, "]")
  }

  new_names
}

rename_by_key <- function(x, key) {
  ifelse(names(x) %in% names(key), key[names(x)], names(x))
}

get_coefmat <- function(x, j) {
  x[, j, ]
}

array_to_mat <- function(x, model_name, regression_pars, covariate_names, sigma) {
  out <- mapply(get_coefmat, 1:length(regression_pars), MoreArgs=list(x = x), SIMPLIFY=FALSE)
  out <- do.call(cbind, out)

  colnames_out <- create_names(regression_pars, covariate_names, sigma)
  colnames(out) <- colnames_out

  out
}
