get_coefmat <- function(x, j) {
  x[, j, ]
}

array_to_mat <- function(x, model_name, regression_pars, covariate_names, type="beta") {
  model_pars_order <- switch(
    model_name,
    "orl_regression" = c("Arew", "Apun", "K", "betaF", "betaP"),
    "pvl_decay_regression" = c("A", "alpha", "cons", "lambda"),
    "vpp_regression" = c("A", "alpha", "cons", "lambda", "K", "w", "epP", "epN")
  )

  # Re-order user supplied `regression_pars` to match the order used in Stan code
  regression_pars <- regression_pars[order(match(regression_pars, model_pars_order))]

  out <- mapply(get_coefmat, 1:length(regression_pars), MoreArgs=list(x = x), SIMPLIFY=FALSE)
  out <- do.call(cbind, out)

  colnames_out <- as.character(t(outer(regression_pars, covariate_names, paste, sep="_")))
  if (type == "sigma") {
    colnames_out <- paste0("sigma_", colnames_out)
  }
  colnames(out) <- colnames_out

  out
}
