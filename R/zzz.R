#' @noRd

.onLoad <- function(libname, pkgname) {
  # nocov start
  if (FLAG_BUILD_ALL) {
    modules <- paste0("stan_fit4", names(stanmodels), "_mod")
    for (m in modules) loadModule(m, what = TRUE)
  }
} # nocov end
