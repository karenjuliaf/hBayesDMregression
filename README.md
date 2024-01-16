
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hBayesDMregression

hBayesDMregression is an extension of the
[hBayesDM](https://github.com/CCS-Lab/hBayesDM) package, for performing
hierarchical linear regression of Iowa Gambling Task (IGT) model
parameters on a set of user-supplied covariates.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("adamoshen/hBayesDMregression")
```

A list of included vignettes can be found
[below](https://github.com/adamoshen/hBayesDMregression?tab=readme-ov-file#vignette-index).
If you would like a local copy of the vignettes with your installation,
install the package instead using

``` r
remotes::install_github("adamoshen/hBayesDMregression", build_vignettes=TRUE)
```

## Included objects

### Data

The `igt_example` data set can be attached by calling

``` r
data("igt_example", package="hBayesDMregression")
```

This data set is an augmented version of the IGT example data set found
in the [hBayesDM](https://github.com/CCS-Lab/hBayesDM) package,
containing additional continuous and categorical covariates.

### Model

An example model built using the `igt_example` data set and
`hBayesDMregression::igt_orl_regression()` can be found under
`vignettes/orl_model.rds`. This model can be attached by first saving it
to your computer and calling

``` r
orl_model <- readRDS("orl_model.rds")
```

## Vignette index

[Introduction](https://adamoshen.github.io/hBayesDMregression/introduction.html):
provides the mathematical background for the core functions found in
this package.

[Getting
started](https://adamoshen.github.io/hBayesDMregression/getting-started.html):
a quick demonstration on making diagnostic plots using a basic model.
