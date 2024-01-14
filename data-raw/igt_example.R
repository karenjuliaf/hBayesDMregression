library(tidyverse)

set.seed(25)

igt_example <- system.file("extdata/igt_exampleData.txt", package="hBayesDM") %>%
  read_delim()

igt_example <- igt_example %>%
  group_nest(subjID) %>%
  mutate(
    age = sample(25:65, size=n(), replace=TRUE),
    x1 = runif(n(), min=4, max=6),
    x2 = runif(n(), min=3500, max=4500),
    x3 = runif(n(), min=4, max=6),
    x4 = runif(n(), min=1e4, max=1.5e4),
    sex = factor(sample(c(0, 1), size=n(), replace=TRUE)),
    cond = sample(1:3, size=n(), replace=TRUE)
  ) %>%
  mutate(across(x1:x4, ~(.x - mean(.x))/sd(.x))) %>%
  pivot_wider(
    names_from = cond, values_from = cond,
    names_prefix = "cond", names_sep = "", names_sort = TRUE,
    values_fill = 0, values_fn = ~if_else(.x > 0, 1, 0)
  ) %>%
  unnest(data) %>%
  as.data.frame()

usethis::use_data(igt_example, overwrite=TRUE)
