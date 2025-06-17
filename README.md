# glmmTMBaug

`glmmTMBaug` is an R package that extends the functionality of [`glmmTMB`](https://cran.r-project.org/package=glmmTMB) by providing tools for penalized estimation of random effects and constructing pseudo-data for custom penalization approaches. It enables flexible model fitting while preserving compatibility with the `glmmTMB` interface.

## Features

- Penalized likelihood estimation of random effects for generalizef mixed effects models (limited to gaussian, binomial and Poisson models)
- Seamless integration with `glmmTMB` formula syntax
- Focused penalization of specific levels of random effects
