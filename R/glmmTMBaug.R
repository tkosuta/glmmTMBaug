#' Fit Penalized Generalized Linear Mixed Model via Data Augmentation
#'
#' This function extends `glmmTMB` to fit generalized linear mixed models (GLMMs)
#' with a penalty on the random effects covariance matrix via a data augmentation approach.
#'
#' @param formula A model formula specifying fixed and random effects, as in \code{glmmTMB}.
#' The current implementation supports penalization of a single grouping factor
#' (i.e., one level of random effects).
#' By default, the first random effect term specified in the formula is penalized.
#' @param data A data frame containing the variables used in the model.
#' @param family a family function, a character string naming a family function, or the result of a call to a family function (variance/link function) information. Only \code{binomial}, \code{poisson} and \code{gaussian} are currently supported.
#' @param penOpt A named list of penalty options used to control the penalized likelihood fit. If \code{psi=NULL} and \code{nu=NULL}, a data-driven approach is used to determine psi. If \code{tau=NULL} a data driven approach is used to determine tau.
#' If both \code{psi} and \code{nu} are specified, a prior with those parameters is used instead. In this case, the parameters \code{tau}, \code{trunc} and \code{alpha} are ignored.
#'   The following penalty options are recognized (with defaults):
#'
#'   \itemize{
#'     \item \code{tau} (default: \code{NULL}) — Numeric value in the interval \[0, 1\]. If \code{NULL}, a data-driven method is used to estimate it.
#'     \item \code{trunc} (default: \code{c(1e-4, 1e4)}) — Lower and upper truncation bounds of eigenvalues in data-drive approach to determine \code{psi}.
#'     \item \code{alpha} (default: \code{0.05}) — Significance level used in the data-driven procedure for estimating \code{tau}.
#'     \item \code{psi} (default: \code{NULL}) — A positive-definite matrix specifying the prior scale matrix for the inverse Wishart (variance) or Wishart (precision) distribution. Required if using a fixed predefined penalty.
#'     \item \code{nu} (default: \code{NULL}) — Degrees of freedom for the prior distribution. Required if using a fixed predefined penalty. Must be an integer and compatible with the chosen \code{param}. For example, \code{(nu + q + 1)/q} must be an integer if \code{param = "variance"} and \code{(nu - q - 1)/q} must be an integer if \code{param = "precision"}.
#'     \item \code{const} (default: \code{1e8}) — A constant used in data augmentation to set the value of precision weights for pseudo-observations. Larger values may provide better approximation of the penalty, but may introduce numerical instability. Small values may result in poor approximation of the penalty.
#'     \item \code{param} (default: \code{"variance"}) — Specifies the parametrization of the prior: either \code{"variance"} for an inverse Wishart prior on the variance-covariance matrix or \code{"precision"} for a Wishart prior on the precision matrix. For univariate random effect, \code{"log variance"} can also be specified.
#'   }
#' @param verbose Logical. If TRUE, prints details about the penalty used for fitting penalized model and messages about the penalized model fitting.
#' @param ... Additional arguments passed to `glmmTMB()`.
#'
#' @details
#' When the random effect structure is univariate, the prior distribution on the variance is Inverse Gamma distribution that equals to univariate Inverse Wishart distribution.
#' To implement an Inverse Gamma penalty (or some other parametrization; see \code{penOpt} parameter \code{param} for details) with specified shape (\eqn{\alpha}) and scale (\eqn{\beta}),
#' set \code{nu=}\eqn{2\alpha} and \code{psi= }\eqn{2\beta}.
#'
#' Augmented approach in Poisson model with log link is implemented by the offset approach to improve numerical stability. Poisson models with different link functions
#' are implemented using weights. For these models the recommended value is \code{cons=10^6}.
#'
#' @return The function returns a list with elements:
#'
#' \describe{
#'   \item{non_pen}{The non-penalized `glmmTMB` model.}
#'   \item{pen}{The penalized `glmmTMB` model fit with augmented data.}
#'   \item{tau}{Estimated or supplied shrinkage strength towards mean eigenvalue. If the penalty was prespecified the value is \code{NULL}.}
#'   \item{error_pen}{If penalized fitting fails, the function returns the original unpenalized glmmTMB fit with pen = NULL and a flag error_pen = TRUE to indicate the failure.}
#' }
#'
#'
#' @references
#'  To be added after the article is published.
#'
#' @examples
#' \dontrun{
#'
#' library(blmeco)
#' data(ellenberg)
#' ellenberg$gradient  <-  paste(ellenberg$Year, ellenberg$Soil)
#' ellenberg$water.z <- as.numeric(scale(ellenberg$Water))
#' glmmTMBaug(log(Yi.g) ~ water.z + I(water.z^2) +
#'           (water.z + I(water.z^2)|Species) + (1|gradient), data=ellenberg, family="gaussian")
#'
#'}
#' @importFrom glmmTMB glmmTMB glmmTMBControl
#' @importFrom stats as.formula formula getCall model.frame model.matrix
#'  model.offset model.response model.weights qchisq uniroot
#' @importFrom utils modifyList
#' @importFrom lme4 findbars nobars VarCorr ranef

#' @export

glmmTMBaug <- function(formula, data, family,
                       penOpt = list(),
                       verbose=TRUE,
                       ...) {

  dots <- list(...)
  if ("ziformula" %in% names(dots) && !isTRUE(all.equal(dots$ziformula, ~0))) {
    stop("Non-default ziformula is not supported.")
  }
  if ("dispformula" %in% names(dots) && !isTRUE(all.equal(dots$dispformula, ~1))) {
    stop("Non-default dispformula is not supported.")
  }

  penOpt <- modifyList(
    list(tau = NULL, trunc=c(10^-4, 10^4), alpha = 0.05, psi = NULL, nu = NULL, const = 1e8, param = "variance"),
    penOpt
  )

  if (!is.null(penOpt$tau) && (!is.numeric(penOpt$tau) || length(penOpt$tau) != 1 || penOpt$tau < 0 || penOpt$tau > 1)) {
    stop("'tau' must be either NULL or a number on [0,1] interval.")
  }

  if (is.null(penOpt$const) && (!is.numeric(penOpt$const) || length(penOpt$const) != 1 || penOpt$const < 0)){
    stop("'const' has to be a positive number.")
  }

  if (penOpt$const < 1e4) {
    warning("Smaller value of 'const' may result in poor penalty approximation. The recomended value of const is 10^6 or higher.")
  }

  tau_spec <- penOpt$tau

  model <- withCallingHandlers(
    tryCatch({
      glmmTMB::glmmTMB(formula = formula, data = data, family = family, ...)
    }, error = function(e) {
      stop("Non-penalized fitting failed: ", conditionMessage(e))
      NULL
    }),
    warning = function(w) {
      message("Non-penalized fitting warning: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )


  allowed_families <- c("binomial", "poisson", "gaussian")
  if (!(model$modelInfo$family$family %in% allowed_families)) {
    warning(paste0("The ", model$modelInfo$family$family,
                   " family is not supported for penalized fitting. Returning the unpenalized model."))
    return(list(non_pen = model, pen = NULL, tau = NULL))
  }

  data_driven <- is.null(penOpt$psi) && is.null(penOpt$nu)

  pen_fit <- withCallingHandlers(
    tryCatch({
  if (!data_driven) {
    fit <- fit_augmented(model, data_driven = data_driven, penOpt = penOpt, ...)
    tau <- NULL
  } else {
    if (!is.null(penOpt$tau)) {
      fit <- fit_augmented(model, data_driven = data_driven, penOpt = penOpt, ...)
      tau <- penOpt$tau
    } else {
      # Try initial tau = 1
      penOpt_temp <- penOpt
      penOpt_temp$tau <- 1
      fit_tau1 <- fit_augmented(model, data_driven = data_driven, penOpt = penOpt_temp, ...)

      if (is.null(penOpt$alpha) && (!is.numeric(penOpt$alpha) || length(penOpt$alpha) != 1 || penOpt$alpha < 0 || penOpt$alpha > 1)) {
        stop("'alpha' must be a number on [0,1] interval. It should be specified for data-driven approach to determine 'tau'. ")
      }

      if(ncol(get_recovmat(model)[[1]])>1){
        crit <- abs(get_margLik_glmmtmb(model, fit_tau1) - model$fit$objective) >
          qchisq(1 - penOpt$alpha, df = 1) / 2
        }else{
        crit <- FALSE
      }

      if (crit) {
        tau_finder <- function(tau) {
          penOpt_local <- penOpt
          penOpt_local$tau <- tau
          fit <- fit_augmented(model, data_driven = TRUE, penOpt = penOpt_local, ...)
          abs(get_margLik_glmmtmb(model, fit) - model$fit$objective) -
            qchisq(1 - penOpt$alpha, df = 1) / 2
        }

        opt_tau <- try(uniroot(f = tau_finder, interval = c(0, 1)), silent = TRUE)

        if (inherits(opt_tau, "try-error")) {
          penOpt$tau <- 1
        } else {
          penOpt$tau <- opt_tau$root
        }
        fit <- fit_augmented(model, data_driven = data_driven, penOpt = penOpt, ...)

        tau <- penOpt$tau
      } else {
        fit <- fit_tau1
        tau <- 1
        if(ncol(get_recovmat(model)[[1]])==1) tau <- 0
      }
    }
  }

    list(fit=fit, tau=tau, error_pen=FALSE)

    }, error = function(e) {
      message("Penalized fitting failed: ", conditionMessage(e))
      list(fit=NULL, tau=NULL, error_pen=TRUE)
  }),
  warning = function(w) {
    if (!grepl("non-integer counts in a", conditionMessage(w))) {
      message("Penalized fitting warning: ", conditionMessage(w))
    }
    invokeRestart("muffleWarning")
  }
  )

  if (verbose & !pen_fit$error_pen) {
    msg <- paste0(
      "Penalized estimation summary:\n",
      "  - prior source: ",
      if (!data_driven) {
        paste0("predefined", ifelse(penOpt$param=="precision", " (Wishart)", " (Inverse Wishart)"))
      } else{
        paste0(paste0("data-driven", ifelse(penOpt$param=="precision", " (Wishart)", " (Inverse Wishart)")),
               " \n  - tau: ",
               ifelse(is.null(tau_spec),
                      paste0("estimated to ", format(pen_fit$tau, digits = 4)),
                      paste0("predefined to ", format(pen_fit$tau, digits = 4))))
      }, "\n",
      "  - parametrization: ", penOpt$param,
      "\n",
      "  - constant: ", penOpt$const
    )
    message(msg)
  }

  return(list(non_pen = model, pen = pen_fit$fit, tau=pen_fit$tau, error_pen = pen_fit$error_pen))

}

