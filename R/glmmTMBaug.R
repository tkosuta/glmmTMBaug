#' Fit Penalized Generalized Linear Mixed Model via Data Augmentation
#'
#' This function extends `glmmTMB` to fit generalized linear mixed models (GLMMs)
#' with a penalty on the random effects covariance matrix via a data augmentation approach.
#'
#' @param formula A model formula specifying fixed and random effects, as in `glmmTMB`.
#' @param data A data frame containing the variables used in the model.
#' @param family family a family function, a character string naming a family function, or the result of a call to a family function (variance/link function) information. Only binomial(), poison() and gaussian() are currently supproted.
#' @param penOpt A named list of penalty options used to control the penalized likelihood fit. If \code{psi=NULL} and \code{nu=NULL}, a data-driven approach is used to determine psi. If \code{tau=NULL} a data driven approach is used to determine tau.
#' If both \code{psi} and \code{nu} are specified, a prior with those parameters is used instead. In this case, the parameters \code{tau}, \code{trunc} and \code{alpha} are ignored.
#'   The following penalty options are recognized (with defaults):
#'   \itemize{
#'     \item \code{tau} (default: \code{NULL}) — Numeric value in the interval \[0, 1\]. If NULL, a data-driven method is used to estimate it.
#'     \item \code{trunc} (default: \code{c(1e-4, 1e4)}) — Lower and upper truncation bounds of eigenvalues in data-drive approach to determine \code{psi}.
#'     \item \code{alpha} (default: \code{0.05}) — Significance level used in the data-driven procedure for estimating \code{tau}.
#'     \item \code{psi} (default: \code{NULL}) — A positive-definite matrix specifying the prior scale matrix for the inverse Wishart (variance) or Wishart (precision) distribution. Required if using a fixed prior.
#'     \item \code{nu} (default: \code{NULL}) — Degrees of freedom for the prior distribution. Must be an integer and compatible with the chosen \code{param}. For example, \code{(nu + q + 1)/q} must be an integer if \code{param = "variance"} and \code{(nu - q - 1)/q} must be an integer if \code{param = "precision"}.
#'     \item \code{const} (default: \code{1e6}) — A constant used in data augmentation to set the value of precision weights for pseudo-observations.
#'     \item \code{param} (default: \code{"variance"}) — Specifies the parametrization of the prior: either \code{"variance"} for an inverse Wishart prior on the variance-covariance matrix or \code{"precision"} for a Wishart prior on the precision matrix.
#'   }
#' @param verbose Logical. If TRUE, prints details about the penalty used for fitting penalized model.
#' @param ... Additional arguments passed to `glmmTMB()`.
#'
#' @return A list with elements:
#'
#' \describe{
#'   \item{non_pen}{The unpenalized `glmmTMB` model.}
#'   \item{pen}{The penalized `glmmTMB` model fit with augmented data.}
#'   \item{tau}{Estimated or supplied shrinkage strength towards mean eigenvalue.}
#'   \item{error_pen}{If penalized fitting fails (e.g., due to prior mis-specification), the function returns the original unpenalized glmmTMB fit with pen = NULL and a flag error_pen = TRUE to indicate the failure.}
#' }
#'
#' @examples
#' \dontrun{
#' data(sleepstudy, package = "lme4")
#' glmmTMBaug(Reaction ~ Days + (Days|Subject), data = sleepstudy, family = gaussian())
#' }
#'
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

  rand_terms <- lme4::findbars(formula)
  if (length(rand_terms) != 1) {
    stop("Only models with a single random effect term are supported currently.")
  }

  penOpt <- modifyList(
    list(tau = NULL, trunc=c(10^-4, 10^4), alpha = 0.05, psi = NULL, nu = NULL, const = 1e6, param = "variance"),
    penOpt
  )

  if (!is.null(penOpt$tau) && (!is.numeric(penOpt$tau) || length(penOpt$tau) != 1 || penOpt$tau < 0 || penOpt$tau > 1)) {
    stop("'tau' must be either NULL or a number on [0,1] interval.")
  }

  if (!is.null(penOpt$alpha) && (!is.numeric(penOpt$alpha) || length(penOpt$alpha) != 1 || penOpt$alpha < 0 || penOpt$alpha > 1)) {
    stop("'alpha' must be either NULL or a number on [0,1] interval.")
  }

  if (is.null(penOpt$const) && (!is.numeric(penOpt$const) || length(penOpt$const) != 1 || penOpt$const < 0)){
    stop("'const' has to be a positive number.")
  }

  if (penOpt$const < 1e4) {
    warning("Smaller value of 'const' may result in poor penalty approximation. The recomended value of const is 10^6.")
  }

  tau_spec <- penOpt$tau

  model <- glmmTMB::glmmTMB(formula = formula, data = data, family = family, ...)

  allowed_families <- c("binomial", "poisson", "gaussian")
  if (!(model$modelInfo$family$family %in% allowed_families)) {
    warning(sprintf(
      "The '%s' family is not supported for penalized fitting. Returning the unpenalized model.",
      model$modelInfo$family$family
    ))
    return(list(non_pen = model, pen = NULL, tau = NULL))
  }

  data_driven <- is.null(penOpt$psi) && is.null(penOpt$nu)

  pen_fit <- tryCatch({
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

      crit <- abs(get_margLik_glmmtmb(model, fit_tau1) - model$fit$objective) >
        qchisq(1 - penOpt$alpha, df = 1) / 2

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
      }
    }
  }

    list(fit=fit, tau=tau, error_pen=FALSE)

    }, error = function(e) {
    if (verbose) {
      message("Penalized fitting failed: ", conditionMessage(e))
    }
     list(fit=NULL, tau=NULL, error_pen=TRUE)
  })

  if (verbose) {
    msg <- paste0(
      "Penalized estimation summary:\n",
      "  - prior source: ",
      if (!data_driven) {
        paste0("predefined", ifelse(penOpt$param=="precision", " (Wishart)", " (Inverse Wishart)"))
      } else{
        paste0(paste0("data-driven", ifelse(penOpt$param=="precision", " (Wishart)", " (Inverse Wishart)")),
               " \n  - tau: ",
               ifelse(is.null(tau_spec),
                      paste0("estimated to ", format(tau, digits = 4)),
                      paste0("predefined to ", format(tau, digits = 4))))
      }, "\n",
      "  - parametrization: ", penOpt$param,
      "\n",
      "  - constant: ", penOpt$const
    )
    message(msg)
  }

  return(list(non_pen = model, pen = pen_fit$fit, tau=pen_fit$tau, error_pen = pen_fit$error_pen))

}
