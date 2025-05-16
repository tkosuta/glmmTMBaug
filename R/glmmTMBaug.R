#' Fit Penalized Generalized Linear Mixed Model via Data Augmentation
#'
#' This function extends `glmmTMB` to fit generalized linear mixed models (GLMMs)
#' with a penalty on the random effects covariance matrix via a data augmentation approach.
#'
#' @param formula A model formula specifying fixed and random effects, as in `glmmTMB`.
#' @param data A data frame containing the variables used in the model.
#' @param family family a family function, a character string naming a family function, or the result of a call to a family function (variance/link function) information. Only binomial(), poison() and gaussian() are currently supproted
#' @param data_driven Logical; if `TRUE`, estimates the penalty strength (`tau`) from the data.
#' @param penOpt A list of penalty options. Common elements include `tau`, `alpha`, `psi`, `nu`, `const`, `param`.
#' @param verbose Logical. If TRUE, prints details about the penalty used for fitting penalized model.
#' @param ... Additional arguments passed to `glmmTMB()`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{non_pen}{The unpenalized `glmmTMB` model.}
#'   \item{pen}{The penalized `glmmTMB` model fit with augmented data.}
#'   \item{tau}{Estimated or supplied penalty strength.}
#' }
#'
#' @examples
#' \dontrun{
#' data(sleepstudy, package = "lme4")
#' glmmTMB_aug(Reaction ~ Days + (Days|Subject), data = sleepstudy, family = gaussian())
#' }
#'
#' @importFrom glmmTMB glmmTMB glmmTMBControl
#' @importFrom stats as.formula formula getCall model.frame model.matrix
#'  model.offset model.response model.weights qchisq uniroot
#' @importFrom utils modifyList
#' @importFrom lme4 findbars nobars VarCorr ranef

#' @export

glmmTMBaug <- function(formula, data, family,
                       data_driven = TRUE,
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
    list(tau = NULL, alpha = 0.05, psi = NULL, nu = NULL, const = 1e6, param = "variance"),
    penOpt
  )
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
        fit <- fit_augmented(model, data_driven = TRUE, penOpt = penOpt, ...)

        tau <- penOpt$tau
      } else {
        fit <- fit_tau1
        tau <- 1
      }
    }
  }


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

  return(list(non_pen = model, pen = fit, tau=tau))
}
