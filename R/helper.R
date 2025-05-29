get_fixed_formula <- function(model) {
  fixed_form <- nobars(formula(model))
  fixed_form_no_offset <- as.formula(trimws(gsub("\\+?\\s*offset\\([^)]*\\)", "", deparse(fixed_form))))
  return(fixed_form_no_offset)
}

get_random_formula <- function(model) {

  random_form <- findbars(formula(model))

  if (!is.null(random_form)) {
    re_list <- NULL; j=1
    for(i in random_form){
      random_parts <- i
      bar <- random_parts[[1]]
      lhs <- random_parts[[2]]
      rhs <- random_parts[[3]]

      if(length(as.character(rhs))>1) stop("Only simple random effects of type (REexpr1 | factor1) + (REexpr2 | factor2) + ... are supported. Expressions of type (REexpr | factor1:factor2) or (REexpr | factor1/factor2) are not supported.")

      res <- list(expr=as.formula(paste0("~", deparse(lhs))),
                  gr=as.character(rhs))

      re_list[[j]] <- res
      j=j+1

    }
    return(re_list)
  } else {
    return(NULL)
  }
}

get_recovmat <- function(model){

  q <- sapply(ranef(model)$cond, ncol)
  mat <- lapply(1:length(ranef(model)$cond), function(x) matrix(VarCorr(model)$cond[[x]], ncol=q[x], nrow=q[x]))

  return(mat)
}

get_psi <- function(D_est, tau, trunc){
  if (!is.numeric(trunc) || length(trunc) != 2 || any(trunc <= 0)) {
    stop("'trunc' must be a numeric vector of two positive numbers.")
  }

  q <- ncol(D_est)
  ee <- eigen(D_est)
  min_trunc <- min(trunc)
  max_trunc <- max(trunc)
  ee$values[ee$values < min_trunc] <- min_trunc
  ee$values[ee$values > max_trunc] <- max_trunc
  lm <- mean(ee$values)
  li <- ee$values + tau * (lm - ee$values)
  psi <- ee$vectors %*% diag(li, ncol=q, nrow=q) %*% t(ee$vectors) * 3 * q
}

make_pseudo_data <- function(model_list, psi, nu, const=1e8, param="variance", linkinv=function(x) x){

  d <- ncol(model_list$X)

  if(ncol(model_list$Z1)>1){
    if (ncol(model_list$Z1)!=ncol(psi)) stop("'psi' does not match dimension of random effects")

    if (is.null(match.arg(param,c("precision", "variance")))) stop("'param' needs to be one of: precision,variance")

    q<-ncol(psi)
    if (param=="precision") {
      cc<-(nu-q-1)/q

      if(nu <  2*q+1) stop(paste0("Increase the value of 'nu'. Minimal 'nu' that can be implemented is: ", 2*q+1))
      if (!isTRUE(all.equal(cc %% 1, 0))) stop(paste0("'nu' must be such that (nu-",q,"-1)/",q," is an integer"))
      }
    if (param=="variance") {
      cc<-(nu+q+1)/q

      if(nu <  2*q-1) stop(paste0("Increase the value of 'nu'. Minimal 'nu' that can be implemented is: ", 2*q-1))
      if (!isTRUE(all.equal(cc %% 1, 0))) stop(paste0("'nu' must be such that (nu+",q,"+1)/",q," is an integer"))
    }

    true<-psi/cc
    ee<-eigen(true,TRUE)

    if(!all(ee$values>0)) stop("'psi' needs to be positive definite.")

    ui<-list()
    for (j in 1:q){
      ui[[j]]<-sqrt(ee$values[j])*ee$vectors[,j]
    }

    pi<-list()

    for (j in 1:length(ui)){
      I<-diag(rep(1,length(ui[[j]])))
      pi[[j]]<-linkinv(I%*%matrix(ui[[j]],ncol=1))
    }
    Y<-unlist(pi)

    id<- rep(1:q,each=q)

    Zi<-matrix(0,ncol=q,nrow=q)

    for (j in 1:q){
      Zi[j,j]<-1
    }
    for (j in 1:q){
      if (j==1) Z=Zi else Z<-rbind(Z,Zi)
    }

    fact<-cc
    if (fact>1){
      Y<-rep(Y,fact)
      id<-rep(1:(q*fact),each=q)
      for (j in 1:(q*fact)){
        if (j==1) Z=Zi else Z<-rbind(Z,Zi)
      }
    }
  }else{
    if (!(is.numeric(psi) && length(psi) == 1)) stop("'psi' must be a numeric scalar or a 1x1 matrix in the case of univariate random effects.")

    if (is.null(match.arg(param,c("precision","variance","log variance")))) stop("'param' needs to be one of: precision, variance or log variance")

    if (param=="precision") N<-max(c(floor(2*((nu/2)-1)),1))
    if (param=="variance") N<-max(c(floor(2*((nu/2)+1)),1))
    if (param=="log variance") N<-max(c(floor(2*(nu/2)),1))

    var.int <- (psi/2)*2/N
    fact <- N

    true <- matrix(var.int, ncol=1, nrow=1)

    ee <- eigen(true,TRUE)
    u1 <- sqrt(ee$values[1])*ee$vectors[,1]
    pi0 <- exp(u1[1])

    Y <- rep(c(pi0),fact)
    id <- c(1:fact)
    Z <- matrix(rep(1,fact),ncol=1)
  }

  pseudo_list <- list(Y=Y,
                      X=matrix(0, ncol=d, nrow=length(Y)),
                      offset=rep(0, length(Y)),
                      prec_weights=rep(const,length(id)),
                      freq_weights=rep(1, length(Y)),
                      Z1=Z,
                      gr1=max(model_list$gr1) + id
                      )

  num_re <- sum(grepl("Z", names(model_list)))

  if(num_re>1){
    for(i in 2:num_re){
      lml <- length(pseudo_list)
      pseudo_list <- append(pseudo_list,
                           eval(parse(text=paste0("list(Z",i,"=matrix(0, ncol=ncol(model_list$Z",i,"), nrow=length(Y)), gr",i,"=max(model_list$gr",i,") + id)"))),
                           after=lml)
    }
  }
  return(pseudo_list)
}

combine_two_lists <- function(x, y) {
  if (is.matrix(x) && is.matrix(y)) {
    return(rbind(x, y))
  } else {
    return(c(x, y))
  }
}

get_margLik_glmmtmb <- function(base_model, model){
  marg_llik <- base_model$obj$fn #take the negative log likelihood function (the objective function) from the model that was fitted only on the original data
  return(marg_llik(model$fit$par)) #evaluate the function at the parameters obtained from maximization of the (penalized) log likelihood
}

fit_augmented <- function(model, data_driven, penOpt = list(tau, psi, nu, const, param), ...) {

  mf <- model.frame(model)
  rand_formula <- get_random_formula(model)
  fix_formula <- get_fixed_formula(model)
  family <- model$modelInfo$family
  linkinv <- family$linkinv

  Y <- model.response(mf)
  X <- model.matrix(fix_formula, mf)
  for(i in 1:length(rand_formula)){
    eval(parse(text=paste0("Z",i, "<- model.matrix(rand_formula[[",i,"]][[1]], mf)")))
    eval(parse(text=paste0("gr",i, "<- as.numeric(factor(mf[, rand_formula[[",i,"]][[2]]]))")))
  }

  if (is.matrix(Y) && ((family$family != "binomial" && ncol(Y) > 1) || (ncol(Y) > 2))) {
    stop("Can't handle matrix-valued responses.")
  }

  mf$`(offset)` <- if (!is.null(model.offset(mf))) model.offset(mf) else 0

  if (is.null(model.weights(mf))) {
    mf$`(weights)` <- 1
    mf$prec_weights <- 1
    mf$freq_weights <- 1
  } else {
    if (any(family$family %in% c("binomial", "poisson"))) {
      mf$prec_weights <- mf$`(weights)`
      mf$freq_weights <- mf$`(weights)`
    } else {
      mf$prec_weights <- 1
      mf$freq_weights <- mf$`(weights)`
    }
  }

  model_list <- list(
    Y = Y,
    X = X,
    weights = mf$`(weights)`,
    offset = mf$`(offset)`,
    prec_weights = mf$prec_weights,
    freq_weights = mf$freq_weights
  )

  for(i in 1:length(rand_formula)){
    lml <- length(model_list)
    model_list <- append(model_list,
                         eval(parse(text=paste0("list(Z",i,"=Z",i,", gr",i,"=gr",i,")"))),
                         after=lml)
  }

  if (data_driven) {
    D_est <- get_recovmat(model)[[1]]
    psi <- get_psi(D_est=D_est, tau=penOpt$tau, trunc=penOpt$trunc)
    q <- ncol(psi)
    nu <- 2 * q - 1
    const <- penOpt$const
    param <- penOpt$param
  } else {
    psi <- penOpt$psi
    nu <- penOpt$nu
    const <- penOpt$const
    param <- penOpt$param
  }

  if(is.null(psi)) stop("Argument 'psi' is missing. If using the data-driven approach, set psi = NULL and nu = NULL. If specifying the penalty parameters manually, provide both 'psi' and 'nu'.")
  if(is.null(nu))  stop("Argument 'nu' is missing. If using the data-driven approach, set psi = NULL and nu = NULL. If specifying the penalty parameters manually, provide both 'psi' and 'nu'.")

  pseudo_list <- make_pseudo_data(
    model_list = model_list,
    psi = psi,
    nu = nu,
    const = const,
    param = param,
    linkinv = linkinv
  )

  if (is.matrix(Y) && family$family == "binomial") {
    n_events <- round(pseudo_list$Y * pseudo_list$prec_weights)
    pseudo_list$Y <- cbind(n_events, pseudo_list$prec_weights - n_events)
  }

  augmented_list <- lapply(names(pseudo_list), function(x) combine_two_lists(model_list[[x]], pseudo_list[[x]]))
  names(augmented_list) <- names(pseudo_list)

  call <- getCall(model)
  call_pen <- call

  form <- "Y~-1+X"
  for(i in 1:length(rand_formula)){
    form <-  paste0(form, "+(-1+Z", i, "|gr",i,")")
  }

  call_pen$formula <- as.formula(form)
  call_pen$data <- as.symbol("augmented_list")
  if(!is.null(model.offset(mf))) call_pen$offset <-  as.symbol("offset")
  if(any(family$family %in% c("binomial", "poisson"))){
    if(is.matrix(Y)){
      call_pen$weights <- as.symbol("freq_weights")
    }else{
      call_pen$weights <- as.symbol("prec_weights")
    }
  }

  # if (family$family == "poisson" && family$link == "log") {
  #   augmented_list$Y <- augmented_list$Y*augmented_list$prec_weights
  #   augmented_list$offset <- log(augmented_list$prec_weights)+augmented_list$offset
  #   call_pen$weights <- as.symbol("freq_weights")
  # }

  if (family$family == "gaussian") {
    call_pen$weights <- as.symbol("freq_weights")
    call_pen$dispformula <- as.formula("~ offset(-log(prec_weights))")
  }

  temp_env <- new.env(parent = parent.frame())
  temp_env$augmented_list <- augmented_list

  pen_model <- eval(call_pen, envir = temp_env)

  #pen_model$fit$objective <- model$obj$fn(pen_model$fit$par)
  #pen_model$frame <- model$frame
  #pen_model$modelInfo <- model$modelInfo
  #pen_model$call <- call
  #pen_model$obj$env$data <- model$obj$env$data

  return(pen_model)
}
