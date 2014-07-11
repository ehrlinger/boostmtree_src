##------------------------------------------------------
## prediction mode for boostmtree
##
## if test y is supplied then an optimal M is determined
## returned objects (such as VIMP) will depend on this
##
## obj         boosting object
## x           data frame of test set x values
## tm          time values (optional: but required if y is supplied)
## id          test set subject id (optional: but required if y is supplied)
## y           test set outcome (optional)
## M           (optional) fix the value of M (over-rides the optimized M if supplied)
## importance  (logical) calculate variable importance (VIMP)?
## proximity   (logical) determine proximity of test data to training data
## verbose     (logical) terminal output?
## eps         tolerance value for determining the optimal M
##------------------------------------------------------

generic.predict.boostmtree <- function(object,
                                       x,
                                       tm,
                                       id,
                                       y,
                                       M,
                                       importance = TRUE,
                                       proximity = FALSE,
                                       verbose = TRUE,
                                       eps = 1e-4,
                                       ...)
{

  ##------------------------------------------------------
  ## preliminary checks: all are fatal
  ##------------------------------------------------------

  if (missing(object)) {
    stop("object is missing!")
  }
  if (sum(inherits(object, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)'")
  }

  if (verbose) cat("  running prediction mode for multivariate boosting\n")

  ##------------------------------------------------------
  ## check if this is a partial plot call
  ##------------------------------------------------------
  user.option <- match.call(expand.dots = TRUE)
  partial <- is.hidden.partial(user.option)

  ##------------------------------------------------------
  ## parse the data: processing/initialization
  ## depends if a partial plot call was initiated
  ##------------------------------------------------------

  if (!partial) {## regular call

    ## the grow data is used if x is missing
    ## no y is assumed
    if (missing(x)) {
      X <- object$x
      n <- nrow(X)
      D <- object$D
      tm <- object$time
      tm.unq <- sort(unique(unlist(tm)))
      testFlag <- FALSE
    }
    else {
      ## if x is provided but no other values
      ## we assume id's are 1:1 and swap in the original unique time values
      ## no y is assumed
      if (!missing(x) && (missing(id) || missing(tm))) {
        X <- x
        n <- nrow(X)
        tm.unq <- sort(unique(unlist(object$time)))
        tm <- lapply(1:n, function(i){tm.unq})
        testFlag <- FALSE
      }
      ## we now assume all test set information is provided: fails otherwise
      else{
        if (missing(id)) {
          stop("test set id values are missing\n")
        }
        id.unq <- sort(unique(id))
        n <- length(id.unq)
        if (missing(x)) {
          stop("test set x values are missing\n")
        }
        X <- do.call(rbind, lapply(1:n, function(i) {
          x[id == id.unq[i],, drop = FALSE][1,, drop = FALSE]}))
        if (missing(tm)) {
          stop("test set time values are missing\n")
        }
        tm.unq <- sort(unique(tm))
        if (!missing(y)) {
          tm <- lapply(1:n, function(i) {tm[id == id.unq[i]]})
          Y <- lapply(1:n, function(i) {y[id == id.unq[i]]})
          testFlag <- TRUE
        }
        else {
          testFlag <- FALSE
        }
      }
    }
  }

  else {## partial plot call

      X <- x
      n <- nrow(X)
      tm.unq <- tm
      tm <- lapply(1:n, function(i){tm.unq})
      testFlag <- FALSE

  }


  ##------------------------------------------------------
  ## construct the test set bspline basis functions
  ##------------------------------------------------------

  if (object$d > 0) {
    if (length(tm.unq) > 1) {
      D <- cbind(1, bs(tm.unq, knots = attr(object$D, "knots"),
                       Boundary.knots = attr(object$D, "Boundary.knots"), degree = object$d))
    }
    else {
          stop("only one unique time point")
        }
  }
  else {
    D <- cbind(rep(1, length(tm.unq)))
  }


  ##------------------------------------------------------
  ## set dimensions
  ## extract values from the object
  ##------------------------------------------------------

  ## basic parameters
  if (missing(M)) {
    M <- object$M
    Mflag <- FALSE
  }
  else {
    M <- max(1, min(M, object$M))
    Mflag <- TRUE
  }
  K <- object$K
  nu <- object$nu
  ntree <- object$ntree
  p <- ncol(X)
  df.D <- ncol(D)
  xvar.names <- colnames(X)

  ## regularization details
  ## allows different nu values for b0 and b1
  nu.vec <- c(nu[1], rep(nu[2], df.D - 1))

  ## objects needed for updating beta
  Ymean <- object$ymean
  Ysd <- object$ysd
  gamma <- object$gamma
  baselearner <- object$baselearner

  ## initialize beta
  if (ntree == 1) {
    beta <- matrix(0, n, df.D)
    beta.list <-  beta.vimp <- vector("list", M)
  }
  else {
    beta.list <- NULL
  }

  ## vimp
  vimpFlag <- testFlag && importance
  vimp <- NULL

  ## proximity
  prox <- NULL

  ##############################################################################
  ##
  ##  Iterate over the M iterations to recursively define the test predictor
  ##  At this point the algorithm diverges depending upon ntree
  ##
  if (ntree == 1) {#####SINGLE TREE BASE LEARNER#####
  ##
  ##
  #############################################################################

    ## we make the predict call in parallel for speed
    ## to execute a parallel call we must over-ride rf.cores
    ## we must be careful not to over-tax mc.cores
    rf.cores.old <- getOption("rf.cores")
    mc.cores.old <- getOption("mc.cores")
    options(rf.cores = 0, mc.cores = max(1, detectCores() / 2))

    ## membership from original grow tree
    membership <- mclapply(1:M, function(m) {

      options(rf.cores = 0, mc.cores = 1)
      c(predict.rfsrc(baselearner[[m]],
                      newdata = X,
                      outcome.target = 1, ## the target is irrelevant
                      outcome = "train",
                      membership = TRUE,
                      ptn.count = K,
                      importance = "none")$ptn.membership)

    })

    ## membership for noised up data (if requested)
    if (vimpFlag) {
      membershipNoise <- mclapply(1:M, function(m) {


        Xnoise <- do.call(rbind, lapply(1:p, function(k) {
          X.k <- X
          X.k[, k] <- sample(X.k[, k])
          X.k
        }))

        options(rf.cores = 0, mc.cores = 1)
        c(predict.rfsrc(baselearner[[m]],
                        newdata = Xnoise,
                        outcome.target = 1, ## the target is irrelevant
                        outcome = "train",
                        membership = TRUE,
                        ptn.count = K,
                        importance = "none")$ptn.membership)
      })
    }

    ## proximity of test data to training data
    if (proximity) {

      prox <- Reduce("+", mclapply(1:M, function(m) {

        options(rf.cores = 0, mc.cores = 1)
        if (!testFlag) {
          predict.rfsrc(baselearner[[m]],
                        newdata = X,
                        outcome.target = 1, ## the target is irrelevant
                        membership = TRUE,
                        ptn.count = K,
                        importance = "none",
                        proximity = TRUE)$proximity

        }
        else {
          ntest <- nrow(X)
          predict.rfsrc(baselearner[[m]],
                        newdata = rbind(X, object$x),
                        outcome.target = 1, ## the target is irrelevant
                        membership = TRUE,
                        ptn.count = K,
                        importance = "none",
                        proximity = TRUE)$proximity[1:ntest, -(1:ntest)]
        }

      })) / M

    }

    ## reset the rf.cores/mc.cores option to its original value
    if (!is.null(rf.cores.old)) options(rf.cores = rf.cores.old)
    if (!is.null(mc.cores.old)) options(mc.cores = mc.cores.old)

    for (m in 1:M) {

      ##---------------------------------------------------------
      ## pass the x-data down the mth tree
      ## acquire terminal node membership (tnm)
      ## update beta using gamma making use of tnm
      ## rescale by Ysd and add Ymean
      ##---------------------------------------------------------

      orgMembership <- gamma[[m]][, 1]
      beta <- t(t(gamma[[m]][match(membership[[m]], orgMembership), -1, drop = FALSE]) * nu.vec) * Ysd

      if (m == 1) {
        beta[, 1] <- beta[, 1] + Ymean
        beta.list[[m]] <- beta
      }
      else {
        beta.list[[m]] <- beta.list[[m-1]] + beta
      }

      ##---------------------------------------------------------
      ## vimp calculation: update "vimp beta" from the perturbed X
      ##---------------------------------------------------------

      if (vimpFlag) {

        beta.vimp[[m]] <- lapply(1:p, function(k) {
          membership.k <- membershipNoise[[m]][((k-1) * n + 1):(k * n)]
          beta.update.k <- t(t(gamma[[m]][match(membership.k, orgMembership), -1, drop = FALSE]) * nu.vec) * Ysd
          if (m == 1) {
            beta.m <- beta.update.k
            beta.m[, 1] <- beta.m[, 1] + Ymean##add the overall mean
            beta.m
          }
          else {
            beta.vimp[[m-1]][[k]] + beta.update.k
          }
        })

      }

    }##loop is complete

    
  }

  ##############################################################################
  ##
  ##
  else{ #####FOREST BASE LEARNER ####
  ##
  ##
  #############################################################################

    stop("TBD TBD TBD")

  }

  ##---------------------------------------------------------
  ## determine a cumulative error rate
  ##---------------------------------------------------------

  muhat <- mclapply(1:M, function(m) {
    ## extract those mu values corresponding to the observed time points
    ## we allow replicated time measurements for an individual
    ## therefore the match of tm to tm.unq is delicate
    DbetaT.m <- D %*% t(beta.list[[m]])
    lapply(1:n, function(i) {
      DbetaT.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
    })
  })

  if (testFlag) {
    err.rate <- matrix(unlist(lapply(1:M, function(m) {
      c(l1Dist(Y, muhat[[m]]), l2Dist(Y, muhat[[m]]))
    })), ncol = 2, byrow = TRUE)
    colnames(err.rate) <- c("l1", "l2")
  }
  else {
    err.rate <- NULL
  }

  ##---------------------------------------------------------
  ## determine the optimal M: uses L2 error rate
  ##---------------------------------------------------------

  if (!Mflag && testFlag) {
    diff.err <- abs(err.rate[, "l2"] - min(err.rate[, "l2"], na.rm = TRUE))
    diff.err[is.na(diff.err)] <- 1
    if (sum(diff.err < eps) > 0) {
      Mopt <- min(which(diff.err < eps))
    }
    else {
      Mopt <- M
    }
  }
  else {
    Mopt <- M
  }


  ##---------------------------------------------------------
  ## variable importance
  ##---------------------------------------------------------

  if (vimpFlag) {
    vimp <- rep(0, p)
    vimp <- unlist(mclapply(1:p, function(k) {
      DbetaT.k <- D %*% t(beta.vimp[[Mopt]][[k]])
      ## extract those mu values corresponding to the observed time points
      ## we allow replicated time measurements for an individual
      ## therefore the match of tm to tm.unq is delicate
      muhat.k <- lapply(1:n, function(i) {
        DbetaT.k[, i][match(tm[[i]], tm.unq, tm[[i]])]
      })
      l2Dist(Y, muhat.k) - err.rate[Mopt, "l2"]
    }))
    names(vimp) <- xvar.names
  }

  ##---------------------------------------------------------
  ##return the promised object
  ##---------------------------------------------------------

  pobj <- list(
       boost.obj = object,
       x = X,
       time = tm,
       time.unq = tm.unq,
       y = if (testFlag) Y else NULL,
       beta = beta.list[[Mopt]],
       mu = muhat[[Mopt]],
       err.rate = err.rate,
       mse = err.rate[Mopt, "l2"],
       vimp = vimp,
       proximity = prox,
       Mopt = if (testFlag) Mopt else NULL)

  class(pobj) <- c("boostmtree", "predict", class(object)[3])

  invisible(pobj)


}

