##------------------------------------------------------
## load required libraries here
##------------------------------------------------------

require(parallel)
require(Matrix)
require(splines)
require(nlme)
require(lme4)
#require(randomForestSRCM)

##------------------------------------------------------
## main function: boostmtree
##
## x             data frame containing the x-values
## tm            time values
## id            subject id
## y             outcome
## M             number of boosting iterations
## M.burn        number of burn-in iterations (applies only when lambda is estimated)
## nu            boosting regularization parameter (must be in [0,1])
##               can be a vector of length two (for b0 and b1)
## K             desired number of terminal nodes
## nknots        number of knots used
## d             degree of the piecewise B-spline polynomial (d=0 or d<1 kills the time effect)
## lambda        penalty parameter (if missing, or non-positive, estimated using mixed models)
## lambda.iter   number of iterations for iterative lambda estimation
## pen.ord       differencing order used to define the penalty
## svd.tol       tolerance used in svd of penalty matrix
## lambda.max    tolerance used for adaptively estimated lambda (caps it)
## forest.tol    tolderance used for forest weighted least squares solution
## verbose       (logical) terminal output?
##------------------------------------------------------

boostmtree <- function(x,
                       tm,
                       id,
                       y,
                       M = 200,
                       M.burn = 0,
                       nu = 0.01,
                       K = 3,
                       nknots = 10,
                       d = 3,
                       lambda,
                       lambda.iter = 5,
                       pen.ord = 2,
                       svd.tol = 1e-6,
                       lambda.max = 1e3,
                       forest.tol = 1e-3,
                       verbose = TRUE,
                       ...)
{

  ##------------------------------------------------------
  ## parse the data: processing/initialization
  ## set various dimensions/terms required later
  ## center y-values for numerical stability
  ## the y-centered mean must be added back to b0
  ## the data is sorted on id: thus we sort(id)
  ##------------------------------------------------------

  id.unq <- sort(unique(id))
  n <- length(id.unq)
  p <- ncol(x)
  xvar.names <- colnames(x)
  X <- do.call(rbind, lapply(1:n, function(i) {
    x[id == id.unq[i],, drop = FALSE][1,, drop = FALSE]}))
  x <- do.call(rbind, lapply(1:n, function(i) {x[id == id.unq[i],, drop = FALSE]}))
  Ymean <- mean(y, na.rm = TRUE)
  Y <- lapply(1:n, function(i) {y[id == id.unq[i]] - Ymean})
  ni <- unlist(lapply(1:n, function(i) {sum(id == id.unq[i])}))
  id <- sort(id)
  user.option <- match.call(expand.dots = TRUE)

  ##------------------------------------------------------
  ## define the Di matrix
  ##------------------------------------------------------

  ## define the unique time points
  tm.unq <- sort(unique(tm))
  n.tm <- length(tm.unq)
  tm.id <- lapply(1:n, function(i) {
    tm.i <- tm[id == id.unq[i]]
    match(tm.i, tm.unq)
  })

  ## piecewise B-spline polynomials
  if (nknots < 0) {
    warning("bsplines require a positive number of knots: eliminating b-spline fitting")
    d <- 0
  }
  if (d > 0) {
    if (n.tm > 1) {
      bs.tm <- bs(tm.unq, df = nknots + d, degree = d)
      X.tm <- cbind(1, bs.tm)
      attr(X.tm, "knots") <- attr(bs.tm, "knots")
      attr(X.tm, "Boundary.knots") <- attr(bs.tm, "Boundary.knots")
    }
    else {
      X.tm <- cbind(1, cbind(tm.unq))
    }
  }
  ## there is no time-componentt!!
  else {
    X.tm <- cbind(rep(1, n.tm))
    lambda <- 0
  }


  ## crucial dimension: the number of columns in the Di matrix used for modeling time
  ## if df.D == 1 there is no time component
  df.D <- ncol(X.tm)

  ## for each i extract the Di matrix
  D <- lapply(1:n, function(i) {
    ## we allow replicated time measurements for an individual
    ## therefore the match of tm.i to tm.unq is delicate
    cbind(X.tm[tm.id[[i]],, drop = FALSE])
  })

  ##------------------------------------------------------
  ## regularization details
  ## allows different nu values for b0 and b1
  ##------------------------------------------------------

  nu <- {if (length(nu) > 1) nu else rep(nu, 2)}
  if (sum(!(0 <= nu & nu <= 1)) > 0) {
    stop("regularization parameter (nu) must be in [0,1]")
  }
  nu.vec <- c(nu[1], rep(nu[2], df.D - 1))

  ##------------------------------------------------------
  ## tree/forest base-learner parameters used for rfsrc
  ##------------------------------------------------------

  ntree <- is.hidden.ntree(user.option)
  bootstrap <- is.hidden.bootstrap(user.option)
  if (ntree == 1) {
    nodesize = max(1, round(n/(2 * K)))
    mtry <- df.D + p
  }
  else {
    nodesize <- max(1, round(n/K))
    mtry <- NULL
    if (missing(lambda) || lambda < 0) {
      lambda <- 0
    }
  }

  ##------------------------------------------------------
  ## define the learner (used for setting the class)
  ##------------------------------------------------------

  if (ntree > 1) {
    learnerUsed <- "mforest"
  }
  else {
    if (df.D == 1) {
      learnerUsed <- "mtree"
    }
    else {
      learnerUsed <- "mtree-Pspline"
    }
  }


  ##------------------------------------------------------
  ## penalization/lambda details
  ##------------------------------------------------------

  ## determine if lambda can be estimated
  ## define the penalty matrices
  lambda.est.flag <- FALSE

  pen.lsq.matx <- penBSderiv(df.D - 1, pen.ord)
  if (ntree == 1 && (missing(lambda) || lambda < 0)) {
    if (df.D >= (pen.ord + 2)) {
      lambda.est.flag <- TRUE
      ## use svd to get the square root and inverse square root of P
      ## there may be numerical issues with the d entries: caution
      pen.mix.matx <- penBS(df.D - 1, pen.ord)
      svd.pen <- svd(pen.mix.matx)
      d.zap <- svd.pen$d < svd.tol
      d.sqrt <- sqrt(svd.pen$d)
      d.sqrt[d.zap] <- 0
      d.inv.sqrt <- 1 / sqrt(svd.pen$d)
      d.inv.sqrt[d.zap] <- 0
      pen.inv.sqrt.matx <- svd.pen$v %*% (t(svd.pen$v) * d.inv.sqrt)
    }
    else {
      warning("not enough degrees of freedom to estimate lambda: setting lambda to zero\n")
      lambda <- 0
    }
  }

  ## use burn-in?
  if (!lambda.est.flag) {
    M.burn <- 0
  }

  ##------------------------------------------------------
  ## initialization
  ## set various dimensions/terms that will be required
  ## final parameter checks
  ## multivariate tree formula/details
  ##------------------------------------------------------

  mu <- lapply(1:n, function(i) {rep(0, ni[i])})
  beta <- matrix(0, n, df.D)
  if (ntree == 1) {
    baselearner <- membership.list <- gamma.list <- vector("list", length = M)
  }
  else {
   membership.list <- gamma.list <- NULL
   baselearner <- vector("list", length = M)
  }
  sigma.vec <- phi.vec <- rho.vec <- rep(0, M)
  sigma.initial <- var(unlist(Y), na.rm = TRUE)
  phi <- 1

  ## rho initialization
  rho.fit.flag <- is.hidden.rho(user.option)
  if (rho.fit.flag == TRUE) {
    rho <- 0
  }
  else {
    rho <- rho.fit.flag
    rho.fit.flag <- FALSE
    if (rho < 0 || rho > 1) {
      stop("user specified rho is not valid:", rho)
    }
  }

  ## initialize sigma (lambda * phi * (1-rho))
  sigma <- var(unlist(Y), na.rm = TRUE)
  if (!lambda.est.flag) {
    sigma <- lambda * var(unlist(Y), na.rm = TRUE) * (1 - rho)
  }

  Y.names <- paste("Y", 1:df.D, sep = "")
  rfsrc.f <- as.formula(paste("Multivar(", paste(Y.names, collapse = ","), paste(") ~ ."), sep = ""))

  ##------------------------------------------------------
  ## MAIN LOOP
  ##
  ## gradient boosting loop
  ##
  ## MAIN LOOP
  ##------------------------------------------------------

  if (verbose) cat("  implementing multivariate boosting...\n")

  for (m in 1:(M.burn + M)) {

    if (verbose) {
      if (m <= M.burn) {
        cat("\t-- burn-in iteration:", m, "\n")
      }
      else {
        cat("\t-- iteration:", m - M.burn, "\n")
      }
    }

    ##---------------------------------------------------------
    ## step 3: calculate the negative gradient
    ## we discard the unnecessary (phi*(1-rho))^{-1} constant
    ##---------------------------------------------------------

    gm <- t(matrix(unlist(lapply(1:n, function(i) {
      rmi <- rho /(1 + (ni[i] - 1) * rho)
      cmi <- rmi * sum(Y[[i]] - mu[[i]], na.rm = TRUE)
      t(D[[i]]) %*% (Y[[i]] - mu[[i]] - cbind(rep(cmi, ni[i])))
    })), nrow = df.D))

    ##---------------------------------------------------------
    ##step 4(a): fit a K-terminal node multivariate regression tree

    ## gm    n x df.D  gradient-matrix
    ## X     n x p data frame containing x values
    ## K     desired number of terminal nodes

    ## Currently pruning is only available in predict mode.
    ## To recover the pruned membership of the training data,
    ## we must use the option "ptn.count"

    ## Thus we first we grow a tree with > K terminal nodes.  We then
    ## prune the tree back to K terminal nodes.  For convenience, we
    ## roughly calculate the initial nodesize such that the number of
    ## terminal nodes is 2 * K.

    ##
    ## If ntree > 1, we don't worry about pruning and grow a forest with
    ## roughly K terminal nodes
    ##
    ##---------------------------------------------------------

    incoming.data <- cbind(gm, X)
    names(incoming.data) = c(Y.names, names(X))

    ## multivariate forest learner
    if (ntree > 1) {
      rfsrc.obj <- rfsrc(rfsrc.f,
                         data = incoming.data,
                         mtry = mtry,
                         nodesize = nodesize,
                         importance = "none",
                         bootstrap = bootstrap,
                         ntree = ntree,
                         forest.wt = TRUE)

      Kmax <- n


    }

    ## multivariate tree learner
    else {
      rfsrc.obj <- rfsrc(rfsrc.f,
                         data = incoming.data,
                         ntree = 1,
                         mtry = mtry,
                         nodesize = nodesize,
                         importance = "none",
                         bootstrap = bootstrap)

      ## the outcome.target is irrelevant
      result.pred <- predict.rfsrc(rfsrc.obj,
                                   membership = TRUE,
                                   ptn.count = K,
                                   importance = "none")


      ## The pruned membership is an n x 1 vector of non-consecutive
      ## integers.  They map immutably to NODE_ID of the tree.
      membership <- membership.org <- c(result.pred$ptn.membership)
      if (m > M.burn) {
        membership.list[[m - M.burn]] <- membership.org
      }
      ## Recode the membership to be sequential
      membership <- as.numeric(factor(membership))
      ptn.id <- unique(membership)
      ## Kmax is the number of pseudo-terminal nodes, where Kmax <= K
      Kmax <-  length(ptn.id)

    }


    ## save the base learner
    if (m > M.burn) {
      baselearner[[m - M.burn]] <- rfsrc.obj
    }


    ##############################################################################
    ##
    ##  At this point the algorithm diverges depending upon ntree
    ##
    if (ntree == 1) {#####SINGLE TREE BASE LEARNER#####
    ##
    ##
    #############################################################################



      ##---------------------------------------------------------
      ## ESTIMATE LAMBDA USING MIXED MODELS: ntree = 1
      ## step 4(b)
      ##
      ## Ynew       transformed Y
      ## Xnew       transformed fixed effects design matrix
      ## Znew       transformed random effects desgin matrix
      ## gamma      gamma solution
      ##---------------------------------------------------------

      if (lambda.est.flag) {

        ## ---------------------------------------------------------
        ## transform the X, Y, Z values
        ## ---------------------------------------------------------

        transf.data <- mclapply(1:n, function(i) {
          if (ni[i] > 1) {
            ci <- const.sqrt(ni[i], rho)
            ##!!!!!!!!!!THIS CAN CREATE INSTABILITY!!!!!!!!!
            R.inv.sqrt <- (diag(1, ni[i]) - matrix(ci, ni[i], ni[i])) / sqrt(1 - rho)
          }
          else {
            R.inv.sqrt <- cbind(1)
          }
          Ynew <- R.inv.sqrt %*% (Y[[i]] - mu[[i]])
          Xnew <- R.inv.sqrt %*% D[[i]][, 1, drop = FALSE]
          Znew <- R.inv.sqrt %*% D[[i]][, -1, drop = FALSE] %*% pen.inv.sqrt.matx
          list(Ynew = Ynew, Xnew = Xnew, Znew = Znew)
        })

        ## ----------------------------------------------------------------------
        ## iterate to estimate gamma and sigma (= lambda * phi)
        ## ----------------------------------------------------------------------

        ## initialize sigma
        sigma.hat <- sigma.initial

        ## iterative loop
        for (k in 1:lambda.iter) {

          ## BLUP for random effects (u_k) and fixed effects (alpha_k) conditional on sigma
          blup.obj <-  blup.solve(transf.data, membership, sigma.hat, Kmax)


          ## sigma method of moments estimator conditional on BLUP
          ## robust rss calculation used to avoid deflating sigma estimate
          sigma.obj <- lapply(1:Kmax, function(k) {
            pt.k <- (membership == k)
            Z <- do.call(rbind, lapply(which(pt.k), function(j) {transf.data[[j]]$Znew}))
            X <- do.call(rbind, lapply(which(pt.k), function(j) {transf.data[[j]]$Xnew}))
            Y <- unlist(lapply(which(pt.k), function(j) {transf.data[[j]]$Ynew}))
            ZZ <- t(Z) %*% Z
            rss <- (Y - X %*% c(blup.obj[[k]]$fix.eff))^2
            robust.pt <- (rss <= quantile(rss, .99, na.rm = TRUE))
            rss <- sum(rss[robust.pt], na.rm = TRUE)
            resid <- (Y - X %*% c(blup.obj[[k]]$fix.eff) - Z %*% c(blup.obj[[k]]$rnd.eff))^2
            resid <- resid[robust.pt]
            return(list(trace.Z = sum(diag(ZZ)), rss = rss, resid = resid))
          })
          num <- sum(unlist(lapply(1:Kmax, function(k) {sigma.obj[[k]]$trace.Z})), na.rm = TRUE)
          den <- sum(unlist(lapply(1:Kmax, function(k) {sigma.obj[[k]]$rss})), na.rm = TRUE)
          N <- sum(unlist(lapply(1:Kmax, function(k) {sigma.obj[[k]]$resid})), na.rm = TRUE)

          if (!is.na(den) && den > (.99 * N)) {## ensure that denominator is larger than degrees of freedom
            sigma.hat <- num / (den - .99 * N)
          }
          else {
            sigma.hat <- min(sigma.hat, lambda.max * phi)
          }

          ## cap sigma.hat
          sigma.hat <- min(sigma.hat, lambda.max * phi)

          ## update BLUP on last iteration, conditional on sigma
          if (k == lambda.iter) {
            blup.obj <-  blup.solve(transf.data, membership, sigma.hat, Kmax)
          }

        }

        ## update gamma (don't forget to transform the random effects!)
        gamma <- lapply(1:Kmax, function(k) {
          c(blup.obj[[k]]$fix.eff, pen.inv.sqrt.matx %*% blup.obj[[k]]$rnd.eff)
        })

      }

      ##---------------------------------------------------------
      ## LEAST SQUARES SOLUTION FOR GAMMA: ntree = 1
      ## applies only when mixed models are not used
      ## step 4(b): solve for gamma
      ##
      ## Xnew      pseudo x matrix
      ## YnewSum   summed pseudo y value
      ## XnewSum   summed pseudo x values
      ## gamma     weighted least squares solution
      ##---------------------------------------------------------

      if (!lambda.est.flag) {
        Xnew <- mclapply(1:n, function(i) {
          rmi <- rho / (1 + (ni[i] - 1) * rho)
          Wi <- diag(1, ni[i]) - matrix(rmi, ni[i], ni[i])
          t(D[[i]]) %*% Wi %*% D[[i]]
        })
        gamma <- lapply(1:Kmax, function(k) {
          pt.k <- (membership == k)
          ##sum the pseudo y's over a given terminal node
          YnewSum <- colSums(gm[pt.k,, drop = FALSE])
          ##sum the pseudo x's over a given terminal node
          XnewSum <- Reduce("+", lapply(which(pt.k), function(j) {Xnew[[j]]}))
          ## add penalization
          XnewSum <- XnewSum + sigma * pen.lsq.matx
          ## solve using QR
          qr.obj <- tryCatch({qr.solve(XnewSum, YnewSum)}, error = function(ex){NULL})
          if (!is.null(qr.obj)) {
            qr.obj
          }
          else {
            rep(0, df.D)
          }
        })
      }


      ##---------------------------------------------------------
      ## step 4(c): save gamma: ntree = 1
      ## needed for prediction
      ## save as a matrix and include original grow terminal node membership
      ##---------------------------------------------------------

      gamma.matx <- matrix(0, Kmax, df.D+1)
      gamma.matx[, 1] <- sort(unique(membership.org))
      gamma.matx[, 2:(df.D+1)] <- matrix(unlist(gamma), ncol = df.D, byrow = TRUE)

      if (m > M.burn) {
        gamma.list[[m - M.burn]] <- gamma.matx
      }

      ##---------------------------------------------------------
      ## step 5: update beta and mu: ntree = 1
      ##---------------------------------------------------------

      ## reset beta if a burn-in is used
      ## shut-off lamda adaptivity
      if (lambda.est.flag && m == (M.burn + 1)) {
        beta <- matrix(0, n, df.D)
        if (M.burn > 0) {
          lambda <- sigma.hat / phi ####TBD
          lambda.est.flag <- FALSE
        }
      }

      ## beta update
      beta.update <- matrix(unlist(lapply(1:n, function(i) {
        gamma[[membership[i]]]})), nrow = df.D)
      beta <- beta + t(beta.update * nu.vec)
      ## mu update
      mu <- lapply(1:n, function(i) {D[[i]] %*% beta[i, ]})

    }


    ##############################################################################
    ##
    ##
    else{ #####FOREST BASE LEARNER ####
    ##
    ##
    #############################################################################

      ##---------------------------------------------------------
      ## WEIGHTED LEAST SQUARES SOLUTION FOR BETA WHERE WEIGHTS
      ## ARE OBTAINED FROM THE FOREST
      ##
      ## Xnew      pseudo x matrix
      ## YnewSum   summed pseudo y value
      ## XnewSum   summed pseudo x values
      ## bhat      weighted least squares solution
      ##---------------------------------------------------------


      ## extract the forest weight
      forest.wt <- rfsrc.obj$forest.wt

      ## define Xnew
      Xnew <- mclapply(1:n, function(i) {
        rmi <- rho / (1 + (ni[i] - 1) * rho)
        Wi <- diag(1, ni[i]) - matrix(rmi, ni[i], ni[i])
        t(D[[i]]) %*% Wi %*% D[[i]]
      })

      ## iterate over cases i to get bhat
      ## for speed we eliminate cases with small forest weights

      bhat <- do.call("cbind", mclapply(1:n, function(i) {
        fwt.i <- forest.wt[i, ]
        fwt.i[fwt.i <= forest.tol] <- 0
        if (all(fwt.i == 0)) {
          fwt.i <- forest.wt[i, ]
        }
        fwt.i <- fwt.i / sum(fwt.i)
        pt.i <- (fwt.i != 0)
        ##sum the pseudo y's over non-zero weights: scale by the weights
        YnewSum <- colSums(fwt.i[pt.i] * gm[pt.i,, drop = FALSE])
        ##sum the pseudo x's over non-zero weights: scale by the weights
        XnewSum <- Reduce("+", lapply(which(pt.i), function(j) {fwt.i[j] * Xnew[[j]]}))
        ## add penalization
        XnewSum <- XnewSum + sigma * pen.lsq.matx
        ## solve using QR
        qr.obj <- tryCatch({qr.solve(XnewSum, YnewSum)}, error = function(ex){NULL})
        if (!is.null(qr.obj)) {
          qr.obj
        }
        else {
          rep(0, df.D)
        }
      }))

      ##---------------------------------------------------------
      ## update beta and  mu
      ##---------------------------------------------------------

      beta <- beta + t(bhat * nu.vec)
      mu <- lapply(1:n, function(i) {D[[i]] %*% beta[i, ]})

    }

    ##############################################################################
    ##
    ##  The algorithm is now the same
    ##
    #############################################################################


    ##---------------------------------------------------------
    ## step 6 and 7: update phi and rho using REML
    ## we make a call to lme (actually gls since there are no r.eff)
    ## we use a compound symmetric correlation matrix, but this can be generalized
    ## note: tm and x are used in the REML call -- although the theory does not call for it
    ##---------------------------------------------------------
    resid.data <- data.frame(y  = unlist(lapply(1:n, function(i) {Y[[i]] - mu[[i]]})),
                             x  = x,
                             tm = unlist(lapply(1:n, function(i) {tm[id == id.unq[i]]})),
                             id = unlist(lapply(1:n, function(i) {rep(id.unq[i], ni[i])})))
    gls.obj <- tryCatch({gls(y ~ ., data = resid.data, corr = corCompSymm(form = ~ 1 | id))},
                        error = function(ex){NULL})
    if (is.null(gls.obj)) {
      gls.obj <- tryCatch({gls(y ~ 1, data = resid.data, corr = corCompSymm(form = ~ 1 | id))},
                          error = function(ex){NULL})
    }
    if (!is.null(gls.obj)) {
      phi <- gls.obj$sigma^2
      if (rho.fit.flag) {
        rho <- as.numeric(coef(gls.obj$modelStruct$corStruc, unconstrained = FALSE))
        rho <- max(min(0.999, rho, na.rm = TRUE), -0.999)
      }
    }

    if (m > M.burn) {
      phi.vec[m - M.burn] <- phi
      rho.vec[m - M.burn] <- rho
    }

    if (verbose) {
      cat("phi   :", phi, "\n")
      cat("rho   :", rho, "\n")
    }

    ##---------------------------------------------------------
    ## step 8: update sigma (if lambda is not adaptively estimated)
    ##---------------------------------------------------------

    if (!lambda.est.flag) {
      sigma <- lambda * phi * (1 - rho)
      if (verbose) {
        cat("sigma :", sigma, "\n")
        cat("lambda:", sigma / (phi * (1-rho)), "\n")
      }
      if (m > M.burn) {
        sigma.vec[m - M.burn] <- sigma
      }
    }

    if (lambda.est.flag && m > M.burn) {
      if (verbose) {
        cat("sigma :", sigma.hat, "\n")
        cat("lambda:", sigma.hat / phi, "\n")
      }
      sigma.vec[m - M.burn] <- sigma.hat * (1 - rho)
    }


  }##end boosting iteration

  ##---------------------------------------------------------
  ## add the Y mean back
  ##---------------------------------------------------------

  Y <- lapply(1:n, function(i) {y[id == id.unq[i]]})
  if (ntree == 1) {
    beta[, 1] <- beta[, 1] + Ymean
  }
  mu <- lapply(1:n, function(i) {c(mu[[i]] + Ymean)})

  ##---------------------------------------------------------
  ## return the promised object
  ##---------------------------------------------------------

  obj <- list(x = X,
              time = lapply(1:n, function(i) {tm[id == id.unq[i]]}),
              id = id,
              y = Y,
              ymean = Ymean,
              gamma = gamma.list,
              beta = beta,
              mu = mu,
              lambda = (sigma.vec / (phi.vec * (1 - rho.vec))),
              phi = phi.vec,
              rho = rho.vec,
              baselearner = baselearner,
              membership = membership.list,
              D = X.tm,
              d = d,
              K = K,
              M = M,
              nu = nu,
              ntree = ntree)

  class(obj) <- c("boostmtree", "grow", learnerUsed)

  invisible(obj)

}

