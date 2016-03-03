##------------------------------------------------------
## main function: boostmtree
##
## x             data frame containing the x-values
## tm            time values
## id            subject id
## y             outcome
## M             number of boosting iterations
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
## forest.tol    tolerance used for forest weighted least squares solution
## verbose       (logical) terminal output?
##------------------------------------------------------

boostmtree <- function(x,
                       tm,
                       id,
                       y,
                       M = 200,
                       nu = 0.05,
                       K = 5,
                       nknots = 10,
                       d = 3,
                       pen.ord = 3,
                       lambda,
                       lambda.max = 1e6,
                       lambda.iter = 2,
                       importance = FALSE,
                       svd.tol = 1e-6,
                       forest.tol = 1e-3,
                       verbose = TRUE,
                       ...)
{

  ##------------------------------------------------------
  ## is this a univariate setting?
  ## if so, flag it and make a zero time vector
  ## we also set d < 0 to trigger no time-covariate fitting
  ##------------------------------------------------------
  univariate <- FALSE
  id.unq <- sort(unique(id))
  n <- length(id.unq)
  if (length(id.unq) == length(id)) {
    univariate <- TRUE
    tm <- rep(0, n)
    d <- -1
  }
  
  ##------------------------------------------------------
  ## parse the data: processing/initialization
  ## set various dimensions/terms required later
  ## center and scale y-values for numerical stability
  ## the y-centered mean and scaling values must be corrected later
  ## data is sorted on the id
  ##------------------------------------------------------
  user.option <- list(...)
  if (any(is.na(id)) || any(is.na(y)) || any(is.na(x)) || any(is.na(tm))) {
    stop("missing values encountered: remove observations with missing values")
  }
  x <- as.data.frame(x)
  p <- ncol(x)
  xvar.names <- colnames(x)
  X <- do.call(rbind, lapply(1:n, function(i) {
    x[id == id.unq[i],, drop = FALSE][1,, drop = FALSE]}))
  x <- do.call(rbind, lapply(1:n, function(i) {x[id == id.unq[i],, drop = FALSE]}))
  Ymean <- mean(y, na.rm = TRUE)
  Ysd <- sd(y, na.rm = TRUE)
  if (Ysd < 1e-6) {
    Ysd <- 1
  }
  Yorg <- lapply(1:n, function(i) {y[id == id.unq[i]]})
  Y <- lapply(1:n, function(i) {(y[id == id.unq[i]] - Ymean) / Ysd})
  ni <- unlist(lapply(1:n, function(i) {sum(id == id.unq[i])}))
  id <- sort(id)
  
  ##------------------------------------------------------
  ## define the Di matrix
  ##------------------------------------------------------

  ## define the unique time points
  ## we allow replicated time measurements for an individual
  ## therefore the match of tm.i to tm.unq is delicate
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
  if (d >= 1) {
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
    cbind(X.tm[tm.id[[i]],, drop = FALSE])
  })

  ##------------------------------------------------------
  ## regularization details
  ## allows different nu values for b0 and b1
  ##------------------------------------------------------
  nu <- {if (length(nu) > 1) nu else rep(nu, 2)}
  if (sum(!(0 < nu & nu <= 1)) > 0) {
    stop("regularization parameter (nu) must be in (0,1]")
  }
  nu.vec <- c(nu[1], rep(nu[2], df.D - 1))

  ##------------------------------------------------------
  ## tree/forest base-learner parameters used for rfsrc
  ##------------------------------------------------------
  ntree <- is.hidden.ntree(user.option)
  bootstrap <- is.hidden.bootstrap(user.option)
  if (ntree == 1) {
    nodesize <- max(1, round(n/(2 * K)))
    mtry <- df.D + p
  }
  else {
    nodedepth <- max(0, log(max(0, K - 1), base = 2))
    nodesize <- 1
    mtry <- NULL
    if (missing(lambda) || lambda < 0) {
      lambda <- 0
    }
  }
  
  ##------------------------------------------------------
  ## error.rate/vimp details
  ## TBD TBD TBD
  ##------------------------------------------------------
  importance <- FALSE
  vimpFlag <- bootstrap == "by.root" && importance && ntree == 1
  vimp <- NULL
    

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
  if (!univariate && ntree == 1 && (missing(lambda) || lambda < 0)) {
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


  ##------------------------------------------------------
  ## initialization
  ## set various dimensions/terms that will be required
  ## final parameter checks
  ## multivariate tree formula/details
  ##------------------------------------------------------
  mu <- lapply(1:n, function(i) {rep(0, ni[i])})
  beta <- matrix(0, n, df.D)
  if (vimpFlag) {
    vimp <- matrix(0, M, p)
    colnames(vimp) <- xvar.names
  }
  if (ntree == 1) {
    baselearner <- membership.list <- gamma.list <- vector("list", length = M)
  }
  else {
   membership.list <- gamma.list <- NULL
   baselearner <- vector("list", length = M)
  }
  lambda.vec <- phi.vec <- rho.vec <- rep(0, M)
  lambda.initial <- var(unlist(Y), na.rm = TRUE)

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

  ## initialize sigma (=lambda * (1-rho)) and phi
  sigma <- phi <- 1
  if (!lambda.est.flag) {
    sigma <- lambda * (1 - rho)
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

  for (m in 1:M) {

    if (verbose) {
      cat("\t-- iteration:", m, "\n")
    }

    ##---------------------------------------------------------
    ## step 3: calculate the negative gradient
    ## we discard the unnecessary (1-rho)^{-1} constant
    ##---------------------------------------------------------
    gm <- t(matrix(unlist(lapply(1:n, function(i) {
      rmi <- rho.inv(ni[i], rho)##this function controls instability in R^{-1}
      cmi <- rmi * sum(Y[[i]] - mu[[i]], na.rm = TRUE)
      t(D[[i]]) %*% (Y[[i]] - mu[[i]] - cbind(rep(cmi, ni[i])))
    })), nrow = df.D))

    ##---------------------------------------------------------
    ##step 4(a): fit a K-terminal node multivariate regression tree
    ##
    ## gm    n x df.D  gradient-matrix
    ## X     n x p data frame containing x values
    ## K     desired number of terminal nodes
    ##
    ## Currently pruning is only available in predict mode.
    ## To recover the pruned membership of the training data,
    ## we must use the option "ptn.count"
    ##
    ## Thus we first we grow a tree with > K terminal nodes.  We then
    ## prune the tree back to K terminal nodes.  For convenience, we
    ## roughly calculate the initial nodesize such that the number of
    ## terminal nodes is 2 * K.
    ##
    ##
    ## If ntree > 1, we don't worry about pruning and grow a forest with
    ## roughly K terminal nodes.  We do this by constraining the node depth
    ## to log_2(K-1) and use a minimal nodesize of 1.
    ##
    ##---------------------------------------------------------
    incoming.data <- cbind(gm, X)
    names(incoming.data) = c(Y.names, names(X))

    ## multivariate forest learner
    if (ntree > 1) {

      rfsrc.obj <- rfsrc(rfsrc.f,
                         data = incoming.data,
                         mtry = mtry,
                         nodedepth = nodedepth,
                         nodesize = nodesize,
                         importance = "none",
                         bootstrap = bootstrap,
                         ntree = ntree,
                         forest.wt = TRUE)##TBD TBD, use OOB forest weights?

      ## Kmax is the maximum number of terminal nodes
      Kmax <- max(rfsrc.obj$leaf.count, na.rm = TRUE)
      
      ## save the base learner
      baselearner[[m]] <- list(forest = rfsrc.obj)


    }

    ## multivariate tree base learner 
    else {
      
      rfsrc.obj <- rfsrc(rfsrc.f,
                         data = incoming.data,
                         ntree = 1,
                         mtry = mtry,
                         nodesize = nodesize,
                         importance = "none",
                         bootstrap = bootstrap)

      ## save the base learner
      baselearner[[m]] <- rfsrc.obj
      
      ## the outcome.target is irrelevant
      result.pred <- predict.rfsrc(rfsrc.obj,
                                   membership = TRUE,
                                   ptn.count = K,
                                   importance = "none")


      ## The pruned membership is an n x 1 vector of non-consecutive
      ## integers.  They map immutably to NODE_ID of the tree.
      membership <- membership.org <- c(result.pred$ptn.membership)
      membership.list[[m]] <- membership.org
      
      ## Recode the membership to be sequential
      membership <- as.numeric(factor(membership))
      ptn.id <- unique(membership)
      ## Kmax is the number of pseudo-terminal nodes, where Kmax <= K
      Kmax <-  length(ptn.id)

      ## membership for noised up data (if vimp requested)
      if (vimpFlag) {

        ## record the OOB data
        ## NOTE: OOB is subject/id specific, exactly what we want in longitudinal settings
        oob <- which(rfsrc.obj$inbag == 0)
        n.oob <- length(oob)
        
        ## construct the noise matrix: restricted to OOB data 
        Xnoise <- do.call(rbind, lapply(1:p, function(k) {
          X.k <- X[oob,, drop = FALSE]
          X.k[, k] <- sample(X.k[, k])
          X.k
        }))

        ##determine the noised up membership
        membershipNoise <- c(predict.rfsrc(rfsrc.obj,
                        newdata = Xnoise,
                        membership = TRUE,
                        ptn.count = K,
                        importance = "none")$ptn.membership)

        ## Recode the membership to be sequential: needs to match the full sample order
        membershipNoise <- as.numeric(factor(membershipNoise, levels = levels(factor(membership.org))))
        
      }

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
            ci <- rho.inv.sqrt(ni[i], rho)##this function controls instability in R^{-1/2}
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
        ## iterate to estimate lambda 
        ## ----------------------------------------------------------------------

        ## initialize lambda.hat
        lambda.hat <- lambda.initial

        ## iterative loop
        for (k in 1:lambda.iter) {

          ## BLUP for random effects (u_k) and fixed effects (alpha_k) conditional on lambda
          blup.obj <-  blup.solve(transf.data, membership, lambda.hat, Kmax)

          ## lambda method of moments estimator conditional on BLUP
          ## robust rss calculation used to avoid deflating lambda estimate
          lambda.obj <- lapply(1:Kmax, function(k) {
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
          num <- sum(unlist(lapply(1:Kmax, function(k) {lambda.obj[[k]]$trace.Z})), na.rm = TRUE)
          den <- sum(unlist(lapply(1:Kmax, function(k) {lambda.obj[[k]]$rss})), na.rm = TRUE)
          N <- sum(unlist(lapply(1:Kmax, function(k) {lambda.obj[[k]]$resid})), na.rm = TRUE)

          if (!is.na(den) && den > (.99 * N)) {## ensure that denominator is larger than degrees of freedom
            lambda.hat <- num / (den - .99 * N)
          }
          else {
            lambda.hat <- min(lambda.hat, lambda.max)
          }

          ## cap lambda.hat
          lambda.hat <- min(lambda.hat, lambda.max)

        }

        ## update lambda, sigma
        lambda <- lambda.hat 
        sigma <- lambda * (1 - rho) 

      }

      ##---------------------------------------------------------
      ## LEAST SQUARES SOLUTION FOR GAMMA: ntree = 1
      ## step 4(b)
      ##
      ## Xnew      WLS x-matrix
      ## YnewSum   summed pseudo y-value
      ## XnewSum   summed pseudo WLS x-values
      ## gamma     weighted least squares solution
      ##---------------------------------------------------------

      Xnew <- mclapply(1:n, function(i) {
        rmi <- rho.inv(ni[i], rho)##this function controls instability in R^{-1}
        Wi <- diag(1, ni[i]) - matrix(rmi, ni[i], ni[i])
        t(D[[i]]) %*% Wi %*% D[[i]]
      })
      gamma <- lapply(1:Kmax, function(k) {
        pt.k <- (membership == k)
        ## the following is only relevant if pt.k is non-NULL (can happen in cv setting)
        if (sum(pt.k) > 0) {
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
        }
        else {## NULL case
          rep(0, df.D)
        }
      })

      ##---------------------------------------------------------
      ## step 4(c): save gamma: ntree = 1
      ## needed for prediction
      ## save as a matrix and include original grow terminal node membership
      ##---------------------------------------------------------
      gamma.matx <- matrix(0, Kmax, df.D+1)
      gamma.matx[, 1] <- sort(unique(membership.org))
      gamma.matx[, 2:(df.D+1)] <- matrix(unlist(gamma), ncol = df.D, byrow = TRUE)

      gamma.list[[m]] <- gamma.matx

      ##---------------------------------------------------------
      ## step 5: update beta and mu: ntree = 1
      ##---------------------------------------------------------

      ## beta update
      beta.update <- matrix(unlist(lapply(1:n, function(i) {
        gamma[[membership[i]]]})), nrow = df.D)
      beta.old <- beta
      beta <- beta.old + t(beta.update * nu.vec)
      ## mu update
      mu.old <- mu
      mu <- lapply(1:n, function(i) {D[[i]] %*% beta[i, ]})


      ##---------------------------------------------------------
      ## step 5': vimp
      ## TBD TBD TBD: currently does not work properly
      ##---------------------------------------------------------      
      if (vimpFlag) {

        ## standardize the current beta coefficient: restricted to OOB data
        ## calculate the tree contributed oob error
        beta.std <- beta[oob,, drop = FALSE] * Ysd
        beta.std[, 1] <- beta.std[, 1] + Ymean
        err.oob <- l2Dist(Yorg[oob], lapply(1:n.oob, function(i) {D[[oob[i]]] %*% beta.std[i, ]}))

        ## the vimp oob requires an updated vimped beta
        ## calculate the tree contributed noised up error
        beta.update.vimp <- lapply(1:p, function(k) {
          membership.k <- membershipNoise[((k-1) * n.oob + 1):(k * n.oob)]
          matrix(unlist(lapply(1:n.oob, function(i) {
            gamma[[membership.k[i]]]})), nrow = df.D)
        })
        err.vimp <- sapply(1:p, function(k) {
          beta.vimp.k <- beta.old[oob,, drop = FALSE] + t(beta.update.vimp[[k]] * nu.vec)
          beta.vimp.k <- beta.vimp.k * Ysd
          beta.vimp.k[, 1] <- beta.vimp.k[, 1] + Ymean
          l2Dist(Yorg[oob], lapply(1:n.oob, function(i) {D[[oob[i]]] %*% beta.vimp.k[i, ]}))
        })

        ## save the *standardized* vimp
        vimp[m, ] <- 100 * (err.vimp - err.oob) / Ysd
        
      }

    }
   

    ##############################################################################
    ##
    ##
    else{#####FOREST BASE LEARNER ######
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
        rmi <- rho.inv(ni[i], rho)##this function controls instability in R^{-1}
        Wi <- diag(1, ni[i]) - matrix(rmi, ni[i], ni[i])
        t(D[[i]]) %*% Wi %*% D[[i]]
      })

      ## iterate over cases i to get bhat
      ## for speed we eliminate cases with small forest weights

      bhat <- do.call("cbind", mclapply(1:n, function(i) {
        fwt.i <- forest.wt[i, ]
        fwt.i[fwt.i <= forest.tol] <- 0
        pt.i <- (fwt.i != 0)
        ## the following is only relevant if pt.i is non-NULL
        if (sum(pt.i) > 0) {
          fwt.i <- fwt.i / sum(fwt.i)
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
        }
        else {##NULL case
          rep(0, df.D)
        }
      }))

      ##---------------------------------------------------------
      ## update beta and  mu
      ##---------------------------------------------------------
      beta <- beta + t(bhat * nu.vec)
      mu <- lapply(1:n, function(i) {D[[i]] %*% beta[i, ]})

      ##--------------------------------------------------------------------
      ## baselearner needs certain objects in predict mode: append them here
      ##--------------------------------------------------------------------
      baselearner[[m]] <- c(baselearner[[m]],
                            list(gm = gm),
                            list(Xnew = Xnew),
                            list(pen = sigma * pen.lsq.matx))


    }

    ##############################################################################
    ##
    ##  The algorithm is now the same
    ##
    #############################################################################

    if (!univariate) {

      ##---------------------------------------------------------
      ## step 6 and 7: update phi and rho using REML
      ## we make a call to lme (actually gls since there are no r.eff)
      ## we use a compound symmetric correlation matrix, but this can be generalized
      ## note: tm and x are used in the REML call -- although the theory does not call for it
      ##---------------------------------------------------------
      resid.data <- data.frame(y  = unlist(lapply(1:n, function(i) {Y[[i]] - mu[[i]]})),
                               x,
                               tm = unlist(lapply(1:n, function(i) {tm[id == id.unq[i]]})),
                               id = unlist(lapply(1:n, function(i) {rep(id.unq[i], ni[i])})))
      gls.obj <- tryCatch({gls(y ~ ., data = resid.data,
                               correlation = corCompSymm(form = ~ 1 | id))},
                          error = function(ex){NULL})
      if (is.null(gls.obj)) {
        gls.obj <- tryCatch({gls(y ~ 1, data = resid.data,
                                 correlation = corCompSymm(form = ~ 1 | id))},
                            error = function(ex){NULL})
      }
      if (!is.null(gls.obj)) {
        phi <- gls.obj$sigma^2
        if (rho.fit.flag) {
          rho <- as.numeric(coef(gls.obj$modelStruct$corStruc, unconstrained = FALSE))
          rho <- max(min(0.999, rho, na.rm = TRUE), -0.999)
        }
      }
      
      phi.vec[m] <- phi * Ysd ^ 2##scale by the overall variance
      rho.vec[m] <- rho
      
      if (verbose) {
        cat("phi   :", phi.vec[m], "\n")
        cat("rho   :", rho.vec[m], "\n")
      }
    
    
      ##---------------------------------------------------------
      ## step 8: update sigma, save lambda
      ##---------------------------------------------------------
      
      sigma <- lambda * (1 - rho) 
      lambda.vec[m] <- lambda
      if (verbose) {
        cat("lambda:", lambda.vec[m], "\n")
      }

    }

  }##end boosting iteration

  ##---------------------------------------------------------
  ## rescale by the std and add the Y mean back
  ##---------------------------------------------------------
  beta <- beta * Ysd
  beta[, 1] <- beta[, 1] + Ymean
  mu <- lapply(1:n, function(i) {c(mu[[i]] * Ysd + Ymean)})
  y <- lapply(1:n, function(i) {y[id == id.unq[i]]})

  ##---------------------------------------------------------
  ## return the promised object
  ##---------------------------------------------------------

  obj <- list(x = X,
              xvar.names = xvar.names,
              time = lapply(1:n, function(i) {tm[id == id.unq[i]]}),
              id = id,
              y = Yorg,
              ymean = Ymean,
              ysd = Ysd,
              gamma = gamma.list,
              beta = beta,
              mu = mu,
              lambda = lambda.vec,
              phi = phi.vec,
              rho = rho.vec,
              baselearner = baselearner,
              membership = membership.list,
              vimp = if (!is.null(vimp)) colSums(vimp, na.rm = TRUE) else NULL,
              D = X.tm,
              d = d,
              pen.ord = pen.ord,
              K = K,
              M = M,
              nu = nu,
              ntree = ntree)

  class(obj) <- c("boostmtree", "grow", learnerUsed)

  invisible(obj)

}

