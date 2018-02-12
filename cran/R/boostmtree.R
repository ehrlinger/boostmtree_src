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
                       svd.tol = 1e-6,
                       forest.tol = 1e-3,
                       verbose = TRUE,
                       cv.flag = FALSE,
                       eps = 1e-5,
                       importance = FALSE,
                       mod.grad = TRUE,
                       ...)
{
  
  ##------------------------------------------------------
  ## custom mclapply/lapply switch
  ##------------------------------------------------------
  if (grepl("Debian", Sys.info()["version"])) {
    papply <- lapply
  }
  else {
    papply <- mclapply
  }

  ##------------------------------------------------------
  ## is this a univariate setting?
  ## if so, flag it and make a zero time vector
  ## we also set d < 0 to trigger no time-covariate fitting
  ##------------------------------------------------------
  univariate <- FALSE
  if (missing(tm)) {
    id <- 1:nrow(x)
  }
  id.unq <- sort(unique(id))
  n <- length(id.unq)
  if (length(id.unq) == length(id)) {
    univariate <- TRUE
    tm <- rep(0, n)
    d <- -1
  }
  if (univariate) {
    mod.grad <- FALSE
    rho <- 0
    lambda.vec <- phi.vec <- rho.vec <- NULL
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
    nodedepth <- max(0, log(max(0, K), base = 2))
    nodesize <- 1
    mtry <- NULL
    if (missing(lambda) || lambda < 0) {
      lambda <- 0
    }
  }

  ##------------------------------------------------------
  ## define the learner (used for setting the class)
  ##------------------------------------------------------
  if (ntree > 1) {
    if (univariate) {
      learnerUsed <- "forest learner"
    }
    else {
      learnerUsed <- "mforest learner"
    }
  }
  else {
    if (df.D == 1) {
      learnerUsed <- "tree learner"
    }
    else {
      learnerUsed <- "mtree-Pspline learner"
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
      ## use svd to get the square root and inverse square root of penalty
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
  if (ntree == 1) {
    baselearner <- membership.list <- gamma.list <- vector("list", length = M)
  }
  else {
   membership.list <- gamma.list <- NULL
   baselearner <- vector("list", length = M)
  }
  if (!univariate) {
    lambda.vec <- phi.vec <- rho.vec <- rep(0, M)
  }
  lambda.initial <- var(unlist(Y), na.rm = TRUE)

  ## rho initialization
  rho.fit.flag <- TRUE
  rho.tree.grad <- 0
  rho.hide <- is.hidden.rho(user.option)
  ## did user provide a fixed rho value?
  if (!is.null(rho.hide) && (rho.hide >= 0 && rho.hide < 1)) {
    rho.fit.flag <- FALSE
    rho <- rho.hide
  }
  else {
    rho <- 0
  }
  
  ## initialize sigma and phi
  ## we use a robust version of sigma for penalization
  sigma <- phi <- 1
  if (!lambda.est.flag) {
    sigma <- sigma.robust(lambda, rho)
  }

  ## formula
  Y.names <- paste("Y", 1:df.D, sep = "")
  rfsrc.f <- as.formula(paste("Multivar(", paste(Y.names, collapse = ","), paste(") ~ ."), sep = ""))

  ##------------------------------------------------------
  ## cross-validation details
  ## only applies if ntree = 1
  ## parse for hidden options
  ## we allow various combinations of CV lambda,rho
  ## vimp not currently implemented
  ##------------------------------------------------------
  cv.flag <- cv.flag && (ntree == 1)
  cv.lambda.flag <- cv.flag && is.hidden.CVlambda(user.option) && lambda.est.flag
  cv.rho.flag <- cv.flag && is.hidden.CVrho(user.option) && rho.fit.flag
  vimp.flag <- importance && cv.flag
  vimp.flag <- FALSE##TBD TBD TBD implement vimp

  if (cv.flag) {
    ## assign lists/vectors
    mu.cv.list <- vector("list", M)
    mu.cv <- lapply(1:n, function(i) {rep(0, ni[i])})
    mu.i <- lapply(1:n, function(i) {
      lapply(1:n, function(j) {rep(0, ni[j])})
    })
    err.rate <- matrix(NA, M, 2)
    colnames(err.rate) <- c("l1", "l2")
    ## hold out mean and std
    Ymean.i <- sapply(1:n, function(i) {
      mean(unlist(Yorg[-i]), na.rm = TRUE) 
    })
    Ysd.i <- sapply(1:n, function(i) {
      sd.i <- sd(unlist(Yorg[-i]), na.rm = TRUE)
      if (sd.i < 1e-6) {
        1
      }
      else {
        sd.i
      }
    })
    ## vimp matrix
    if (vimp.flag) {
      vimp <- matrix(0, M, p)
      colnames(vimp) <- xvar.names
    }
    else {
      vimp <- NULL
    }
  }
  else {
    err.rate <- rmse <- Mopt <- vimp <- NULL
  }


  
  ##------------------------------------------------------
  ## MAIN LOOP
  ##
  ## gradient boosting loop
  ##
  ## MAIN LOOP
  ##------------------------------------------------------
  if (verbose) pb <- txtProgressBar(min = 0, max = M, style = 3)

  for (m in 1:M) {

    if (verbose) setTxtProgressBar(pb, m)
    if (verbose && m == M) cat("\n")


    ##---------------------------------------------------------
    ## step 3: calculate the negative gradient
    ## we discard the unnecessary (1-rho)^{-1} constant
    ##
    ## NEW ADDITION
    ## do we use a modified gradient for the tree growing? 
    ## if so, replace rho with rho.tree.grad=0
    ##---------------------------------------------------------
    if (mod.grad == FALSE) {
      gm <- gm.mod <- t(matrix(unlist(lapply(1:n, function(i) {
        rmi <- rho.inv(ni[i], rho)##this function controls instability in R^{-1}
        cmi <- rmi * sum(Y[[i]] - mu[[i]], na.rm = TRUE)
        t(D[[i]]) %*% (Y[[i]] - mu[[i]] - cbind(rep(cmi, ni[i])))
      })), nrow = df.D))
    }
    else {
      gm.mod <- t(matrix(unlist(lapply(1:n, function(i) {
        rmi <- rho.inv(ni[i], rho.tree.grad)##this function controls instability in R^{-1}
        cmi <- rmi * sum(Y[[i]] - mu[[i]], na.rm = TRUE)
        t(D[[i]]) %*% (Y[[i]] - mu[[i]] - cbind(rep(cmi, ni[i])))
      })), nrow = df.D))
      gm <- t(matrix(unlist(lapply(1:n, function(i) {
        rmi <- rho.inv(ni[i], rho)##this function controls instability in R^{-1}
        cmi <- rmi * sum(Y[[i]] - mu[[i]], na.rm = TRUE)
        t(D[[i]]) %*% (Y[[i]] - mu[[i]] - cbind(rep(cmi, ni[i])))
      })), nrow = df.D))
    }
    
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
    ## to log_2(K) and use a minimal nodesize of 1.
    ##
    ##---------------------------------------------------------
    incoming.data <- cbind(gm.mod, X)
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
                         forest.wt = TRUE, 
                         memebership = TRUE)

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
                         bootstrap = bootstrap,
                         membership = TRUE)

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

      ## vimp: membership for noised up data 
      if (vimp.flag) {

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

        transf.data <- papply(1:n, function(i) {
          if (ni[i] > 1) {
            ci <- rho.inv.sqrt(ni[i], rho)##this function controls instability in R^{-1/2}
            R.inv.sqrt <- (diag(1, ni[i]) - matrix(ci, ni[i], ni[i])) / sqrt(1 - rho)
          }
          else {
            R.inv.sqrt <- cbind(1)
          }
          if (cv.lambda.flag) {
            Ynew <- R.inv.sqrt %*% (Y[[i]] - mu.cv[[i]])
          }
          else {
            Ynew <- R.inv.sqrt %*% (Y[[i]] - mu[[i]])
          }
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

        ## update lambda, sigma (robust version)
        lambda <- lambda.hat 
        sigma <- sigma.robust(lambda, rho) 

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

      Xnew <- papply(1:n, function(i) {
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
      gamma.matx <- matrix(0, Kmax, df.D + 1)
      gamma.matx[, 1] <- sort(unique(membership.org))
      gamma.matx[, 2:(df.D+1)] <- matrix(unlist(gamma), ncol = df.D, byrow = TRUE)

      gamma.list[[m]] <- gamma.matx

      ##---------------------------------------------------------
      ## step 5: update mu
      ##---------------------------------------------------------

      ## mu update
      bhat <- t(matrix(unlist(lapply(1:n, function(i) {
        gamma[[membership[i]]]})), nrow = df.D) * nu.vec)
      mu <- lapply(1:n, function(i) {mu[[i]] + D[[i]] %*% bhat[i, ]})


      ##---------------------------------------------------------
      ## step 6: in-sample cv
      ## includes err.rate caclulations
      ## includes vimp calculations (if requested)
      ##---------------------------------------------------------
      if (cv.flag) {

        ## iterate over each case, holding it out
        mu.i <- lapply(1:n,function(i) {

          mem.i <- membership[i]  
          mu.ij <- mu.i[[i]]

          ## hold out gradient
          grad.i <- t(matrix(unlist(lapply(1:n, function(i) {
            rmi <- rho.inv(ni[i], rho)##this function controls instability in R^{-1}
            cmi <- rmi * sum(Y[[i]] - mu.ij[[i]], na.rm = TRUE)
            t(D[[i]]) %*% (Y[[i]] - mu.ij[[i]] - cbind(rep(cmi, ni[i])))
          })), nrow = df.D))

          ## hold out gamma
          gamma.i <- lapply(1:Kmax, function(k) {
            pt.k <- (membership == k)
            YnewSum <- colSums(grad.i[pt.k, , drop = FALSE])
            XnewSum <- Reduce("+", lapply(which(pt.k), function(j) {Xnew[[j]]}))
            if (is.null(XnewSum)) {
              XnewSum <- matrix(0,df.D,df.D)
            }
            else {
              XnewSum <- XnewSum
            }
            if (k == mem.i){
              XnewSum <- XnewSum - Xnew[[i]]
              YnewSum <- YnewSum - grad.i[i, ]
            }
            else {
              XnewSum <- XnewSum
              YnewSum <- YnewSum
            }
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

          ## save the hold out gamma value as a matrix
          gamma.matx.i <- matrix(0, Kmax, df.D + 1)
          gamma.matx.i[, 1] <- 1:Kmax
          gamma.matx.i[, 2:(df.D+1)] <- matrix(unlist(gamma.i), ncol = df.D, byrow = TRUE)

          ## return the hold out mu
          lapply(1:n,function(j) {
            which.j <- which(gamma.matx.i[, 1] == membership[j])
            mu.i[[i]][[j]] + c(D[[j]] %*% (gamma.matx.i[which.j, -1] * nu.vec))
          })

        })

        ## update the hold out mu
        mu.cv <- lapply(1:n,function(i){mu.i[[i]][[i]]})
        mu.cv.list[[m]] <- mu.cv

        ## update the scaled-centered hold out mu
        mu.cv.org <- lapply(1:n,function(i){mu.cv[[i]] * Ysd.i[i] + Ymean.i[i]})
        err.rate[m, ] <- c(l1Dist(Yorg, mu.cv.org), l2Dist(Yorg, mu.cv.org))
 
        ## TBD TBD TBD
        ## vimp (if requested)
        ## need to save perturbed hold out mu: seems very time consuming
        if (vimp.flag) {
        
        }

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
      Xnew <- papply(1:n, function(i) {
        rmi <- rho.inv(ni[i], rho)##this function controls instability in R^{-1}
        Wi <- diag(1, ni[i]) - matrix(rmi, ni[i], ni[i])
        t(D[[i]]) %*% Wi %*% D[[i]]
      })

      ## iterate over cases i to get bhat
      ## for speed we eliminate cases with small forest weights

      bhat <- do.call("cbind", papply(1:n, function(i) {
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
      ## update  mu
      ##---------------------------------------------------------
      bhat <- t(bhat * nu.vec)
      mu <- lapply(1:n, function(i) {mu[[i]] + D[[i]] %*% bhat[i, ]})

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
    if (!univariate && rho.fit.flag) {

      ##---------------------------------------------------------
      ## update phi and rho using REML
      ## we make a call to lme (actually gls since there are no r.eff)
      ## we use a compound symmetric correlation matrix, but this can be generalized
      ## note: tm and x are used in the REML call -- although the theory does not call for it
      ##---------------------------------------------------------
      if (cv.rho.flag) {
        resid.data <- data.frame(y  = unlist(lapply(1:n, function(i) {Y[[i]] - mu.cv[[i]]})),
                                 x,
                                 tm = unlist(lapply(1:n, function(i) {tm[id == id.unq[i]]})),
                                 id = unlist(lapply(1:n, function(i) {rep(id.unq[i], ni[i])})))
      }
      else {
        resid.data <- data.frame(y  = unlist(lapply(1:n, function(i) {Y[[i]] - mu[[i]]})),
                                 x,
                                 tm = unlist(lapply(1:n, function(i) {tm[id == id.unq[i]]})),
                                 id = unlist(lapply(1:n, function(i) {rep(id.unq[i], ni[i])})))
      }
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
        rho <- as.numeric(coef(gls.obj$modelStruct$corStruc, unconstrained = FALSE))
        rho <- max(min(0.999, rho, na.rm = TRUE), -0.999)
      }

    }

    ## save rho/phi values
    if (!univariate) {
      phi.vec[m] <- phi * Ysd ^ 2##scale by the overall variance
      rho.vec[m] <- rho
    }
    
    ##---------------------------------------------------------
    ## update sigma, save lambda
    ## use robust sigma function to enforce numerical stability
    ## when penalizing
    ##---------------------------------------------------------
    if (!univariate) {
      sigma <- sigma.robust(lambda, rho)
      lambda.vec[m] <- lambda
    }

  }##end boosting iteration

  ##---------------------------------------------------------
  ## rescale by the std and add the Y mean back
  ##---------------------------------------------------------
  mu <- lapply(1:n, function(i) {c(mu[[i]] * Ysd + Ymean)})
  y <- lapply(1:n, function(i) {y[id == id.unq[i]]})

  ##---------------------------------------------------------
  ## final cv details
  ##---------------------------------------------------------
  if (cv.flag) {

    ## determine the in-sample optimized boosting iterations
    diff.err <- abs(err.rate[, "l2"] - min(err.rate[, "l2"], na.rm = TRUE))
    diff.err[is.na(diff.err)] <- 1
      if (sum(diff.err < Ysd * eps) > 0) {
        Mopt <- min(which(diff.err < eps))
      }
      else {
        Mopt <- M
      }

    ## what is the rmse at Mopt?
    rmse <- err.rate[Mopt, "l2"]

    ## pull the cross-validated mu at Mopt
    mu <- lapply(1:n,function(i){mu.cv.list[[Mopt]][[i]] * Ysd.i[i] + Ymean.i[i]})
    
  }



  
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
              mu = mu,
              lambda = lambda.vec,
              phi = phi.vec,
              rho = rho.vec,
              baselearner = baselearner,
              membership = membership.list,
              D = X.tm,
              d = d,
              pen.ord = pen.ord,
              K = K,
              M = M,
              nu = nu,
              ntree = ntree,
              err.rate = if (!is.null(err.rate)) err.rate / Ysd else NULL,
              rmse = if (!is.null(rmse)) rmse / Ysd else NULL,
              Mopt = Mopt,
              vimp = if (!is.null(vimp)) colSums(vimp, na.rm = TRUE) else NULL,
              forest.tol = forest.tol)

  class(obj) <- c("boostmtree", "grow", learnerUsed)

  invisible(obj)

}

