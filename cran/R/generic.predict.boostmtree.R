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
                                       eps = 1e-5,
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
      tm.unq <- sort(unique(unlist(object$time)))
      testFlag <- FALSE
    }
    else {
      ## if x is provided but no other values
      ## we assume id's are 1:1 and swap in the original unique time values
      if (!missing(x) && (missing(id) || missing(tm))) {
        X <- x
        n <- nrow(X)
        tm.unq <- sort(unique(unlist(object$time)))
        tm <- lapply(1:n, function(i){tm.unq})
        id <- id.unq <- 1:n
        if (!missing(y)) {
          Y <- lapply(1:n, function(i) {y[i]})
          testFlag <- TRUE
        }
        else {
          testFlag <- FALSE
        }
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
  beta <- matrix(0, n, df.D)
  if (ntree == 1) {
    beta.vimp <- beta.cov.vimp <- beta.time.vimp <- vector("list", p)
  }

  ## initialize mu
  mu.list <- vector("list", M)

  ## forest.tolerance
  forest.tol <- object$forest.tol
  
  ## vimp details
  ## TBD TBD TBD VIMP NOT IMPLEMENTED FOR FORESTS TBD TBD TBD
  vimpFlag <- testFlag && importance && ntree == 1
  vimp <- NULL

  ##############################################################################
  ##
  ##  Iterate over the M iterations to recursively define the test predictor
  ##  At this point the algorithm diverges depending upon ntree
  ##
  if (ntree == 1) {#####SINGLE TREE BASE LEARNER#####
  ##
  ##
  #############################################################################

    ## we make the predict call in parallel in R for speed
    ## to execute a parallel call we must over-ride rf.cores
    ## we must be careful not to over-tax mc.cores
    rf.cores.old <- getOption("rf.cores")
    mc.cores.old <- getOption("mc.cores")
    
    ##---------------------------------------------------------
    ## membership
    ## pass the x-data down the mth tree
    ## acquire terminal node membership (tnm)
    ##---------------------------------------------------------
    membership <- papply(1:M, function(m) {

      ## verbose output
      ##if (verbose) {
      ##cat("\t-- iteration:", m, "\n")
      ##}

      ## set the rf.cores/mc.cores to minimal values
      options(rf.cores = 1, mc.cores = 1)
      c(predict.rfsrc(baselearner[[m]],
                      newdata = X,
                      membership = TRUE,
                      ptn.count = K,
                      importance = "none")$ptn.membership)

    })


    ##---------------------------------------------------------
    ## obtain mu by determining gamma making use of tnm
    ##---------------------------------------------------------
    nullObj <- lapply(1:M, function(m) {
      ## obtain updated beta
      orgMembership <- gamma[[m]][, 1]
      beta.m <- t(gamma[[m]][match(membership[[m]], orgMembership), -1, drop = FALSE]) * nu.vec
      ## extract those mu values corresponding to the observed time points
      ## we allow replicated time measurements for an individual
      ## therefore the match of tm to tm.unq is delicate
      Dbeta.m <- D %*% beta.m
      ## update mu
      if (m == 1) {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      else {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          unlist(mu.list[[m-1]][i]) + Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      NULL##memory saving measure
    })
    rm(nullObj)

    ##---------------------------------------------------------
    ## scale mu and add Ymean
    ##---------------------------------------------------------
    mu.list <- lapply(mu.list, function(mlist){  
      lapply(1:n,function(i) {mlist[[i]] * Ysd + Ymean})
    })
    
    ##---------------------------------------------------------
    ## determine a cumulative error rate
    ##---------------------------------------------------------
    if (testFlag) {
      err.rate <- matrix(unlist(lapply(mu.list, function(mlist) {
        c(l1Dist(Y, mlist), l2Dist(Y, mlist)) 
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
      if (sum(diff.err < Ysd * eps) > 0) {
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
    ##
    ## VIMP
    ## computational saving as we iterate only from 1 to Mopt
    ##
    ##---------------------------------------------------------
    if (vimpFlag) {
      membershipNoise <- papply(1:Mopt, function(m) {

        Xnoise <- do.call(rbind, lapply(1:p, function(k) {
          X.k <- X
          X.k[, k] <- sample(X.k[, k])
          X.k
        }))

        ## set the rf.cores/mc.cores to minimal values
        options(rf.cores = 1, mc.cores = 1)
        c(predict.rfsrc(baselearner[[m]],
                        newdata = Xnoise,
                        membership = TRUE,
                        ptn.count = K,
                        importance = "none")$ptn.membership)
      })

      ## reset the rf.cores/mc.cores option to its original value
      if (!is.null(rf.cores.old)) options(rf.cores = rf.cores.old)
      if (!is.null(mc.cores.old)) options(mc.cores = mc.cores.old)
      
      nullObj <- lapply(1:Mopt, function(m) {

        ##---------------------------------------------------------
        ## the unperturbed beta
        ##---------------------------------------------------------
        orgMembership <- gamma[[m]][, 1]
        beta.m <- t(t(gamma[[m]][match(membership[[m]], orgMembership), -1, drop = FALSE]) * nu.vec) * Ysd

        if (m == 1) {
          beta.m[, 1] <- beta.m[, 1] + Ymean
          beta <<- beta.m
        }
        else {
          beta <<- beta + beta.m
        }

        ##---------------------------------------------------------
        ## vimp calculation using the perturbed X
        ##---------------------------------------------------------
        beta.vimp <<- lapply(1:p, function(k) {
          membership.k <- membershipNoise[[m]][((k-1) * n + 1):(k * n)]
          beta.m.k <- t(t(gamma[[m]][match(membership.k, orgMembership), -1, drop = FALSE]) * nu.vec) * Ysd
          if (m == 1) {
            beta.m.k[, 1] <- beta.m.k[, 1] + Ymean##add the overall mean
            beta.m.k
          }
          else {
            beta.vimp[[k]] + beta.m.k
          }
        })
        ## break vimp into covariate only, covariate-time effects
        ## applies only if splines are fit
        if (df.D > 1) {
          beta.cov.vimp <<- lapply(1:p, function(k) {
            b.c.v <- beta.vimp[[k]]
            b.c.v[, -1] <- beta[, -1]
            b.c.v
          })
          beta.time.vimp <<- lapply(1:p, function(k) {
            b.t.v <- beta.vimp[[k]]
            b.t.v[, 1] <- beta[, 1]
            b.t.v
          })
        }
        
        NULL##memory saving measure
        
      })##loop is complete
      rm(nullObj)
    
    }## end VIMP end VIMP end VIMP end

    ##---------------------------------------------------------
    ##
    ## no vimp requested
    ## just calculate beta
    ##
    ##---------------------------------------------------------

    else {

      nullObj <- lapply(1:Mopt, function(m) {

        orgMembership <- gamma[[m]][, 1]
        beta.m <- t(t(gamma[[m]][match(membership[[m]], orgMembership), -1, drop = FALSE]) * nu.vec) * Ysd

        if (m == 1) {
          beta.m[, 1] <- beta.m[, 1] + Ymean
          beta <<- beta.m
        }
        else {
          beta <<- beta + beta.m
        }

        NULL
      })

    }
    
  }## end ntree=1 end ntree =1 end 

  ##############################################################################
  ##
  ##
  else{#####FOREST BASE LEARNER ####
  ##
  ##
  #############################################################################

    
    ##--------------------------------------------------------------
    ## WE NEED TO REPRODUCE THE WEIGHTED LEAST SQUARES BETA SOLUTION 
    ## WEIGHTS ARE OBTAINED FROM THE FOREST USING THE TEST DATA
    ##
    ## gm        training gradient
    ## Xnew      training pseudo x matrix
    ## pen       training penalty matrix
    ##---------------------------------------------------------

    nullObj <- lapply(1:M, function(m) {

      ##---------------------------------------------------------      
      ## extract gm, Xnew and pen
      ##---------------------------------------------------------
      gm <- baselearner[[m]]$gm
      Xnew <- baselearner[[m]]$Xnew
      pen <- baselearner[[m]]$pen
      
      ##---------------------------------------------------------
      ## determine the test set forest weights
      ##---------------------------------------------------------
      forest.wt <- predict.rfsrc(baselearner[[m]]$forest,
                                 newdata = X,
                                 importance = "none",
                                 forest.wt = TRUE)$forest.wt
                      
      ##---------------------------------------------------------
      ## iterate over cases i to get beta; eliminate cases with small forest weights
      ##---------------------------------------------------------
      beta.m.org <- do.call("cbind", papply(1:n, function(i) {
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
          XnewSum <- XnewSum + pen
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
      ## update beta
      ##---------------------------------------------------------
      beta.m <- t(beta.m.org * nu.vec * Ysd)
      if (m == 1) {
        beta.m[, 1] <- beta.m[, 1] + Ymean
        beta <<- beta.m
      }
      else {
        beta <<- beta + beta.m 
      }

      ##---------------------------------------------------------
      ## update mu
      ##---------------------------------------------------------
      Dbeta.m <- D %*% (beta.m.org * nu.vec)
      if (m == 1) {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      else {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          unlist(mu.list[[m-1]][i]) + Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      NULL##memory saving measure
    })

    ##---------------------------------------------------------
    ## scale mu and add Ymean
    ##---------------------------------------------------------
    mu.list <- lapply(mu.list, function(mlist){  
      lapply(1:n,function(i) {mlist[[i]] * Ysd + Ymean})
    })
    
    ##---------------------------------------------------------
    ## determine a cumulative error rate
    ##---------------------------------------------------------
    if (testFlag) {
      err.rate <- matrix(unlist(lapply(mu.list, function(mlist) {
        c(l1Dist(Y, mlist), l2Dist(Y, mlist))
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
      if (sum(diff.err < Ysd * eps) > 0) {
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
    ## vimp calculation: TBD TBD TBD 
    ##---------------------------------------------------------
    

  }


  ##############################################################################
  ##
  ##  The algorithm is now the same
  ##
  ##############################################################################


  ##---------------------------------------------------------
  ## predicted mu at all unique time points
  ##---------------------------------------------------------
  DbetaT <- D %*% t(beta)
  muhat <- lapply(1:n, function(i) {DbetaT[, i]})
  
  ##---------------------------------------------------------
  ## variable importance
  ##---------------------------------------------------------

  if (vimpFlag) {
    ## vimp when there is no time effect
    if (df.D <= 1) {
      vimp <- unlist(papply(1:p, function(k) {
        DbetaT.k <- D %*% t(beta.vimp[[k]])
        mu.k <- lapply(1:n, function(i) {
          DbetaT.k[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
        l2Dist(Y, mu.k) - err.rate[Mopt, "l2"]
      }))
      names(vimp) <- xvar.names
    }
    ## break vimp into covariate and covariate-time effects
    if (df.D > 1) {
      vimp.cov <- unlist(papply(1:p, function(k) {
        DbetaT.k <- D %*% t(beta.cov.vimp[[k]])
        mu.k <- lapply(1:n, function(i) {
          DbetaT.k[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
        l2Dist(Y, mu.k) - err.rate[Mopt, "l2"]
      }))
      vimp.cov.time <- unlist(papply(1:p, function(k) {
        DbetaT.k <- D %*% t(beta.time.vimp[[k]])
        mu.k <- lapply(1:n, function(i) {
          DbetaT.k[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
        l2Dist(Y, mu.k) - err.rate[Mopt, "l2"]
      }))
      ## modified sapply to lapply here as requested by Amol
      ## previously an error occured for balanced designs with importance = TRUE
      ## the use of sapply was creating a matrix instead of a list
      mu.time <- lapply(1:n, function(i) {
        DbetaT[, i][sample(match(tm[[i]], tm.unq, tm[[i]]))]
      })
      vimp.time <- l2Dist(Y, mu.time) - err.rate[Mopt, "l2"]
      vimp <- c(vimp.cov, vimp.cov.time, vimp.time)
      names(vimp) <- c(xvar.names, paste(xvar.names, "time", sep=":"), "time")
    }
  }

  ##---------------------------------------------------------
  ##return the promised object
  ##---------------------------------------------------------

  ## !!! memory savings !!!!
  ## remove base learner and other objects unnecessary for predict objects
  object$baselearner <- object$membership <- object$gamma <- NULL
  
  pobj <- list(
               boost.obj = object,
               x = X,
               time = tm,
               time.unq = tm.unq,
               y = if (testFlag) Y else NULL,
               mu = mu.list[[Mopt]],
               muhat = muhat,
               phi = object$phi[Mopt],
               rho = object$rho[Mopt],
               err.rate = if (!is.null(err.rate)) err.rate / Ysd else NULL,
               rmse = if (!is.null(err.rate)) err.rate[Mopt, "l2"] / Ysd else NULL,
               vimp = if (!is.null(vimp)) vimp / err.rate[Mopt, "l2"] else NULL,
               Mopt = if (testFlag) Mopt else NULL)

  class(pobj) <- c("boostmtree", "predict", class(object)[3])

  invisible(pobj)


}

