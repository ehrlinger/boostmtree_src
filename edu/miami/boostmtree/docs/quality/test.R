## #########################################################################
##
##  Outputs:
##    test.txt
##    pred.txt
##    size.txt
##    baseline.txt
##
## #########################################################################

options(object.size=Inf,expressions=100000,memory=Inf,width=150)

output.path = paste(getwd(), "/output/", sep="")

if (!file.exists(path = output.path)) {
  dir.create(path = output.path, showWarnings = TRUE)
}

###########################################################################
#  Load the package and its dependencies.
###########################################################################
use.package = gsub("-use.package=","",grep("-use.package=",commandArgs(),value=T));

if ( length(use.package) == 0) {
  use.package = FALSE
}
if (use.package == "TRUE") {
  use.package = TRUE
} else {
  use.package = FALSE
}

if (use.package) {
  library("_PROJECT_PACKAGE_NAME_")
} else {
  source("load.R")
}

update.benchmark = gsub("-update.benchmark=","",grep("-update.benchmark=",commandArgs(),value=T));

if ( length(update.benchmark) == 0) {
  update.benchmark = FALSE
}
if (update.benchmark == "TRUE") {
  update.benchmark = TRUE
} else {
  update.benchmark = FALSE
}

cat("\n use.package = ", use.package, "\n")
cat("\n update.benchmark = ", update.benchmark, "\n")

###########################################################################

##----------------------------------------------------------------------------
## sim1: iid no time effect low noise quadratic effect
##----------------------------------------------------------------------------
M <- 100
n <- 100
ntest <- 100
N <- 15
phi <- 0.1
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, 3 * N))
  eps <- rnorm(Ni, sd = sqrt(phi))
  x <- rnorm(2)
  y <- 1.5 + 2.5 * x[1] - 3.5 * x[2]^2  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        1:Ni, rep(i, Ni), y)
})))
colnames(dta) <- c("x1", "x2", "time", "id", "y")
trn <- c(1:sum(dta$id <= n))

boost.grow <- boostmtree(dta[trn, 1:2], dta[trn, 3], dta[trn, 4], dta[trn, 5], rho = 0,
                          M = M, K = 75, nu = .1, ntree = 1, lambda = -1, d = 0)
plot(boost.grow)

boost.pred <- predict(boost.grow, dta[-(trn), 1:2], dta[-(trn), 3], dta[-(trn), 4], dta[-(trn), 5])
plot(boost.pred)

## error rates
cat("true model mse:", boostmtree:::gls.mse("y ~ x1 + I(x2^2)", dta, trn), "\n")
cat("linear gls mse:", boostmtree:::gls.mse("y ~ x1 + x2", dta, trn), "\n")
cat("boostmtree mse:", boost.pred$mse, "\n")

##----------------------------------------------------------------------------
## sim2: high correlation no time effect low noise quadratic effect
##----------------------------------------------------------------------------
M <- 100
n <- 100
ntest <- 100
N <- 5
rho <- .80
phi <- .1
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, 3 * N))
  corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
  diag(corr.mat) <- 1
  eps <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
  x <- rnorm(2)
  y <- 1.5 + 2.5 * x[1] - 3.5 * x[2]^2  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        1:Ni, rep(i, Ni), y)
})))
colnames(dta) <- c("x1", "x2", "time", "id", "y")
trn <- c(1:sum(dta$id <= n))

##adaptive rho
boost.grow.1 <- boostmtree(dta[trn, 1:2], dta[trn, 3], dta[trn, 4], dta[trn, 5],
                          M = M, K = 5, nu = .1, ntree = 1, lambda = -1, d = 0)
boost.pred.1 <- predict(boost.grow.1, dta[-(trn), 1:2], dta[-(trn), 3], dta[-(trn), 4], dta[-(trn), 5])

##fixed rho = 0
boost.grow.2 <- boostmtree(dta[trn, 1:2], dta[trn, 3], dta[trn, 4], dta[trn, 5], rho = 0,
                          M = M, K = 5, nu = .1, ntree = 1, lambda = -1, d = 0)
boost.pred.2 <- predict(boost.grow.2, dta[-(trn), 1:2], dta[-(trn), 3], dta[-(trn), 4], dta[-(trn), 5])

##fixed rho = .99
boost.grow.3 <- boostmtree(dta[trn, 1:2], dta[trn, 3], dta[trn, 4], dta[trn, 5], rho = .99,
                          M = M, K = 5, nu = .1, ntree = 1, lambda = -1, d = 0)
boost.pred.3 <- predict(boost.grow.3, dta[-(trn), 1:2], dta[-(trn), 3], dta[-(trn), 4], dta[-(trn), 5])

## error rates
cat("true model mse            :", boostmtree:::gls.mse("y ~ x1 + time * I(x2^2)", dta, trn), "\n")
cat("linear gls mse            :", boostmtree:::gls.mse("y ~ x1*time + x2*time", dta, trn), "\n")
cat("boostmtree mse            :", boost.pred.1$mse, "\n")
cat("boostmtree mse (rho = 0)  :", boost.pred.2$mse, "\n")
cat("boostmtree mse (rho = .99):", boost.pred.3$mse, "\n")

##----------------------------------------------------------------------------
## sim3: high correlation quadratic time with quadratic interaction
##       setting N low gives advantage to correct rho model
##----------------------------------------------------------------------------
M <- 100
n <- 100
ntest <- 100
N <- 5
rho <- .80
phi <- 1
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, 3 * N))
  corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
  diag(corr.mat) <- 1
  eps <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
  x1 <- rnorm(1)
  x2 <- runif(1, 2, 3)
  x <- c(x1, x2)
  tm <- sample((1:(3 * N))/N, size = Ni, replace = TRUE)
  y <- 1.5 + 2.5 * x1 - .65 * (tm ^ 2) * (x2 ^ 2)  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        tm, rep(i, Ni), y)
})))
colnames(dta) <- c("x1", "x2", "time", "id", "y")
trn <- c(1:sum(dta$id <= n))

##adaptive rho
boost.grow.1 <- boostmtree(dta[trn, 1:2], dta[trn, 3], dta[trn, 4], dta[trn, 5],
                          M = M, K = 5, nu = .1, ntree = 1, lambda = -1)
boost.pred.1 <- predict(boost.grow.1, dta[-(trn), 1:2], dta[-(trn), 3], dta[-(trn), 4], dta[-(trn), 5])

##fixed rho = 0
boost.grow.2 <- boostmtree(dta[trn, 1:2], dta[trn, 3], dta[trn, 4], dta[trn, 5], rho = 0,
                          M = M, K = 5, nu = .1, ntree = 1, lambda = -1)
boost.pred.2 <- predict(boost.grow.2, dta[-(trn), 1:2], dta[-(trn), 3], dta[-(trn), 4], dta[-(trn), 5])

##fixed rho = .99
boost.grow.3 <- boostmtree(dta[trn, 1:2], dta[trn, 3], dta[trn, 4], dta[trn, 5], rho = .99,
                          M = M, K = 5, nu = .1, ntree = 1, lambda = -1)
boost.pred.3 <- predict(boost.grow.3, dta[-(trn), 1:2], dta[-(trn), 3], dta[-(trn), 4], dta[-(trn), 5])


## error rates
cat("true model mse            :", boostmtree:::gls.mse("y ~ x1 + I(x2^2) * I(time^2)", dta, trn), "\n")
cat("linear gls mse            :", boostmtree:::gls.mse("y ~ x1 + x2 + x1*time + x2*time", dta, trn), "\n")
cat("boostmtree mse            :", boost.pred.1$mse, "\n")
cat("boostmtree mse (rho = 0)  :", boost.pred.2$mse, "\n")
cat("boostmtree mse (rho = .99):", boost.pred.3$mse, "\n")


## partial plots
partialPlot(boost.grow.1, "x1")
partialPlot(boost.grow.1, "x2")

##----------------------------------------------------------------------------
## sim4 model: same as sim3 but with noisy variables
##----------------------------------------------------------------------------
M <- 100
n <- 100
ntest <- 100
N <- 5
rho <- .80
phi <- 1
p <- 5
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, 3 * N))
  corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
  diag(corr.mat) <- 1
  eps <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
  x1 <- rnorm(1)
  x2 <- runif(1, 2, 3)
  xnoise <- rnorm(p - 2)
  x <- c(x1, x2, xnoise)
  tm <- sample((1:(3 * N))/N, size = Ni, replace = TRUE)
  y <- 1.5 + 2.5 * x1 - .65 * (tm ^ 2) * (x2 ^ 2)  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        tm, rep(i, Ni), y)
})))
colnames(dta) <- c(paste("x", 1:p, sep = ""), "time", "id", "y")
trn <- c(1:sum(dta$id <= n))

##adaptive rho
boost.grow.1 <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3],
                          M = M, K = 5, nu = .1, ntree = 1, lambda = 1000)
boost.pred.1 <- predict(boost.grow.1,
                        dta[-(trn), 1:p], dta[-(trn), p+1], dta[-(trn), p+2], dta[-(trn), p+3])

##fixed rho = 0
boost.grow.2 <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3], rho = 0,
                          M = M, K = 5, nu = .1, ntree = 1, lambda = 1000)
boost.pred.2 <- predict(boost.grow.2,
                        dta[-(trn), 1:p], dta[-(trn), p+1], dta[-(trn), p+2], dta[-(trn), p+3])

##fixed rho = .99
boost.grow.3 <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3], rho = .99,
                          M = M, K = 5, nu = .1, ntree = 1, lambda = 1000)
boost.pred.3 <- predict(boost.grow.3, 
                        dta[-(trn), 1:p], dta[-(trn), p+1], dta[-(trn), p+2], dta[-(trn), p+3])


cat("true model mse            :", boostmtree:::gls.mse("y ~ x1 + I(x2^2) * I(time^2)", dta, trn), "\n")
cat("linear gls mse            :", boostmtree:::gls.mse("y ~ x1 + x2 + x1*time + x2*time", dta, trn), "\n")
cat("boostmtree mse            :", boost.pred.1$mse, "\n")
cat("boostmtree mse (rho = 0)  :", boost.pred.2$mse, "\n")
cat("boostmtree mse (rho = .99):", boost.pred.3$mse, "\n")


##----------------------------------------------------------------------------
## spirometry data
##----------------------------------------------------------------------------

## get the data
## make some useable names
load("./data/spirometry.RData")
testFrac <- c(0, 0.1, 0.2, 0.5)[3]

## training/testing split
if (testFrac > 0 && testFrac < 1) {
  idUnq <- sort(unique(spirometry$id))
  train <- which(is.element(spirometry$id, sample(idUnq, length(idUnq)  * (1 - testFrac), replace = FALSE)))
  test <- setdiff(1:nrow(spirometry$features), train)
} else {
  test <- train <- 1:nrow(spirometry)
}

## calls
boost.grow <- boostmtree(spirometry$features[train,], spirometry$time[train], spirometry$id[train], spirometry$y[train],
                         M = 200, nu = .05, ntree = 1, K = 15, nknots = 5, d = 3, lambda = -1)
boost.pred <- predict(boost.grow, spirometry$features[test,], spirometry$time[test],
                      spirometry$id[test], spirometry$y[test])

plot(boost.grow)
plot(boost.pred)

## performance
fmla <- as.formula(paste("y ~ ", paste(colnames(spirometry$features)[c(1:7,21)], collapse= "+")))
cat("linear gls mse  :", boostmtree:::gls.mse(fmla,
       data.frame(spirometry$features,time=spirometry$time,id=spirometry$id,y=spirometry$y), train), "\n")
cat("boostmtree error:", boost.pred$mse, "\n")

## partial plot of age versus lung treatment
mu.0 <- partialPlot(boost.grow, "AGE", tm = 1:10, plot.it = FALSE, subset=boost.grow$x$DOUBLE==0)$l.obj[[1]]
mu.1 <- partialPlot(boost.grow, "AGE", tm = 1:10, plot.it = FALSE, subset=boost.grow$x$DOUBLE==1)$l.obj[[1]]
mu <- partialPlot(boost.grow, "AGE", tm = 1:10, plot.it = FALSE)$l.obj[[1]]
tm.idx <- 5

pdf("boostmtree_spirometry.pdf",width=10,height=10)
par(mar=c(6,4,2,2))
plot(boost.pred)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(9,10,4,2), mgp = c(6,2,0))
plot(range(mu[, 1]), range(c(mu.0[, -1], mu.1[, -1], mu[, -1])),
     xlab = "age", ylab = "predicted y (adjusted)", type = "n")
lines(mu.0[, 1], mu.0[, -1][, tm.idx], lty = 1, lwd = 6, col = "blue")
lines(mu.1[, 1], mu.1[, -1][, tm.idx], lty = 1, lwd = 6, col = "red")
lines(mu[, 1], mu[, -1][, tm.idx], lty = 1, lwd = 6, col = "black")
legend("topright", legend = c("SLTx", "DLTx", "combined"), lty = 1,
       fill = c(4, 2, 1), cex = 2.5)
dev.off()

