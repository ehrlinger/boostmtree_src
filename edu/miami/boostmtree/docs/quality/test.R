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

library(mvtnorm)

##----------------------------------------------------------------------------
## internal functions
##----------------------------------------------------------------------------

## dgm: test-set error
dgm <- function(f, dta, trn) {
  f <- as.formula(f)
  gls.grow <- gls(f, data = dta[trn, ], corr = corCompSymm(form = ~ 1 | id))
  gls.pred  <- tapply(model.matrix(f, dta[-trn,]) %*% gls.grow$coef,
                 dta[-trn, "id"], function(x) {x})
  y.test <- tapply(dta[-trn, "y"], dta[-trn, "id"], function(x) {x})
  l2Dist(gls.pred, y.test)
}


##----------------------------------------------------------------------------
## sim1: linear model: iid no time effect
##----------------------------------------------------------------------------
n <- 100
ntest <- 100
N <- 15
phi <- .25
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, N))
  eps <- rnorm(Ni, sd = sqrt(phi))
  x <- rnorm(2)
  y <- 1.5 + 2.5 * x[1] - 3.5 * x[2]^2  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        1:Ni, rep(i, Ni), y)
})))
colnames(dta) <- c("x1", "x2", "time", "id", "y")
trn <- c(1:sum(dta$id <= 100))

boost.grow <- boostmtree(dta[trn, 1:2], dta[trn, 3], dta[trn, 4], dta[trn, 5], rho = 0,
                          M = 100, K = 75, nu = .1, ntree = 1, lambda = -1)

pdf(paste(output.path, "sim1.grow.pdf", sep=""), width=10, height=8)
plot(boost.grow)
dev.off()

boost.pred <- predict(boost.grow, dta[-(trn), 1:2], dta[-(trn), 3], dta[-(trn), 4], dta[-(trn), 5])
pdf(paste(output.path, "sim1.pred.pdf", sep=""), width=10, height=8)
plot(boost.pred)
dev.off()

cat("DGM error       :", dgm("y ~ x1 + I(x2^2)", dta, trn), "\n")
cat("linear gls error:", dgm("y ~ x1 + x2", dta, trn), "\n")
cat("boostmtree error:", boost.pred$err.rate[boost.pred$Mopt, ], "\n")


##----------------------------------------------------------------------------
## sim2: linear model: iid with time effect
##----------------------------------------------------------------------------
n <- 100
ntest <- 100
N <- 25
phi <- .25
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, N))
  eps <- rnorm(Ni, sd = sqrt(phi))
  x <- rnorm(2)
  tm <- 1:Ni
  y <- 1.5 + 2.5 * x[1] - 0.5 * tm * x[2]^2  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        tm, rep(i, Ni), y)
})))
colnames(dta) <- c("x1", "x2", "time", "id", "y")
trn <- c(1:sum(dta$id <= 100))

boost.grow <- boostmtree(dta[trn, 1:2], dta[trn, 3], dta[trn, 4], dta[trn, 5],
                          M = 100, K = 75, nu = .1, ntree = 1, lambda = -1)

pdf(paste(output.path, "sim2.grow.pdf", sep=""), width=10, height=8)
plot(boost.grow)
dev.off()

boost.pred <- predict(boost.grow, dta[-(trn), 1:2], dta[-(trn), 3], dta[-(trn), 4], dta[-(trn), 5])
pdf(paste(output.path, "sim2.pred.pdf", sep=""), width=10, height=8)
plot(boost.pred)
dev.off()

cat("DGM error       :", dgm("y ~ x1 + time * I(x2^2)", dta, trn), "\n")
cat("linear gls error:", dgm("y ~ x1*time + x2*time", dta, trn), "\n")
cat("boostmtree error:", boost.pred$err.rate[boost.pred$Mopt, ], "\n")

##----------------------------------------------------------------------------
## sim3: linear model: iid with time effect + noise
##----------------------------------------------------------------------------
n <- 100
ntest <- 100
N <- 25
phi <- .25
p <- 5
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, N))
  eps <- rnorm(Ni, sd = sqrt(phi))
  x <- rnorm(p)
  tm <- 1:Ni
  y <- 1.5 + 2.5 * x[1] - 0.5 * tm * x[2]^2  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        tm, rep(i, Ni), y)
})))
colnames(dta) <- c(paste("x", 1:p, sep = ""), "time", "id", "y")
trn <- c(1:sum(dta$id <= 100))

boost.grow <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3],
                          M = 100, K = 50, nu = .1, ntree = 1, lambda = -1)

pdf(paste(output.path, "sim3.grow.pdf", sep=""), width=10, height=8)
plot(boost.grow)
dev.off()

boost.pred <- predict(boost.grow,
                      dta[-(trn), 1:p], dta[-(trn), p+1], dta[-(trn), p+2], dta[-(trn), p+3])

pdf(paste(output.path, "sim3.pred.pdf", sep=""), width=10, height=8)
plot(boost.pred)
dev.off()

cat("DGM error       :", dgm("y ~ x1 + time * I(x2^2)", dta, trn), "\n")
cat("linear gls error:", dgm("y ~ x1*time + x2*time", dta, trn), "\n")
cat("boostmtree error:", boost.pred$err.rate[boost.pred$Mopt, ], "\n")

pdf(paste(output.path, "sim3.partial.pdf", sep=""), width=10, height=8)
partialPlot(boost.grow, "x2")
dev.off()

##----------------------------------------------------------------------------
## sim4 model: correlation with time effect + noisy variables
##----------------------------------------------------------------------------
n <- 100
ntest <- 100
N <- 25
phi <- .25
rho <- .35
p <- 5
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, N))
  corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
  diag(corr.mat) <- 1
  eps <- c(rmvnorm(n = 1, mean = rep(0, Ni), sigma = phi * corr.mat))
  x <- rnorm(p)
  tm <- 1:Ni
  y <- 1.5 + 2.5 * x[1] - 0.5 * tm * x[2]^2  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        tm, rep(i, Ni), y)
})))
colnames(dta) <- c(paste("x", 1:p, sep = ""), "time", "id", "y")
trn <- c(1:sum(dta$id <= 100))

boost.grow <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3],
                          M = 500, K = 5, nu = .1, ntree = 1, lambda = -1)

pdf(paste(output.path, "sim4.grow.pdf", sep=""), width=10, height=8)
plot(boost.grow)
dev.off()

boost.pred <- predict(boost.grow,
             dta[-(trn), 1:p], dta[-(trn), p+1], dta[-(trn), p+2], dta[-(trn), p+3])

pdf(paste(output.path, "sim4.pred.pdf", sep=""), width=10, height=8)
plot(boost.pred)
dev.off()

pdf(paste(output.path, "sim4.partial.pdf", sep=""), width=10, height=8)
boost.partial <- partialPlot(boost.grow, "x2", plot.it = FALSE)
dev.off()

cat("DGM error       :", dgm("y ~ x1 + time * I(x2^2)", dta, trn), "\n")
cat("linear gls error:", dgm("y ~ x1*time + x2*time", dta, trn), "\n")
cat("boostmtree error:", boost.pred$err.rate[boost.pred$Mopt, ], "\n")

pdf(paste(output.path, "sim4.pred.nonlinear_correlation.pdf", sep =""), width=10,height=8)
par(mar=c(6,4,2,2))
plot(boost.pred)
par(cex.axis=3.5,cex.lab=3.5,cex.main=3.5,mar=c(9,10,4,2), mgp = c(6,2,0))
 matplot(boost.partial$l.obj[[1]][, 1],boost.partial$l.obj[[1]][, -1],
         type="l",lty=1,col=1,lwd=3,xlab=expression(x[2]),ylab="predicted y (adjusted)")
dev.off()


##----------------------------------------------------------------------------
## sim5 model: correlation with time effect + noisy variables
## compare results under different rho settings
##----------------------------------------------------------------------------
n <- 1000
ntest <- 5000
N <- 25
phi <- .25
rho <- .75
p <- 5
M <- 300
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, N))
  corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
  diag(corr.mat) <- 1
  eps <- c(rmvnorm(n = 1, mean = rep(0, Ni), sigma = phi * corr.mat))
  x <- rnorm(p)
  tm <- 1:Ni
  y <- 1.5 + 2.5 * x[1] - 0.5 * tm * x[2]^2  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        tm, rep(i, Ni), y)
})))
colnames(dta) <- c(paste("x", 1:p, sep = ""), "time", "id", "y")
trn <- c(1:sum(dta$id <= 100))

boost.grow.1 <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3],
                          M = M, K = 5, nu = .05, ntree = 1, lambda = -1)
boost.grow.2 <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3],
                          rho = 0, M = M, K = 5, nu = .05, ntree = 1, lambda = -1)
boost.grow.3 <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3],
                          rho = 0.75, M = M, K = 5, nu = .05, ntree = 1, lambda = -1)

boost.pred.1 <- predict(boost.grow.1,
             dta[-(trn), 1:p], dta[-(trn), p+1], dta[-(trn), p+2], dta[-(trn), p+3])
boost.pred.2 <- predict(boost.grow.2,
             dta[-(trn), 1:p], dta[-(trn), p+1], dta[-(trn), p+2], dta[-(trn), p+3])
boost.pred.3 <- predict(boost.grow.3,
             dta[-(trn), 1:p], dta[-(trn), p+1], dta[-(trn), p+2], dta[-(trn), p+3])

cat("DGM error               :", dgm("y ~ x1 + time * I(x2^2)", dta, trn), "\n")
cat("linear gls error        :", dgm("y ~ x1*time + x2*time", dta, trn), "\n")
cat("boostmtree error        :", boost.pred.1$err.rate[boost.pred.1$Mopt, "l2"], "\n")
cat("boostmtree error rho=0  :", boost.pred.2$err.rate[boost.pred.2$Mopt, "l2"], "\n")
cat("boostmtree error rho=.75:", boost.pred.3$err.rate[boost.pred.3$Mopt, "l2"], "\n")


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
                         M = 500, nu = .015, ntree = 1, K = 15, nknots = 5, d = 3, lambda = -1)
boost.pred <- predict(boost.grow, spirometry$features[test,], spirometry$time[test],
                      spirometry$id[test], spirometry$y[test])

pdf(paste(output.path, "spirometry.grow.pdf", sep=""), width=10,height=8)
plot(boost.grow)
dev.off()

pdf(paste(output.path, "spirometry.pred.pdf", sep=""), width=10,height=8)
plot(boost.pred)
dev.off()

## performance
fmla <- as.formula(paste("y ~ ", paste(colnames(spirometry$features)[c(1:7,21)], collapse= "+")))
cat("linear gls error:", dgm(fmla,
       data.frame(spirometry$features,time=spirometry$time,id=spirometry$id,y=spirometry$y), train), "\n")
cat("boostmtree error:", boost.pred$err.rate[boost.pred$Mopt, ], "\n")

## partial plot of age versus lung treatment
mu.0 <- partialPlot(boost.grow, "AGE", tm = 1:10, plot.it = FALSE, subset=boost.grow$x$DOUBLE==0)$l.obj[[1]]
mu.1 <- partialPlot(boost.grow, "AGE", tm = 1:10, plot.it = FALSE, subset=boost.grow$x$DOUBLE==1)$l.obj[[1]]
mu <- partialPlot(boost.grow, "AGE", tm = 1:10, plot.it = FALSE)$l.obj[[1]]
tm.idx <- 5


pdf(paste(output.path, "spirometry.pred.2.pdf", sep=""), width=10,height=8)
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

