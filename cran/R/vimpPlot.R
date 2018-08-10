vimpPlot <- function(object,
                     xvar.names = NULL,
                     cex.xlab = NULL,
                     ymaxlim = 0,
                     ymaxtimelim = 0,
                     subhead.cexval = 1,
                     yaxishead = NULL,
                     xaxishead = NULL,
                     main = "Variable Importance (%)",
                     col = grey(.80),
                     cex.lab = 1.5,
                     subhead.labels = c("Time-Interactions Effects","Main Effects"),
                     ylbl = FALSE,
                     seplim = NULL
)
{

  if(is.null(object$vimp) ){
    stop("vimp is not present in the object")
  }

  vimp <- object$vimp
  if(is.null(xvar.names)){
    xvar.names <- colnames(object$x)
  }
  p <- ncol(object$x)
  n.vimp <- length(vimp)
  if(n.vimp == p ){
    univariate <- TRUE
  }else
  {
    univariate <- FALSE
  }

  if(univariate){
    vimp <- vimp*100
  }else
  {
    vimp <- (vimp[-n.vimp])*100
    n.vimp <- length(vimp)
  }

  if(univariate){

    ylim <- range(vimp) + c(0,ymaxlim)
    yaxs <- pretty(ylim)
    yat <- abs(yaxs)
    bp <- barplot(as.matrix(vimp),beside=T,col=col,ylim=ylim,yaxt="n",main = main,cex.lab=cex.lab)
    text(c(bp), pmax(as.matrix(vimp),0), rep(xvar.names, 3),srt=90,adj=-0.5,cex= if(!is.null(cex.xlab)) cex.xlab else 1 )
    axis(2,yaxs,yat)
  }else
  {
    vimp.x <- vimp[1:p]
    vimp.time <- vimp[-c(1:p)]
    ylim <- max(c(vimp.x,vimp.time)) * c(-1, 1) + c(-ymaxtimelim,ymaxlim)
    if(ylbl){
      ylbl <- paste("Time-Interactions", "Main Effects", sep = if(!is.null(seplim)) seplim else "                   " )
    }else
    {
      ylbl <- NULL
    }

    yaxs <- pretty(ylim)
    yat <- abs(yaxs)

    if(is.null(yaxishead)){
      yaxishead <- c(-ylim[1],ylim[2])
    }

    if(is.null(xaxishead)){
      xaxishead <- c(floor(n.vimp/4),floor(n.vimp/4))
    }

    bp1 <- barplot(pmax(as.matrix(vimp.x),0),beside=T,col=col,ylim=ylim,yaxt="n",ylab = ylbl,cex.lab=cex.lab,
                   main = main)
    text(c(bp1), pmax(as.matrix(vimp.x),0), rep(xvar.names, 3),srt=90,adj=-0.5,cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[2],yaxishead[2],labels = subhead.labels[2],cex = subhead.cexval)
    bp2 <- barplot(-pmax(as.matrix(vimp.time),0),beside=T,col=col,add=TRUE,yaxt="n")
    text(c(bp2), -pmax(as.matrix(vimp.time),0), rep(xvar.names, 3),srt=90,adj=1.5,yaxt="n",cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[1],-yaxishead[1],labels = subhead.labels[1],cex = subhead.cexval)
    axis(2,yaxs,yat)
  }

}
