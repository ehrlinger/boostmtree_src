vimpPlot <- function(vimp,
                     Time_Interaction = TRUE,
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
                     seplim = NULL,
                     eps = 0.1,
                     Width_Bar = 1
)
{
  if(is.null(vimp) ){
    stop("vimp is not present in the object")
  }
  if(Time_Interaction){
    p <- floor(length(vimp)/2)
  }else
  {
    p <- length(vimp)  
  }
  
  if(is.null(xvar.names)){
    xvar.names <- paste("x",1:p,sep="")
  }
  
  vimp <- vimp*100
  if(!Time_Interaction){
    ylim <- range(vimp) + c(-0,ymaxlim)
    yaxs <- pretty(ylim)
    yat <- abs(yaxs)
    bp <- barplot(pmax(as.matrix(vimp),0),beside=T,width = Width_Bar,col=col,ylim=ylim,yaxt="n",main = main,cex.lab=cex.lab)
    text(c(bp), pmax(as.matrix(vimp),0) + eps, rep(xvar.names, 3),srt=90,adj= 0.0,cex=if(!is.null(cex.xlab)) cex.xlab else 1)
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
      xaxishead <- c(floor(p/4),floor(p/4))
    }
    bp1 <- barplot(pmax(as.matrix(vimp.x),0),width = Width_Bar,horiz = FALSE,beside=T,col=col,ylim=ylim,yaxt="n",ylab = ylbl,cex.lab=cex.lab,
                   main = main,line = 2)
    text(c(bp1), pmax(as.matrix(vimp.x),0) + eps, rep(xvar.names, 3),srt=90,adj= 0.0,cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[2],yaxishead[2],labels = subhead.labels[2],cex = subhead.cexval)
    bp2 <- barplot(-pmax(as.matrix(vimp.time),0) - eps,width = Width_Bar,horiz = FALSE,beside=T,col=col,add=TRUE,yaxt="n")
    #text(c(bp2), -4, rep(xvar.names, 3),srt=270,adj= 0,yaxt="n",cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[1],-yaxishead[1],labels = subhead.labels[1],cex = subhead.cexval)
    axis(2,yaxs,yat)
  }
}