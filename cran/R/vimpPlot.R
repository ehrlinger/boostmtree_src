vimpPlot <- function(vimp,
                     Q_set = NULL,
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
                     Width_Bar = 1,
                     path_saveplot = NULL,
                     Verbose = TRUE
                     )
{
  if(is.null(vimp) ){
    stop("vimp is not present in the object")
  }
  p <- nrow(vimp[[1]])
  if(is.null(xvar.names)){
    xvar.names <- paste("x",1:p,sep="")
  }
  vimp <- lapply(vimp,function(v){
    v*100
  })
  n.Q <- ncol(vimp[[1]])
  if(is.null(Q_set)){
    Q_set <- paste("V",seq(n.Q),sep="")
  }
  for(q in 1:n.Q){
   if(is.null(path_saveplot)){
      path_saveplot <- tempdir()
   }
    Plot_Name <- if(n.Q == 1) "VIMPplot.pdf" else paste("VIMPplot_Prob(y = ",Q_set[q],")",".pdf",sep="")
    pdf(file = paste(path_saveplot,"/",Plot_Name,sep=""),width = 10,height = 10)
    
  if(!Time_Interaction){
    ylim <- range(vimp[[1]][,q]) + c(-0,ymaxlim)
    yaxs <- pretty(ylim)
    yat <- abs(yaxs)
    bp <- barplot(pmax(as.matrix(vimp[[1]][,q]),0),beside=T,width = Width_Bar,col=col,ylim=ylim,yaxt="n",main = main,cex.lab=cex.lab)
    text(c(bp), pmax(as.matrix(vimp[[1]][,q]),0) + eps, rep(xvar.names, 3),srt=90,adj= 0.0,cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    axis(2,yaxs,yat)
  }else
  {
    vimp.x <- vimp[[1]][,q]
    vimp.time <- vimp[[2]][,q]
    ylim <- max(c(vimp.x,vimp.time)) * c(-1, 1) + c(-ymaxtimelim,ymaxlim)
    if(ylbl){
      ylabel <- paste("Time-Interactions", "Main Effects", sep = if(!is.null(seplim)) seplim else "                   " )
    }else
    {
      ylabel <- ""
    }
    yaxs <- pretty(ylim)
    yat <- abs(yaxs)
    if(is.null(yaxishead)){
      yaxishead <- c(-ylim[1],ylim[2])
    }
    if(is.null(xaxishead)){
      xaxishead <- c(floor(p/4),floor(p/4))
    }
    bp1 <- barplot(pmax(as.matrix(vimp.x),0),width = Width_Bar,horiz = FALSE,beside=T,col=col,ylim=ylim,yaxt="n",ylab = ylabel,cex.lab=cex.lab,
                   main = main)
    text(c(bp1), pmax(as.matrix(vimp.x),0) + eps, rep(xvar.names, 3),srt=90,adj= 0.0,cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[2],yaxishead[2],labels = subhead.labels[2],cex = subhead.cexval)
    bp2 <- barplot(-pmax(as.matrix(vimp.time),0) - eps,width = Width_Bar,horiz = FALSE,beside=T,col=col,add=TRUE,yaxt="n")
    #text(c(bp2), -4, rep(xvar.names, 3),srt=270,adj= 0,yaxt="n",cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[1],-yaxishead[1],labels = subhead.labels[1],cex = subhead.cexval)
    axis(2,yaxs,yat)
  }
  dev.off()
   if(Verbose){
       cat("Plot will be saved at:",path_saveplot,sep = "","\n")
    }  
 }
}
