predict.boostmtree <- function(object,
                               x,
                               tm,
                               id,
                               y,
                               M,
                               eps = 1e-5,
                               useCVflag = FALSE,
                               ...)
{
  result.predict <- generic.predict.boostmtree(object = object,
                                               x,
                                               tm,
                                               id,
                                               y,
                                               M,
                                               eps = eps,
                                               useCVflag = useCVflag,
                                               ...)
  return(result.predict)
}
