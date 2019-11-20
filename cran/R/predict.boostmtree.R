predict.boostmtree <- function(object,
                               x,
                               tm,
                               id,
                               y,
                               M,
                               eps = 1e-5,
                               ...)
{

  ## call generic predict 
  result.predict <- generic.predict.boostmtree(object = object,
                                               x,
                                               tm,
                                               id,
                                               y,
                                               M,
                                               eps = eps,
                                               ...)

  return(result.predict)
}
