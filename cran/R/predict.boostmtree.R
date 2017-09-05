predict.boostmtree <-
  function(object,
           x,
           tm,
           id,
           y,
           M,
           importance = TRUE,
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
                                               importance = importance,
                                               eps = eps,
                                               ...)

  return(result.predict)

}
