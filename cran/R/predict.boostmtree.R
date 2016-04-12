predict.boostmtree <-
  function(object,
           x,
           tm,
           id,
           y,
           M,
           importance = TRUE,
           verbose = TRUE,
           eps = 1e-5,
           ...)
{


  result.predict <- generic.predict.boostmtree(object = object,
                                               x = x,
                                               tm = tm,
                                               id = id,
                                               y = y,
                                               M = M,
                                               importance = importance,
                                               verbose = verbose,
                                               eps = eps,
                                               ...)

  return(result.predict)

}
