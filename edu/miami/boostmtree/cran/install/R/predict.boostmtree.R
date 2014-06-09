predict.boostmtree <-
  function(obj,
           x,
           tm,
           id,
           y,
           M,
           importance = TRUE,
           proximity = FALSE,
           verbose = TRUE,
           eps = 1e-4,
           ...)
{


  result.predict <- generic.predict.boostmtree(obj = obj,
                                               x = x,
                                               tm = tm,
                                               id = id,
                                               y = y,
                                               M = M,
                                               importance = importance,
                                               proximity = proximity,
                                               verbose = verbose,
                                               eps = eps,
                                               ...)
  
  return(result.predict)

}
