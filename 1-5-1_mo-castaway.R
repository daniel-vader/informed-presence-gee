# Modify IrregLong mo function to throw out failed outputations.

mo.castaway <- function (noutput, fn, data, weights, singleobs, id, time, keep.first, 
                 var = TRUE, ...) 
{
  dimlast <- as.numeric(var) + 1
  for (it in 1:noutput) {
    if (it == 1) {
      a <- IrregLong:::outputanalfn(fn, data, weights, singleobs, 
                        id, time, keep.first, ...)
      if (is.vector(a)) 
        a <- array(a, dim = c(1, length(a)))
      if (var) 
        ans <- array(dim = c(noutput, dim(a)))
      else ans <- array(dim = c(noutput, length(a), 1))
      ans[1, , ] <- a
    }
    if (it > 1) {
      a <- IrregLong:::outputanalfn(fn, data, weights, singleobs, 
                        id, time, keep.first, ...)
      if (is.vector(a)) 
        a <- array(a, dim = c(1, length(a)))
      ans[it, , ] <- a
    }
  }
  print(ans)
  pooled <- apply(ans, 2:3, mean, na.rm=T)
  if (var) {
    var.within <- pooled[, 2]
    var.between <- apply(array(ans[, , 1], dim = c(nrow(ans), 
                                                   dim(a)[1])), 2, var, na.rm=T)
    var.est <- var.within - ((nrow(ans) - 1)/nrow(ans)) * var.between
    RE.MO <- 1 + (1/nrow(ans)) * (var.between/(var.within - 
                                               var.between))
    return(list(est = pooled[, 1], se = sqrt(var.est), RE.MO = RE.MO))
  }
  else return(list(est = pooled[, 1]))
}
