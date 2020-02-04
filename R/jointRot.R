# Rotate Lambda and eta jointly
# ARGUMENTS: lambda: list of loadings samples
#            eta: list of factor samples

jointRot = function(lambda, eta){
  vari = lapply(lambda, varimax)
  loads = lapply(vari, `[[`, 1)
  rots = lapply(vari, `[[`, 2)
  rotfact = mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)
  
  norms = sapply(loads, norm, "2")
  piv = loads[order(norms)][[round(length(lambda)/2)]]
  
  matches = lapply(loads, msfOUT, piv)
  
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  etaout = mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)
  
  return(list(lambda = lamout, eta = etaout))
}