# Gibbs sampler for factorized covariance estimation
# using dirichlet-laplace prior on factor loadings from Battacharya Dunson (2014)
# same eta, lambda, and sigma updates as MGSP linear

# ARGUMENTS: X: Data matrix (n x p); 
#            nrun: number of iterations;
#            burn: burn-in period;
#            thin: thinning interval;
#            prop: proportion of elements in each column less than epsilon in magnitude cutoff;
#            epsilon: tolerance;
#            k: initial value for the number of factors;
#            output: output type, a vector including some of:
#            c("covMean", "covSamples", "factSamples", "sigSamples");
#            verbose: logical. Show progress barr?
#            dump: Save output object during sampling?
#            filename: if dump, filename for output
#            buffer: if dump, frequency of saving
#            augment: additional sampling steps as an expression

linearDL = function(X, nrun, burn, thin = 1, prop = 1, epsilon = 1e-3,
                    k = NULL, output = c("covMean", "covSamples", 
                                             "factSamples", "sigSamples"), 
                    verbose = TRUE, dump = FALSE, filename = "samps.Rds", 
                    buffer = 10000, augment = NULL){
  
  if(nrun <= burn) stop("nrun must be larger than burn")
  if(!is.matrix(X)) stop("X must be a matrix")
  if(any(is.na(X))) stop("X cannot contain missing data")
  if(!is.null(augment)) if(!is.expression(augment)) stop("augment must be an expression (see expression())")
  
  cm = any(output %in% "covMean")
  cs = any(output %in% "covSamples")
  fs = any(output %in% "factSamples")
  ss = any(output %in% "sigSamples")
  nf = any(output %in% "numFactors")
  
  p = ncol(X)
  n = nrow(X)

  a = 1/2
  as = 1
  bs = 0.3
  
  if(is.null(k)) k = floor(log(p)*3)
  
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  
  VY= apply(X, 2, var)                  # explicitly preserve scale
  scaleMat = sqrt((VY) %*% t(VY))
  X = scale(X)
  
  # --- Initial values --- #
  ps = rgamma(p,as,bs)
  lambda = matrix(0,p,k)
  eta = matrix(rnorm(n*k),n,k)
  
  tau = rgamma(p,p*a,1/2)
  psi = matrix(rexp(p*k,1/2),p,k)
  phi = matrix(0,p,k)
  for(j in 1:p){
    gam = rgamma(k, a)
    phi[j,] = gam /sum(gam)
  }
  Plam = phi*(phi^2)*matrix(rep(tau^2,k),p,k,byrow=F)
  
  if(cm) COVMEAN = matrix(0, nrow = p, ncol = p)
  if(cs) OMEGA = array(dim = c(p, p, sp))
  if(fs) {LAMBDA = list()
          ETA = list()}
  if(ss) SIGMA = array(dim = c(p, sp))
  if(nf) K = rep(NA, sp)
  ind = 1
  
  at = ceiling(nrun/100)
  if(verbose) {
    pb = txtProgressBar(style = 3)
  }
  
  for(i in 1:nrun){
    eta = eta_lin(lambda, ps, k, n, X)
    Plam = plm_dl(psi, phi, tau)
    lambda = lam_lin(eta, Plam, ps, k, p, X)
    psi = psi_dl(lambda, phi, tau)
    tau = tau_dl(lambda, phi, k, p)
    phi = phi_dl(lambda, a, k, p)
    ps = sig_lin(lambda, eta, k, p, n, X, as, bs)
    
    if(!is.null(augment)) eval(augment)
    
    if((i %% thin == 0) & (i > burn)) {
      if(cm | cs) Omega = (tcrossprod(lambda) + diag(1/c(ps))) * scaleMat
      if(cm) COVMEAN = COVMEAN + Omega / sp
      if(cs) OMEGA[,,ind] = Omega
      if(fs) {LAMBDA[[ind]] = lambda
              ETA[[ind]] = eta}
      if(ss) SIGMA[,ind] = 1/ps
      if(nf) K[ind] = k
      ind = ind + 1
    }
    
    if(verbose & (i %% at == 0)) setTxtProgressBar(pb, i / nrun)
    
    if(dump & (i %% buffer == 0)){
      out = list()
      if(cm) out = c(out, list(covMean = COVMEAN))
      if(cs) out = c(out, list(omegaSamps = OMEGA))
      if(fs) out = c(out, list(lambdaSamps = LAMBDA, etaSamps = ETA))
      if(ss) out = c(out, list(sigmaSamps = SIGMA))
      if(nf) out = c(out, list(numFacts = K))
      saveRDS(out, filename, compress = FALSE)
    }
  }
  
  if(verbose) close(pb)
  
  out = list()
  if(cm) out = c(out, list(covMean = COVMEAN))
  if(cs) out = c(out, list(omegaSamps = OMEGA))
  if(fs) out = c(out, list(lambdaSamps = LAMBDA, etaSamps = ETA))
  if(ss) out = c(out, list(sigmaSamps = SIGMA))
  if(nf) out = c(out, list(numFacts = K))
  
  return(out)
}
