# Gibbs sampler for factorized covariance estimation
# using mgps prior on factor loadings from Battacharya Dunson (2011)
# prior hyperparameters informed by Durante (2017)
# version for adaptive number of factors weighted by iteration

# ARGUMENTS: X: Data matrix (n x p); 
#            nrun: number of iterations;
#            burn: burn-in period;
#            thin: thinning interval;
#            prop: proportion of elements in each column less than epsilon in magnitude cutoff;
#            epsilon: tolerance;
#            kinit: initial value for the number of factors;
#            adapt: logical. Whether or not to adapt number of factors across sampling
#            output: output type, a vector including some of:
#            c("covMean", "covSamples", "factSamples", "sigSamples", "numFactors");
#            verbose: logical. Show progress bar?
#            dump: logical. Save output object during sampling?
#            filename: if dump, filename for output
#            buffer: if dump, frequency of saving
#            augment: additional sampling steps as an expression

linearMGSP = function(X, nrun, burn, thin = 1, prop = 1, epsilon = 1e-3,
                      kinit = NULL, adapt = TRUE, output = c("covMean", "covSamples", 
                                                             "factSamples", "sigSamples", 
                                                             "numFactors"), 
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
  
  as = 1                          # gamma hyperparameters for residual precision
  bs = 0.3                        
  df = 3                          # gamma hyperparameters for t_{ij}
  ad1 = 2.1
  bd1 = 1                         # gamma hyperparameters for delta_1
  ad2 = 3.1
  bd2 = 1                         # gamma hyperparameters delta_h, h >= 2
  adf = 1
  bdf = 1 
  b0 = 1
  b1 = 0.0005
  
  if(is.null(kinit)) kinit = floor(log(p)*3)
  
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  
  VY= apply(X, 2, var)                  # explicitly preserve scale
  scaleMat = sqrt((VY) %*% t(VY))
  X = scale(X)
  num = 0
  k=kinit                               # no. of factors to start with
  
  # --- Initial values --- #
  ps = rgamma(p, as, bs)                           # Sigma = diagonal residual covariance
  lambda = matrix(1, nrow = p, ncol = k)           # latent factor distribution = standard normal
  
  psijh = matrix(rgamma(p*k, df/2, df/2), nrow = p, ncol = k)     # local shrinkage coefficients
  delta = c(rgamma(1,ad1,bd1), rgamma(k-1,ad2,bd2))       # gobal shrinkage coefficients multilpliers
  tauh = cumprod(delta)                                       # global shrinkage coefficients
  Plam = t(t(psijh) * (tauh))                                     # precision of loadings rows
  
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
    lambda = lam_lin(eta, Plam, ps, k, p, X)
    psijh = psi_mg(lambda, tauh, ps, k, p, df)
    delta = del_mg(lambda, psijh, tauh, delta, k, p, ad1, bd1, ad2, bd2)
    tauh = cumprod(delta)
    ps = sig_lin(lambda, eta, k, p, n, X, as, bs)
    Plam = plm_mg(psijh, tauh)
    
    if(!is.null(augment)) eval(augment)
    
    if(adapt){
      # ----- make adaptations ----#
      prob = 1/exp(b0 + b1*i)                    # probability of adapting
      uu = runif(1)
      lind = colSums(abs(lambda) < epsilon)/p    # proportion of elements in each column less than eps in magnitude
      vec = lind >= prop
      num = sum(vec)                             # number of redundant columns
      
      if(uu < prob) {
        if((i > 20) & (num == 0) & all(lind < 0.995)) {
          k = k + 1
          lambda = cbind(lambda, rep(0,p))
          eta = cbind(eta,rnorm(n))
          psijh = cbind(psijh, rgamma(p,df/2,df/2))
          delta[k] = rgamma(1, ad2,bd2)
          tauh = cumprod(delta)
          Plam = t(t(psijh) * tauh)
        } else {
          if (num > 0) {
            k = max(k - num,1)
            lambda = lambda[,!vec, drop = F]
            psijh = psijh[,!vec, drop = F]
            eta = eta[,!vec, drop = F]
            delta = delta[!vec]
            tauh = cumprod(delta)
            Plam = t(t(psijh) * tauh)
          }
        }
      }
    }
    
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
