# Gibbs sampler for factorized regression with all 2 way interactions
# using MGSP prior on factor loadings from Battacharya Dunson (2014)
# same eta, lambda, and sigma updates as interactionDL

# ARGUMENTS: Y: Data matrix (n x p); 
#            nrun: number of iterations;
#            burn: burn-in period;
#            thin: thinning interval;
#            prop: proportion of elements in each column less than epsilon in magnitude cutoff;
#            epsilon: tolerance;
#            k: number of factors;
#            output: output type, a vector including some of:
#            c("covMean", "covSamples", "factSamples", "sigSamples", "coefSamples", "errSamples");
#            dump: logical. Save samples while sampling?;
#            filename: if dump optional filename for samples;
#            buffer: if dump how often to save samples;
#            adapt: logical or "burn". Adapt proposal variance in metropolis hastings step? if "burn", 
#                   will adapt during burn in and not after
#            augment: additional sampling steps as an expression

interactionMGSP = function(y, X, nrun, burn, thin = 1, delta_rw = 0.0526749,
                         a = 1/2, k = NULL, output = c("covMean", "covSamples", "factSamples",
                                                           "sigSamples", "coefSamples",
                                                           "errSamples"), verbose = TRUE,
                         dump = FALSE, filename = "samps.Rds", buffer = 10000, adapt = "burn",
                         augment = NULL){
  
  if(nrun <= burn) stop("nrun must be larger than burn")
  if(!is.matrix(X)) stop("X must be a matrix")
  if(any(is.na(X))) stop("X cannot contain missing data")
  if(!is.null(augment)) if(!is.expression(augment)) stop("augment must be an expression (see expression())")
  
  cm = any(output %in% "covMean")
  cs = any(output %in% "covSamples")
  fs = any(output %in% "factSamples")
  ss = any(output %in% "sigSamples")
  cf = any(output %in% "coefSamples")
  es = any(output %in% "errSamples")
  
  p = ncol(X)
  n = nrow(y)
  
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
  
  if(is.null(k)) k = floor(log(p)*3)
  
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  
  VX= apply(X, 2, var)                  # explicitly preserve scale
  scaleMat = sqrt((VX) %*% t(VX))
  X = scale(X)
  
  ssy = 1                     # initial values
  phi = numeric(k)
  ps = rgamma(p,as,bs)
  lambda = matrix(0,p,k)
  eta = matrix(rnorm(n*k),n,k)
  Psi = matrix(0,k,k)
  
  psiMG = matrix(rgamma(p*k, df/2, df/2), nrow = p, ncol = k)     # local shrinkage coefficients
  delta = c(rgamma(1,ad1,bd1), rgamma(k-1,ad2,bd2))       # gobal shrinkage coefficients multilpliers
  tau = cumprod(delta)                                       # global shrinkage coefficients
  Plam = t(t(psiMG) * (tau))                                     # precision of loadings rows

  
  if(cm) COVMEAN = matrix(0, nrow = p, ncol = p)
  if(cs) OMEGA = array(dim = c(p, p, sp))
  if(fs) {LAMBDA = list()
          ETA = list()}
  if(ss) SIGMA = array(dim = c(p, sp))
  if(cf) {
    PHI = array(dim = c(k, sp))
    PSI = array(dim = c(k, k, sp))
    INTERCEPT = numeric(sp)
    MAIN = array(dim = c(p, sp))
    INTERACTION = array(dim = c(p, p, sp))
  }
  if(es) SSY = numeric(sp)
  ind = 1
  acp = numeric(n)
  
  at = ceiling(nrun/100)
  if(verbose) {
    pb = txtProgressBar(style = 3)
  }
  
  for(i in 1:nrun){
    eta = eta_int(lambda, eta, ps, phi, Psi, k, n, y, X, ssy, delta_rw, acp)
    Psi = psi_int(eta, y, phi, ssy, k, n)
    phi = phi_int(eta, y, ssy, Psi, k)
    ssy = ssy_int(eta, phi, Psi, y, n)
    psiMG = psi_mg(lambda, tau, ps, k, p, df)
    delta = del_mg(lambda, psiMG, tau, delta, k, p, ad1, bd1, ad2, bd2)
    tau = cumprod(delta)
    Plam = plm_mg(psiMG, tau)
    lambda = lam_lin(eta, Plam, ps, k, p, X)
    ps = c(sig_lin(lambda, eta, k, p, n, X, as, bs))
    
    if(!is.null(augment)) eval(augment)
    
    if((i %% thin == 0) & (i > burn)) {
      if(cm | cs) Omega = (tcrossprod(lambda) + diag(1/ps)) * scaleMat
      if(cm) COVMEAN = COVMEAN + Omega / sp
      if(cs) OMEGA[,,ind] = Omega
      if(fs) {LAMBDA[[ind]] = lambda
              ETA[[ind]] = eta}
      if(ss) SIGMA[,ind] = 1/ps
      if(cf) {
        SigmaI = diag(ps)
        Lambda.T = t(lambda)
        V_n = solve(Lambda.T%*%SigmaI%*%lambda+diag(k))
        a_n = V_n%*%Lambda.T%*%SigmaI
        dsVX = diag(sqrt(VX))
        dsVX_inv = diag(1/sqrt(VX))
        INTERACTION[,,ind] = dsVX_inv%*%t(a_n)%*%Psi%*%a_n%*%dsVX_inv
        MAIN[,ind] = as.vector(t(phi)%*%a_n)
        INTERCEPT[ind] = sum(diag(Psi%*%V_n))
        PHI[,ind] = phi
        PSI[,,ind] = Psi
      }
      if(es) SSY[ind] = ssy
      ind = ind + 1
    }
    
    if(verbose & (i %% at == 0)) setTxtProgressBar(pb, i / nrun)
    
    if (i%%100==0) {if(adapt == "burn") {if(i <= burn) {
      acp_mean = mean(acp)/100
      if(acp_mean > 0.3){
        delta_rw = delta_rw*2
      }else if(acp_mean < 0.2){
        delta_rw = delta_rw*2/3}
      acp = numeric(n)
    }} else if(adapt) {
      acp_mean = mean(acp)/100
      if(acp_mean > 0.3){
        delta_rw = delta_rw*2
      }else if(acp_mean < 0.2){
        delta_rw = delta_rw*2/3
      }
      acp = numeric(n)
    }}
    
    if(dump & (i %% buffer == 0)){
      out = list()
      if(cm) out = c(out, list(covMean = COVMEAN))
      if(cs) out = c(out, list(omegaSamps = OMEGA))
      if(fs) out = c(out, list(lambdaSamps = LAMBDA,
                               etaSamps = ETA))
      if(ss) out = c(out, list(sigmaSamps = SIGMA))
      if(cf) out = c(out, list(phiSamps = PHI, PsiSamps = PSI, 
                               interceptSamps = INTERCEPT,
                               mainEffectSamps = MAIN,
                               interactionSamps = INTERACTION))
      if(es) out = c(out, list(ssySamps = SSY))
      saveRDS(out, filename, compress = FALSE)
    }
  }
  
  if(verbose) close(pb)
  
  out = list()
  if(cm) out = c(out, list(covMean = COVMEAN))
  if(cs) out = c(out, list(omegaSamps = OMEGA))
  if(fs) out = c(out, list(lambdaSamps = LAMBDA,
                           etaSamps = ETA))
  if(ss) out = c(out, list(sigmaSamps = SIGMA))
  if(cf) out = c(out, list(phiSamps = PHI, PsiSamps = PSI, 
                           interceptSamps = INTERCEPT,
                           mainEffectSamps = MAIN,
                           interactionSamps = INTERACTION))
  if(es) out = c(out, list(ssySamps = SSY))
  
  return(out)
}
