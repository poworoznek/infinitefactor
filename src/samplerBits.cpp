#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;


// [[Rcpp::export]]
Rcpp::NumericMatrix eta_lin(arma::mat lambda, arma::vec ps, int k, int n, arma::mat Y){
  // --- UPDATE eta --- //
  arma::mat Lmsg = lambda.each_col() % ps;
  arma::mat Veta1 = eye<arma::mat>(k,k) + Lmsg.t() * lambda;
  arma::mat S = inv(trimatu(chol(Veta1)));
  arma::mat Veta = S * S.t();
  arma::mat Meta = Y * Lmsg * Veta;
  arma::mat noise = randn<arma::mat>(n, k);
  arma::mat eta = Meta + noise * S.t();
  return Rcpp::wrap(eta);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix lam_lin(arma::mat eta, arma::mat Plam, arma::vec ps, int k, int p, arma::mat Y){
  // --- UPDATE lambda --- //
  arma::mat lambda(p, k);
  arma::mat eta2 = eta.t() * eta;    // prepare eta crossproduct before the loop
  for(int j=0; j < p; ++j) {
    arma::mat Llamt = trimatu(chol(diagmat(Plam.row(j)) + ps(j)*eta2));
    arma::mat Llam = trimatl(Llamt.t());
    lambda.row(j) = (solve(Llamt, randn<arma::vec>(k)) +
      solve(Llamt, solve(Llam, ps(j) * eta.t() * Y.col(j)))).t();
  }
  return Rcpp::wrap(lambda);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix psi_mg(arma::mat lambda, arma::vec tauh, arma::vec ps, int k, int p, double df){
  // --- UPDATE psihj --- //
  arma::mat psijh(p, k);
  arma::mat lambda_sq = square(lambda);
  arma::mat shape = lambda_sq.each_row() % tauh.t();
  for (int l=0; l < p; l++) {
    for (int j=0; j < k; j++) {
      psijh(l,j) = randg(distr_param(df/2 + 0.5, 1 / (df/2 + shape(l,j)/2)));
    }
  }
  return Rcpp::wrap(psijh);
}




// [[Rcpp::export]]
Rcpp::NumericVector del_mg(arma::mat lambda, arma::mat psijh, arma::vec tauh, arma::vec delta, int k, int p, 
                               double ad1, double bd1, double ad2, double bd2){
  // --- UPDATE delta & tauh --- //
  arma::mat matr = psijh % square(lambda);
  double ad = ad1 + 0.5 * p * k;
  double bd = bd1 + 0.5 / delta(0) * sum(tauh.t() % sum(matr, 0));
  delta(0) = randg(distr_param(ad, 1 / bd));
  tauh = cumprod(delta);
  
  for(int h=1; h<k; ++h) {
    ad = ad2 + 0.5 * p *(k-h);
    arma::vec tauh_sub = tauh.subvec(h,k-1);
    bd = bd2 + 0.5 / delta(h) * sum(tauh_sub.t() % sum(matr.cols(h,k-1), 0));
    delta(h) = randg(distr_param(ad, 1 / bd));
    tauh = cumprod(delta);
  }
  return Rcpp::wrap(delta);
}

// [[Rcpp::export]]
Rcpp::NumericVector sig_lin(arma::mat lambda, arma::mat eta, int k, int p, int n, arma::mat Y,
                               double as, double bs){
  // --- UPDATE sigma --- //
  arma::mat Ytil = Y - eta * lambda.t();
  rowvec bsvec =  bs + 0.5 * sum(square(Ytil));
  auto lam = [as, n](double val){return randg<double>(distr_param(as + 0.5*n, 1 / val));};
  return Rcpp::wrap((bsvec.transform(lam)).t());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix plm_mg(arma::mat psijh, arma::vec tauh){
  // --- UPDATE Plam --- //
  return Rcpp::wrap(psijh.each_row() % tauh.t());
}

// [[Rcpp::export]]
double rig(double mu){
  double y = randn<double>();
  y *= y;
  double mu2 = std::pow(mu, 2);
  double quad = 4 * mu * y + mu2 * std::pow(y, 2);
  double x = mu + y * mu2 / 2 - mu / 2  * std::pow(quad, 0.5);
  double u = randu<double>();
  if(u < (mu / (x + mu))) return x;
  else return mu2 / x;
}

 // [[Rcpp::export]]
 double rgig(double lam, double chi, double psi){
   double omega = std::pow(psi*chi, 0.5);
   double alpha = std::pow(chi/psi, 0.5);
   bool flip = false;
   long mode;
   if(lam<0) {lam *= -1; flip = true;}
   if(lam>=1) mode = (std::pow((lam-1)*(lam-1) + omega*omega, 0.5)+(lam-1))/omega;
   else mode = omega/(std::pow((lam-1)*(lam-1) + omega*omega, 0.5)+(lam-1));
   double lm = std::log(mode);
   double cons = 0.5 * (lam-1) * lm - 0.25*omega*(mode + 1/mode);
   long maxp = ((lam+1)+std::pow((lam+1)*(lam+1)+omega*omega,0.5))/omega;
   long eq = 0.5*(lam+1)*std::log(maxp) - 0.25*omega*(maxp + 1/maxp) - cons;
   double ext = std::exp(eq);
   float U, L, prop;
   do{
     U = ext * randu<float>();
     L = randu<float>();
     prop = U/L;
   } while (
       ((std::log(L)) > 
          (0.5*(lam-1)*std::log(prop) - 
          0.25*omega*(prop + 1/prop) - 
          cons)));
   if(flip) return alpha/prop;
   else return prop/alpha;
 }


// [[Rcpp::export]]
Rcpp::NumericMatrix psi_dl(arma::mat lambda, arma::mat phi, arma::vec tau){
  // --- UPDATE psihj --- //
  arma::mat means = phi % pow(abs(lambda), -1);
  means.each_col() %= tau;
  return Rcpp::wrap(means.transform(rig));
}

// [[Rcpp::export]]
Rcpp::NumericVector tau_dl(arma::mat lambda, arma::mat phi, int k, int p){
  arma::vec chis = sum(abs(lambda) % pow(phi, -1), 1);
  arma::vec tau(p);
  for(int i=0;i<p;++i){
    tau(i) = rgig(1-k, 2*chis(i), 1);
  }
  return Rcpp::wrap(tau);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix phi_dl(arma::mat lambda, double a, int k, int p){
  arma::mat chis = abs(lambda);
  arma::mat phi(p,k);
  rowvec ts(k);
  for(int i=0;i<p;++i){
    for(int j=0;j<k;++j){
      ts(j) = rgig(a-1, 2*chis(i,j), 1);
    }
    ts *= 1/accu(ts);
    phi.row(i)=ts;
  }
  return Rcpp::wrap(phi);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix plm_dl(arma::mat psi, arma::mat phi, arma::vec tau){
  // --- UPDATE Plam --- //
  arma::mat prod = psi % square(phi);
  arma::mat prod2 = prod.each_col() % square(tau);
  return Rcpp::wrap(pow(prod2,1));
}

// [[Rcpp::export]]
double mh(arma::vec trial, arma::mat amh, arma::mat lambda, arma::mat psi, arma::vec phi, arma::vec ps, arma::vec xr, double yr, double ssy){
  arma::mat tPt = trial.t() * psi * trial;
  arma::mat res = trial.t() * (amh - 2 * psi * (yr/ssy)) * trial -
    2 * trial.t() * (lambda.t() * diagmat(ps) * xr + phi * (yr/ssy)) +
    2 * trial.t() * phi * tPt(0,0) / ssy +
    square(tPt) / ssy;
  return res(0,0);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix eta_int(arma::mat lambda, arma::mat eta, arma::vec ps, arma::vec phi, arma::mat Psi, 
                            int k, int n, arma::vec Y, arma::mat X, double sigmasq_y, double delta_rw,
                            Rcpp::NumericVector acp){
  // --- UPDATE eta --- //
  arma::mat aMH = phi * phi.t() * (1 / sigmasq_y) + lambda.t()*diagmat(ps)*lambda + diagmat(ones<arma::vec>(k));
  rowvec prop(k);
  arma::vec noise(k);
  double logr, logu;
  
  for(int i=0;i<n;++i){
    prop = eta.row(i) + randn<rowvec>(k)*delta_rw;
    
    logr = mh(prop.t(), aMH, lambda, Psi, phi, ps, (X.row(i)).t(), Y(i), sigmasq_y) -
      mh((eta.row(i)).t(), aMH, lambda, Psi, phi, ps, (X.row(i)).t(), Y(i), sigmasq_y);
    logr *= -0.5;
    
    logu = std::log(randu<float>());
    
    if(logr>logu){
      eta.row(i) = prop;
      acp[i] += 1;
    }
  }
  return Rcpp::wrap(eta);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix psi_int(arma::mat eta, arma::vec y, arma::vec phi, double ssy, int k, int n){
  arma::mat etacopy = eta;
  arma::mat Xreg = square(eta);
  for(int i=0;i<(k-1);++i){
    arma::mat add = eta.each_col() % eta.col(0);
    add.shed_col(0);
    Xreg.insert_cols(Xreg.n_cols, add);
    eta.shed_col(0);
  }
  
  arma::mat ln = Xreg.t() * Xreg / ssy;
  ln.diag() += 1;
  arma::mat s = inv(trimatu(chol(ln)));
  arma::mat v = s * s.t();
  arma::vec m = v * Xreg.t() * (y - etacopy * phi) / ssy;
  arma::vec csi = m + (randn<rowvec>(k*(k-1)/2 + k) * s.t()).t();
  
  arma::mat psi(k, k, fill::zeros);
  int ind = k;
  for(int j=0; j<(k-1); ++j){
    for(int i=(j+1); i<k; ++i){
      psi(i, j) = csi(ind) / 2;
      ind +=1;
    }
  }
  
  psi += psi.t();
  psi.diag() = csi.subvec(0,k-1);
  
  if(k==1){psi(0,0) = csi(0);}
  
  return Rcpp::wrap(psi);
}


// [[Rcpp::export]]
Rcpp::NumericVector phi_int(arma::mat eta, arma::vec y, double ssy, arma::mat psi, int k){
  arma::mat ln = eta.t() * eta / ssy;
  ln.diag() += 0.5;
  arma::mat s = inv(trimatu(chol(ln)));
  arma::mat v = s * s.t();
  arma::vec m = v * eta.t() * (y - diagvec(eta * psi * eta.t())) / ssy;
  return Rcpp::wrap(m + (randn<rowvec>(k) * s.t()).t());
}

// [[Rcpp::export]]
double ssy_int(arma::mat eta, arma::vec phi, arma::mat psi, arma::vec y, int n){
  double an = 0.5 + n/2;
  arma::mat resid = y - eta * phi - diagvec(eta * psi * eta.t());
  arma::mat ss = resid.t() * resid;
  double bn = 0.5 + 0.5 * ss(0,0);
  return 1 / randg<double>(distr_param(an, 1/bn));
}



