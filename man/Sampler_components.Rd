\name{Sampler Components}
\alias{del_mg}
\alias{eta_int}
\alias{eta_lin}
\alias{lam_lin}
\alias{mh}
\alias{phi_dl}
\alias{phi_int}
\alias{plm_dl}
\alias{plm_mg}
\alias{psi_dl}
\alias{psi_int}
\alias{psi_mg}
\alias{rgig}
\alias{rig}
\alias{sig_lin}
\alias{ssy_int}
\alias{tau_dl}

\title{
Sampler Components
}
\description{
These are the component full conditional or Metropolis-Hastings updates coded in C++ used in the samplers in this package. The functions follow naming conventions based on their greek letter notation in their respective original papers, cited below, and the paper they come from. Here _mg refers to a component of the Multiplicative Gamma Shrinkage prior of Bhattacharya and Dunson 2011, _dl refers to a component of the Dirichlet-Laplace shrinkage prior of Bhattacharya et al., _lin refers to a component of a linear factor model as in Bhattacharya and Dunson 2011, and _int refers to a component of a factor model with 2-way interactions as in Ferrari and Dunson 2020.
}
\references{
Bhattacharya, Anirban, and David B. Dunson. "Sparse Bayesian infinite factor models." Biometrika (2011): 291-306.

Bhattacharya, Anirban, et al. "Dirichlet-Laplace priors for optimal shrinkage." Journal of the American Statistical Association 110.512 (2015): 1479-1490.

Ferrari, Federico, and David B. Dunson. "Bayesian Factor Analysis for Inference on Interactions." arXiv preprint arXiv:1904.11603 (2019).
}
\author{
Evan Poworoznek
}
