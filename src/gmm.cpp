// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// #define NDEBUG
// #include <RcppNumerical.h>
// #include <RcppEigen.h>


// typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
// typedef Eigen::Map<Eigen::VectorXd> MapVect;


// GMM objective function
// [[Rcpp::export]]
double gmm_obj(const double&  betal,             // beta_l parameter
               const Eigen::MatrixXd& Z,         // instrument matrix
               const Eigen::VectorXd& y,         // dependent variable
               const Eigen::MatrixXd& endo,      // matrix of endo
               const Eigen::MatrixXd& X_iso,     // covariates for iso
               const Eigen::MatrixXd& X_niso,    // covariates for niso
               const Eigen::MatrixXd& W,         // weighting matrix
               const int& S) {                   // Number of subnets
  int n(endo.rows()), Kx(X_iso.cols());
  
  // 1. psi = [M, X_iso + X_niso / (1 + betal)]
  Eigen::MatrixXd psi(n, 2 + Kx);
  psi << endo, X_iso + X_niso/(1.0 + betal);
  
  // 2. Closed-form GMM parameters: 
  Eigen::MatrixXd Ztpsi(Z.transpose() * psi); // kz x kx  Z'psi
  Eigen::MatrixXd Zty(Z.transpose() * y);  // kz x 1 Z'y
  Eigen::MatrixXd ZtpsiW(Ztpsi.transpose()* W);   //  kx x kz psi'Z W
  Eigen::MatrixXd A(ZtpsiW*Ztpsi); // kx x kx psi'Z W Z' psi
  Eigen::VectorXd b(ZtpsiW * Zty); // kx x 1  psi'Z W Z'Y
  
  Eigen::VectorXd phi(A.colPivHouseholderQr().solve(b)); // (psi'Z W Z' psi)^{-1} psi'Z W Z'Y
  
  // 3. Compute moments: moments = Z * (y - psi * phi)
  Eigen::VectorXd eta(y - psi * phi);
  Eigen::VectorXd m(Z.transpose() * eta / S); // kz x 1
  
  // 4. GMM objective function: g' W g
  return m.dot(W * m);
}

// Optimal Weighting Matrix
// [[Rcpp::export]]
Eigen::MatrixXd W_optimal(const double& betal,              // beta_l parameter
                          const Eigen::MatrixXd& Z,         // instrument matrix (n x k)
                          const Eigen::VectorXd& y,         // dependent variable (n x 1)
                          const Eigen::MatrixXd& endo,      // matrix for iso group
                          const Eigen::MatrixXd& X_iso,     // covariates for iso
                          const Eigen::MatrixXd& X_niso,    // covariates for niso
                          const Eigen::MatrixXd& W,         // weighting matrix
                          const Eigen::ArrayXi& Iso,        // Indice for isolated
                          const Eigen::ArrayXi& nIso,       // Indice for nonisolated
                          const Eigen::VectorXd& cumsn,     // cumulative group indices 
                          const int& dfiso,                 // degree of freedom for isolated
                          const int& dfniso,                // degree of freedom for nonisolated
                          const int& HAC,                   // HAC type
                          const int& S) {                   // Number of subnets
  int n(endo.rows()), Kx(X_iso.cols()),  Kz(Z.cols()), n_iso(Iso.size()), 
  n_niso(n - n_iso);
  
  // 1. psi = [M, X_iso + X_niso / (1 + betal)]
  Eigen::MatrixXd psi(n, 2 + Kx);
  psi << endo, X_iso + X_niso / (1.0 + betal);
  
  // 2. Closed-form GMM parameters: phi = (psi' Z W Z' psi)^{-1} psi' Z W Z' y
  Eigen::MatrixXd Ztpsi(Z.transpose() * psi); // kz x kx  Z'psi
  Eigen::MatrixXd Zty(Z.transpose() * y);  // kz x 1 Z'y
  Eigen::MatrixXd ZtpsiW(Ztpsi.transpose()* W);   //  kx x kz psi'Z W
  Eigen::MatrixXd A(ZtpsiW*Ztpsi); // kx x kx psi'Z W Z' psi
  Eigen::VectorXd b(ZtpsiW * Zty); // kx x 1  psi'Z W Z'Y
  
  Eigen::VectorXd phi(A.colPivHouseholderQr().solve(b)); // (psi'Z W Z' psi)^{-1} psi'Z W Z'Y
  
  // 3. residual
  Eigen::VectorXd eta(y - psi * phi);
  
  // 4.Variance of the moment
  Eigen::MatrixXd Vm(Eigen::MatrixXd::Zero(Kz, Kz));
  if (HAC == 0) { //4.1 iid
    eta(nIso) *= (1 + betal);
    double s2(eta.dot(eta) / (dfiso + dfniso));
    Eigen::MatrixXd Zs2(Z * s2);
    Zs2(nIso, Eigen::all) /= pow(1 + betal, 2);
    Vm = Zs2.transpose() * Z / pow(S, 2);
  } else if (HAC == 1) {// 4.2 iid separately
    double s2iso(eta(Iso).dot(eta(Iso)) / dfiso);
    double s2niso(eta(nIso).dot(eta(nIso)) / dfniso);
    Eigen::MatrixXd Zs2(Z);
    Zs2(Iso, Eigen::all)  *= s2iso;
    Zs2(nIso, Eigen::all) *= s2niso;
    Vm = Zs2.transpose() * Z / pow(S, 2);
  } else if (HAC == 2) { // 4.3 heteroskedasticity
    Eigen::MatrixXd Z_eta(Z.array().colwise() * eta.array());
    Vm = Z_eta.transpose() * Z_eta / pow(S, 2);
  } else {
    for (int s(0); s < S; ++ s) { // 4.4 clustering
      int n1(cumsn(s)), n2(cumsn(s + 1) - 1); 
      Eigen::VectorXd Zes(Z(Eigen::seq(n1, n2), Eigen::all).transpose() * eta(Eigen::seq(n1, n2)));
      Vm += Zes * Zes.transpose();
    }
    Vm /= pow(S, 2);
  }

  // 5 Optimal weighting
  return Vm.inverse();
}


// [[Rcpp::export]]
Rcpp::List compute_estimate(const double& betal,              // beta_l parameter
                            const Eigen::MatrixXd& Z,         // instrument matrix (n x k)
                            const Eigen::VectorXd& y,         // dependent variable (n x 1)
                            const Eigen::MatrixXd& endo,      // matrix for iso group
                            const Eigen::MatrixXd& X_iso,     // covariates for iso
                            const Eigen::MatrixXd& X_niso,    // covariates for niso
                            const Eigen::MatrixXd& W,         // weighting matrix
                            const Eigen::ArrayXi& Iso,        // Indice for isolated
                            const Eigen::ArrayXi& nIso,       // Indice for nonisolated
                            const Eigen::VectorXd& cumsn,     // cumulative group indices 
                            const int& dfiso,                 // degree of freedom for isolated
                            const int& dfniso,                // degree of freedom for nonisolated
                            const int& HAC,                   // HAC type
                            const int& S) {                   // Number of subnets
  int n(endo.rows()), Kx(X_iso.cols()), 
  Kz(Z.cols()), n_iso(Iso.size()), n_niso(n - n_iso);

  // 1. psi = [M, X_iso + X_niso / (1 + betal)]
  Eigen::MatrixXd psi(n, 2 + Kx);
  psi << endo, X_iso + X_niso / (1.0 + betal);
  
  // 2. Closed-form GMM parameters: 
  Eigen::MatrixXd Ztpsi(Z.transpose() * psi); // kz x kx  Z'psi
  Eigen::MatrixXd Zty(Z.transpose() * y);  // kz x 1 Z'y
  Eigen::MatrixXd ZtpsiW(Ztpsi.transpose()* W);   //  kx x kz psi'Z W
  Eigen::MatrixXd A(ZtpsiW*Ztpsi); // kx x kx psi'Z W Z' psi
  Eigen::VectorXd b(ZtpsiW * Zty); // kx x 1  psi'Z W Z'Y
  
  // 3. estimate
  // 3.1 reduced form 
  Eigen::VectorXd phired(3 + Kx);
  phired.tail(2 + Kx) = A.colPivHouseholderQr().solve(b); // (psi'Z W Z' psi)^{-1} psi'Z W Z'Y
  phired(0)           = betal;
  Eigen::VectorXd gamma(phired.tail(Kx));
  // 3.2 structural
  Eigen::VectorXd phistr(3 + Kx); //betal, betah, delta, gamma
  phistr(0)       = phired(0); //betal
  phistr(1)       = phired(2) * (1 + phired(0)) + phired(0); // betah = theta 2 * (1 + betal) + betal
  phistr(2)       = phired(1) * (1 + phired(0)) - phired(0); // delta = theta 1 * (1 + betal) - betal
  phistr.tail(Kx) = phired.tail(Kx); // gamma
  // 4. Covariance matrice for the reduced form parameters
  // 4.0 residual
  Eigen::VectorXd eta(y - psi * phired.tail(2 + Kx));

  // 4.1 Variance of S^0.5 * moment
  Eigen::MatrixXd Vm(Eigen::MatrixXd::Zero(Kz, Kz));
  double s2(std::numeric_limits<double>::quiet_NaN());
  double s2iso(std::numeric_limits<double>::quiet_NaN());
  double s2niso(std::numeric_limits<double>::quiet_NaN());
  if (HAC == 0) { //4.1 iid
    eta(nIso) *= (1 + betal);
    s2         = eta.dot(eta) / (dfiso + dfniso);
    Eigen::MatrixXd Zs2(Z * s2);
    Zs2(nIso, Eigen::all) /= pow(1 + betal, 2);
    Vm                     = Zs2.transpose() * Z / S;
  } else if (HAC == 1) {// 4.2 iid separately
    s2iso  = eta(Iso).dot(eta(Iso)) / dfiso;
    s2niso = eta(nIso).dot(eta(nIso)) / dfniso;
    Eigen::MatrixXd Zs2(Z);
    Zs2(Iso, Eigen::all)  *= s2iso;
    Zs2(nIso, Eigen::all) *= s2niso;
    Vm      = Zs2.transpose() * Z / S;
    s2niso *= pow(1 + betal, 2);
  } else if (HAC == 2) { // 4.3 heteroskedasticity
    Eigen::MatrixXd Z_eta(Z.array().colwise() * eta.array());
    Vm = Z_eta.transpose() * Z_eta / S;
  } else {
    for (int s(0); s < S; ++ s) { // 4.4 clustering
      int n1(cumsn(s)), n2(cumsn(s + 1) - 1); 
      Eigen::VectorXd Zes(Z(Eigen::seq(n1, n2), Eigen::all).transpose() * eta(Eigen::seq(n1, n2)));
      Vm += Zes * Zes.transpose();
    }
    Vm /= S;
  }
  
  // 4.2 Jacobian G = 1/S Z'(X_niso gamma1 / (1.0 + betal)^2, -endo, -X_iso - X_niso / (1.0 + betal))
  Eigen::MatrixXd  G(Kz, 3 + Kx);
  G.col(0) = Z.transpose()* (X_niso * gamma / pow(1.0 + betal, 2)) / S;
  G(Eigen::all, Eigen::seqN(1, 2 + Kx)) = - Z.transpose() * psi / S;
  
  // 4.3 Compute variance: Var = (G'WG)^-1 when W=Omega^-1 (optimal), otherwise
  // Var = (G'WG)^-1 G'W Omega W'G (G'WG)'^-1, Omega=Vm
  Eigen::MatrixXd GtW(G.transpose() * W);
  Eigen::MatrixXd bread_inv((GtW * G).inverse());
  Eigen::MatrixXd meat(GtW * Vm * GtW.transpose());
  Eigen::MatrixXd Vred(bread_inv * meat * bread_inv.transpose() / S);
  
  // 5. Covariance matrice for the structural parameters (DELTA METHOD)
  // Vstr = D Vred D'
  Eigen::MatrixXd D = Eigen::MatrixXd::Identity(3 + Kx, 3 + Kx);
  // derivative for phistr(0) = betal
  D(0, 0) = 1.0;
  
  // derivative for phistr(1) = phired(2) * (1 + phired(0)) + phired(0)
  // d/d phired0 = phired(2) + 1
  // d/d phired2 = 1 + phired(0)
  D(1, 0) = phired(2) + 1.0;
  D(1, 1) = 0;
  D(1, 2) = 1.0 + phired(0);
  
  // derivative for phistr(2) = phired(1) * (1 + phired(0)) - phired(0)
  // d/d phired0 = phired(1) - 1
  // d/d phired1 = 1 + phired(0)
  D(2, 0) = phired(1) - 1.0;
  D(2, 1) = 1.0 + phired(0);
  D(2, 2) = 0;
  
  Eigen::MatrixXd Vstr(D * Vred * D.transpose());  
  
  // 6 J-stat (Sargan overidentification test)
  // overidentification
  Eigen::VectorXd Ze(Z.transpose() * eta.matrix());
  double stat(Ze.dot(Vm.colPivHouseholderQr().solve(Ze)) / S); // because Vm is / S
  
  return Rcpp::List::create(Rcpp::_["redparm"] = phired, // Reduced form parameter 
                            Rcpp::_["strparm"] = phistr, // Structural parameter 
                            Rcpp::_["redcov"]  = Vred,   // Reduced form covariance matrix     
                            Rcpp::_["strcov"]  = Vstr,   // Structural covariance matrix
                            Rcpp::_["s2"]      = s2,
                            Rcpp::_["s2iso"]   = s2iso,
                            Rcpp::_["s2niso"]  = s2niso,
                            Rcpp::_["JStat"]   = stat,
                            Rcpp::_["Jdf"]     = Kz - (3 + Kx));
}

