// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#ifdef _OPENMP 
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
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
  int n(endo.rows()), Kx(X_iso.cols()), Kendo(endo.cols()); 
  
  // 1. psi = [endo, X_iso + X_niso / (1 + betal)]
  Eigen::MatrixXd psi(n, Kendo + Kx);
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



// GMM objective function without spillovers 
// [[Rcpp::export]]
double gmm_obj_nospil(const double&  betal,             // beta_l parameter
               const Eigen::MatrixXd& Z,         // instrument matrix
               const Eigen::VectorXd& y,         // dependent variable
               const Eigen::MatrixXd& endo,      // matrix of endo
               const Eigen::MatrixXd& X_iso,     // covariates for iso
               const Eigen::MatrixXd& X_niso,    // covariates for niso
               const Eigen::MatrixXd& W,         // weighting matrix
               const int& S) {                   // Number of subnets
  int n(endo.rows()), Kx(X_iso.cols()), Kendo(endo.cols()); 
  double theta1(betal / (1 + betal)); // coefficient of ybar
  
  // 1. psi = [ycheck, X_iso + X_niso / (1 + betal)]
  Eigen::MatrixXd psi(n, Kendo-1 + Kx);
  if (Kendo > 1) {
  psi << endo.col(1), X_iso + X_niso/(1.0 + betal);
  } else {
    psi << X_iso + X_niso/(1.0 + betal);
  }
  // 2. Closed-form GMM parameters: 
  Eigen::MatrixXd Ztpsi(Z.transpose() * psi); // kz x (1+kx)  Z'psi
  Eigen::MatrixXd Zty(Z.transpose() * (y - theta1 * endo.col(0)));  // kz x 1 Z'(y - theta1 * ybar)
  Eigen::MatrixXd ZtpsiW(Ztpsi.transpose()* W);   //  (1+kx) x kz psi'Z W
  Eigen::MatrixXd A(ZtpsiW*Ztpsi); // (1+kx) x (1+kx) psi'Z W Z' psi
  Eigen::VectorXd b(ZtpsiW * Zty); // (1+kx) x 1  psi'Z W Z'Y
  
  Eigen::VectorXd phi(A.colPivHouseholderQr().solve(b)); // (psi'Z W Z' psi)^{-1} psi'Z W Z'(Y - theta1 * ybar)
  
  // 3. Compute moments: moments = Z * (y - theta1 * ybar - psi * phi)
  Eigen::VectorXd eta(y - theta1 * endo.col(0) - psi * phi);
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
  int n(endo.rows()), Kx(X_iso.cols()),  Kz(Z.cols()), Kendo(endo.cols());
  
  // 1. psi = [endo, X_iso + X_niso / (1 + betal)]
  Eigen::MatrixXd psi(n, Kendo + Kx);
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


// Optimal Weighting Matrix
// [[Rcpp::export]]
Eigen::MatrixXd W_optimal_nospil(const double& betal,              // beta_l parameter
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
  int n(endo.rows()), Kx(X_iso.cols()),Kz(Z.cols()), Kendo(endo.cols()); 
  double theta1(betal / (1 + betal)); // coefficient of ybar
  
  // 1. psi = [ycheck, X_iso + X_niso / (1 + betal)]
  Eigen::MatrixXd psi(n, Kendo-1 + Kx);
  if (Kendo > 1) {
    psi << endo.col(1), X_iso + X_niso/(1.0 + betal);
  } else {
    psi << X_iso + X_niso/(1.0 + betal);
  }
  // 2. Closed-form GMM parameters: 
  Eigen::MatrixXd Ztpsi(Z.transpose() * psi); // kz x kx  Z'psi
  Eigen::MatrixXd Zty(Z.transpose() * (y - theta1 * endo.col(0)));  // kz x 1 Z'(y - theta1 * ybar)
  Eigen::MatrixXd ZtpsiW(Ztpsi.transpose()* W);   //  kx x kz psi'Z W
  Eigen::MatrixXd A(ZtpsiW*Ztpsi); // kx x kx psi'Z W Z' psi
  Eigen::VectorXd b(ZtpsiW * Zty); // kx x 1  psi'Z W Z'Y
  
  Eigen::VectorXd phi(A.colPivHouseholderQr().solve(b)); // (psi'Z W Z' psi)^{-1} psi'Z W Z'(Y - theta1 * ybar)
  
  // 3. residual
  Eigen::VectorXd eta(y - theta1 * endo.col(0) - psi * phi);
  
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


// Compute parameters and variance, with spillover
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
  int n(endo.rows()), Kx(X_iso.cols()), Kz(Z.cols()), Kendo(endo.cols());
  
  // 1. psi = [M, X_iso + X_niso / (1 + betal)]
  Eigen::MatrixXd psi(n, Kendo + Kx);
  psi << endo, X_iso + X_niso / (1.0 + betal);
  
  // 2. Closed-form GMM parameters: 
  Eigen::MatrixXd Ztpsi(Z.transpose() * psi); // kz x kx  Z'psi
  Eigen::MatrixXd Zty(Z.transpose() * y);  // kz x 1 Z'y
  Eigen::MatrixXd ZtpsiW(Ztpsi.transpose()* W);   //  kx x kz psi'Z W
  Eigen::MatrixXd A(ZtpsiW*Ztpsi); // kx x kx psi'Z W Z' psi
  Eigen::VectorXd b(ZtpsiW * Zty); // kx x 1  psi'Z W Z'Y
  
  // 3. estimate
  // 3.1 reduced form  [betal, theta1, theta2, gamma]
  Eigen::VectorXd phired(1 + Kendo + Kx);
  phired.tail(Kendo + Kx) = A.colPivHouseholderQr().solve(b); // (psi'Z W Z' psi)^{-1} psi'Z W Z'Y
  phired(0)           = betal;
  Eigen::VectorXd gamma(phired.tail(Kx));
  // 3.2 structural [betal, betah, delta, gamma]
  Eigen::VectorXd phistr(1 + Kendo + Kx);
  if(Kendo>1){
  phistr(0)       = phired(0); //betal
  phistr(1)       = phired(2) * (1 + phired(0)) + phired(0); // betah = theta 2 * (1 + betal) + betal
  phistr(2)       = phired(1) * (1 + phired(0)) - phired(0); // delta = theta 1 * (1 + betal) - betal
  } else {
    phistr(0)       = phired(0); //betal
    phistr(1)       = phired(1) * (1 + phired(0)) - phired(0); // delta = theta 1 * (1 + betal) - betal
  }
  phistr.tail(Kx) = phired.tail(Kx); // gamma
  // 4. Covariance matrice for the reduced form parameters
  // 4.0 residual
  Eigen::VectorXd eta(y - psi * phired.tail(Kendo + Kx));
  
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
  Eigen::MatrixXd  G(Kz, 1 + Kendo + Kx);
  G.col(0) = Z.transpose()* (X_niso * gamma / pow(1.0 + betal, 2)) / S;
  G(Eigen::all, Eigen::seqN(1, Kendo + Kx)) = - Z.transpose() * psi / S;
  
  // 4.3 Compute variance: Var = (G'WG)^-1 when W=Omega^-1 (optimal), otherwise
  // Var = (G'WG)^-1 G'W Omega W'G (G'WG)'^-1, Omega=Vm
  Eigen::MatrixXd GtW(G.transpose() * W); 
  Eigen::MatrixXd bread_inv((GtW * G).inverse());
  Eigen::MatrixXd meat(GtW * Vm * GtW.transpose());
  Eigen::MatrixXd Vred(bread_inv * meat * bread_inv.transpose() / S);
  
  // 5. Covariance matrice for the structural parameters (DELTA METHOD)
  // Vstr = D Vred D'
  Eigen::MatrixXd D = Eigen::MatrixXd::Identity(1 + Kendo + Kx, 1 + Kendo + Kx);
  // derivative for phistr(0) = betal
  D(0, 0) = 1.0;
  
  if(Kendo>1){
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
  } else {
  // derivative for  phistr(1) = phired(1) * (1 + phired(0)) - phired(0)
  // d/d phired0 = phired(1) - 1
  // d/d phired1 = 1 + phired(0)
  D(1, 0) = phired(1) - 1.0;
  D(1, 1) = 1.0 + phired(0);
  }
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


// Compute parameters and variance, without spillover
// [[Rcpp::export]]
Rcpp::List compute_estimate_nospil(const double& betal,              // beta_l parameter
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
  int n(endo.rows()), Kx(X_iso.cols()), Kz(Z.cols()), Kendo(endo.cols());
  double theta1(betal / (1 + betal)); // coefficient of ybar
  
  // 1. psi = [ycheck, X_iso + X_niso / (1 + betal)]
  Eigen::MatrixXd psi(n, Kendo-1 + Kx);
  if (Kendo > 1) {
    psi << endo.col(1), X_iso + X_niso/(1.0 + betal);
  } else {
    psi << X_iso + X_niso/(1.0 + betal);
  }
  
  // 2. Closed-form GMM parameters: 
  Eigen::MatrixXd Ztpsi(Z.transpose() * psi); // kz x kx  Z'psi
  Eigen::MatrixXd Zty(Z.transpose() * (y - theta1 * endo.col(0)));  // kz x 1 Z'(y - theta1 * ybar)
  Eigen::MatrixXd ZtpsiW(Ztpsi.transpose()* W);   //  kx x kz psi'Z W
  Eigen::MatrixXd A(ZtpsiW*Ztpsi); // kx x kx psi'Z W Z' psi
  Eigen::VectorXd b(ZtpsiW * Zty); // kx x 1  psi'Z W Z'Y
  
  // 3. estimate
  // 3.1 reduced form  [theta1, theta2, gamma]
  Eigen::VectorXd phired(Kendo + Kx);
  phired.tail(Kendo - 1 + Kx) = A.colPivHouseholderQr().solve(b); // (psi'Z W Z' psi)^{-1} psi'Z W Z'(Y - theta1 * Ybar)
  phired(0)           = theta1;
  Eigen::VectorXd gamma(phired.tail(Kx));
  // 3.2 structural [betal, betah, gamma]
  Eigen::VectorXd phistr(Kendo + Kx); 
  phistr(0)       = betal; //betal
  if (Kendo > 1) {
  phistr(1)       = phired(1) * (1 + betal) + betal; // betah = theta 2 * (1 + betal) + betal
  } 
  phistr.tail(Kx) = gamma; // gamma
  // 4. Covariance matrice for the reduced form parameters
  // 4.0 residual
  Eigen::VectorXd eta(y - phired(0) * endo.col(0) - psi * phired.tail(Kendo -1 + Kx));
  
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
  Eigen::MatrixXd  G(Kz, Kendo + Kx);
  G.col(0) = Z.transpose()* (X_niso * gamma / pow(1.0 + betal, 2)) / S;
  G(Eigen::all, Eigen::seqN(1, Kendo -1 + Kx)) = - Z.transpose() * psi / S;
  
  // 4.3 Compute variance: Var = (G'WG)^-1 when W=Omega^-1 (optimal), otherwise
  // Var = (G'WG)^-1 G'W Omega W'G (G'WG)'^-1, Omega=Vm
  Eigen::MatrixXd GtW(G.transpose() * W); // ((2+Kx) x Kz)x  (Kz x Kz) = (2+Kx) x Kz
  Eigen::MatrixXd bread_inv((GtW * G).inverse()); // ((2+Kx) x Kz) x (Kz x (2+Kx)) = (2+Kx) x (2+Kx)
  Eigen::MatrixXd meat(GtW * Vm * GtW.transpose()); // (2+Kx) x Kz * Kz x Kz * Kz x (2+Kx) = (2+Kx) x (2+Kx)
  Eigen::MatrixXd Vred(bread_inv * meat * bread_inv.transpose() / S); // (2+Kx) x (2+Kx)
  
  // 5. Covariance matrice for the structural parameters (DELTA METHOD)
  // Vstr = D Vred D'
  Eigen::MatrixXd D = Eigen::MatrixXd::Identity(Kendo + Kx, Kendo + Kx);
  // derivative for phistr(0) = betal
  D(0, 0) = 1.0;
  if(Kendo>1){
  // derivative for phistr(1)  = phired(1) * (1 + betal) + betal
  // d/d betal = phired(1) + 1
  // d/d phired1 = 1 + phired(0)
  D(1, 0) = phired(1) + 1.0;
  D(1, 1) = 1.0 + phired(0);
  }
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