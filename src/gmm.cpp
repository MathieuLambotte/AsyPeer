// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
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

///////////////////////////////////////////////////////////////////////
////////////////////////////// Structures /////////////////////////////
///////////////////////////////////////////////////////////////////////
struct strEstim { // output for fAsyGmmEstim and fAsyGmmEstim_nospill
  Eigen::VectorXd theta;
  Eigen::VectorXd eta;
  Eigen::MatrixXd sV;
  Eigen::VectorXd sVphi;
  double objective;
  Eigen::VectorXd gradient;
  int  status;
}; 

struct strCov { // output for fAsyparms and fAsyparms_nospill
  Eigen::VectorXd estimate;
  Eigen::MatrixXd cov;
  double serr;
  double serriso;
  double serrniso;
  double Jstat;
  double Jdf;
  Eigen::ArrayXd TestAsym;                
};

////////////////////////////////////////////////////////////////////////
//////////////// A class for the numerical optimization ////////////////
////////////// It will be used for models with spillover ///////////////
////////////////////////////////////////////////////////////////////////
class ClassAsyGMM: public Numer::MFuncGrad {
private:
  const Eigen::MatrixXd& Z;   // instrument matrix
  const Eigen::VectorXd& y;   // dependent variable
  const Eigen::MatrixXd& V;   // matrix V = [endo, X, X_iso_noncommon]
  const Eigen::MatrixXd& W;   // weighting matrix
  const Eigen::ArrayXi& nIso; // Indices nonisolates
  const int& n;               // Sample size
  const int& Ktheta;          // Number of parameters
public:
  ClassAsyGMM(const Eigen::MatrixXd& Z_,
              const Eigen::VectorXd& y_,
              const Eigen::MatrixXd& V_,
              const Eigen::MatrixXd& W_,
              const Eigen::ArrayXi& nIso_,
              const int& n_,
              const int& Ktheta_) :
  Z(Z_),
  y(y_),
  V(V_),
  W(W_),
  nIso(nIso_),
  n(n_),
  Ktheta(Ktheta_){}
  
  Eigen::VectorXd eta;   // residual will be exported
  Eigen::VectorXd theta; // the full reduced-form parameter
  Eigen::MatrixXd sV;    // Scale V (for non-isolated)
  Eigen::VectorXd sVphi; // sV * phi
  Eigen::VectorXd Grad;  // Gradient to be exported
  
  // function computing the objective function and gradient
  double f_grad(Numer::Constvec& tbetal, Numer::Refvec grad) {
    // 0. betal
    double betal  = tbetal(0);
    if (betal <= -0.5) {
      grad << 1e300;
      Grad = grad; // exported
      return 1e300;
    }
    
    // 1. Scale V for non-isolated
    sV = V;
    sV(nIso, Eigen::all) /= (1 + betal);
    
    // 2. Closed-form GMM parameters:
    Eigen::MatrixXd ZtsV(Z.transpose() * sV); // kz x ksV: Z'sV
    Eigen::VectorXd Zty(Z.transpose() * y);  // kz x 1: Z'y
    Eigen::MatrixXd sVtZW(ZtsV.transpose() * W);   //  ksV x kz: sV'Z W
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Adec(sVtZW * ZtsV); // cholesky decomposition for A = sV'Z W Z' sV
    Eigen::VectorXd b(sVtZW * Zty); // ksV x 1: sV'Z W Z'y
    Eigen::VectorXd phi(Adec.solve(b));
    
    // Full parameter
    theta.resize(Ktheta);
    theta << betal, phi;
    
    // 3. residual (but scaled for nonisolated)
    sVphi = sV * phi;
    eta   = y - sVphi;
    
    // 4. Objective function
    Eigen::VectorXd mom(Z.transpose() * eta);
    Eigen::VectorXd Wmom = W * mom;
    double f = mom.dot(Wmom);
    
    // 5. Gradient with respect to betal
    // 5.1 dV: derivative of V
    Eigen::MatrixXd dsV(Eigen::MatrixXd::Zero(n, Ktheta - 1)); 
    dsV(nIso, Eigen::all) = -sV(nIso, Eigen::all) / (1 + betal); // dsV/dbetal
    
    // 5.2 dphi: derivative of phi
    Eigen::MatrixXd ZtdsV(Z.transpose() * dsV); // kz x ksV: Z'dsV
    Eigen::VectorXd ZtdsVphi(ZtdsV * phi);
    Eigen::MatrixXd WZteta(W * Z.transpose() * eta); // kz x 1: W Z' eta
    Eigen::VectorXd dphi = Adec.solve(ZtdsV.transpose() * WZteta - sVtZW * ZtdsVphi);
    grad = 2 * Wmom.transpose() * (-ZtdsVphi - ZtsV * dphi);
    Grad = grad; // exported
    
    return f;
  }
};


////////////////////////////////////////////////////////////////////////
///////////////// Computes the reduced form parameters /////////////////
/////////////////////////// And residuals //////////////////////////////
////////////////////////////////////////////////////////////////////////
strEstim fAsyGmmEstim(const double betal,             // starting value for betal
                      const Eigen::MatrixXd& Z,       // instrument matrix
                      const Eigen::VectorXd& y,       // dependent variable
                      const Eigen::MatrixXd& V,       // Model's variables
                      const Eigen::MatrixXd& W,       // weighting matrix
                      const Eigen::ArrayXi& nIso,     // Indices for nonisolated
                      const int& n,                   // Sample size
                      const int& S,                   // Number of subnetworks
                      const int& Ktheta,              // Number of parameters
                      const int& maxit,               // optimizer controls
                      const double& eps_f,
                      const double& eps_g){
  
  // 1. Optimization
  double fopt; // objective
  Eigen::VectorXd tbetal(1);
  tbetal << log(betal + 0.5); // optimizer tbetal
  ClassAsyGMM f(Z, y, V, W, nIso, n, Ktheta);
  int status(optim_lbfgs(f, tbetal, fopt, maxit, eps_f, eps_g));
  
  // 2. Output. Put element in the structure Estim
  strEstim out;
  out.theta     = f.theta;
  out.eta       = f.eta;
  out.sV        = f.sV;
  out.sVphi     = f.sVphi;
  out.objective = fopt / pow(S, 2);
  out.gradient  = f.Grad;
  out.status    = status;
  
  return out;
}


////////////////////////////////////////////////////////////////////////
//////////////// A class for the numerical optimization ////////////////
//////////// It will be used for models without spillover //////////////
////////////////////////////////////////////////////////////////////////
class ClassAsyGMM_nospill: public Numer::MFuncGrad {
private:
  const Eigen::MatrixXd& Z;   // instrument matrix
  const Eigen::VectorXd& y;   // dependent variable
  const Eigen::VectorXd& Gy;  // Gy
  const Eigen::MatrixXd& V;   // matrix V = [endo, X, X_iso_noncommon]
  const Eigen::MatrixXd& W;   // weighting matrix
  const Eigen::ArrayXi& nIso; // Indices nonisolates
  const int& n;               // Sample size
  const int& Ktheta;          // Number of parameters 
public:
  ClassAsyGMM_nospill(const Eigen::MatrixXd& Z_,         
                      const Eigen::VectorXd& y_,    
                      const Eigen::VectorXd& Gy_,  
                      const Eigen::MatrixXd& V_,        
                      const Eigen::MatrixXd& W_,  
                      const Eigen::ArrayXi& nIso_,
                      const int& n_,
                      const int& Ktheta_) : 
  Z(Z_),        
  y(y_),    
  Gy(Gy_),
  V(V_),        
  W(W_),
  nIso(nIso_),
  n(n_),
  Ktheta(Ktheta_){}
  
  Eigen::VectorXd eta;   // residual will be exported
  Eigen::VectorXd theta; // the full reduced-form parameter
  Eigen::MatrixXd sV;    // Scale V (for non-isolated)
  Eigen::VectorXd sVphi; // sV * phi
  Eigen::VectorXd Grad;  // Gradient to be exported
  
  // function computing the objective function and gradient
  double f_grad(Numer::Constvec& tbetal, Numer::Refvec grad) {
    // 0. betal and theta0 (coefficient of Gy)
    double betal  = tbetal(0);
    double theta0 = betal / (1.0 + betal);
    if (betal <= -0.5) {
      grad << 1e300;
      Grad = grad; // exported
      return 1e300;
    }
    
    // 1. Scale V for non-isolated
    sV = V;
    sV(nIso, Eigen::all) /= (1.0 + betal);
    
    // 2. Closed-form GMM parameters: 
    Eigen::VectorXd u = y - theta0 * Gy;
    Eigen::MatrixXd ZtsV(Z.transpose() * sV); // kz x ksV: Z'sV
    Eigen::VectorXd Ztu(Z.transpose() * u);  // kz x 1: Z'u
    Eigen::MatrixXd sVtZW(ZtsV.transpose() * W);   //  ksV x kz: sV'Z W
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Adec(sVtZW * ZtsV); // cholesky decomposition for A = sV'Z W Z' sV
    Eigen::VectorXd b(sVtZW * Ztu); // ksV x 1: sV'Z W Z'u
    Eigen::VectorXd phi(Adec.solve(b));
    // Full parameter
    theta.resize(Ktheta);
    theta << betal, phi; 
    
    // 3. residual (but scaled for nonisolated)
    sVphi = sV * phi;
    eta = u - sVphi;
    // 4. Objective function
    Eigen::VectorXd mom(Z.transpose() * eta);
    Eigen::VectorXd Wmom = W * mom;
    double f = mom.dot(Wmom);
    
    // 5. Gradient with respect to betal
    // 5.1 du: derivative of u
    Eigen::VectorXd du(-Gy / pow(1.0 + betal, 2));
    
    // 5.2 dV: derivative of V
    Eigen::MatrixXd dsV(Eigen::MatrixXd::Zero(n, Ktheta - 1)); 
    dsV(nIso, Eigen::all) = -sV(nIso, Eigen::all) / (1 + betal); // dsV/dbetal
    
    // 5.3 dphi: derivative of phi
    Eigen::MatrixXd ZtdsV(Z.transpose() * dsV); // kz x ksV: Z'dsV
    Eigen::VectorXd ZtdudsV((Z.transpose() * du) - (ZtdsV * phi)); // Z' (du - dsV * phi)
    Eigen::VectorXd WZteta(W * (Z.transpose() * eta)); // kz x 1: W Z' eta
    Eigen::VectorXd dphi = Adec.solve(ZtdsV.transpose() * WZteta + sVtZW * ZtdudsV);
    
    grad = 2 * Wmom.transpose() * (ZtdudsV - ZtsV * dphi);
    Grad = grad; // exported
    
    return f;
  }
};

////////////////////////////////////////////////////////////////////////
///////////////// Computes the reduced form parameters /////////////////
//////////////////// And residuals without spillovers///////////////////
////////////////////////////////////////////////////////////////////////
strEstim fAsyGmmEstim_nospill(const double betal,             // starting value for betal
                              const Eigen::MatrixXd& Z,       // instrument matrix
                              const Eigen::VectorXd& y,       // dependent variable
                              const Eigen::VectorXd& Gy,      // Average peer friend
                              const Eigen::MatrixXd& V,       // Model's variables
                              const Eigen::MatrixXd& W,       // weighting matrix
                              const Eigen::ArrayXi& nIso,     // Indices for nonisolated
                              const int& n,                   // Sample size
                              const int& S,                   // Number of subnetworks
                              const int& Ktheta,              // Number of parameters
                              const int& maxit,               // optimizer controls
                              const double& eps_f,
                              const double& eps_g){      
  
  // 1. Optimization 
  double fopt; // objective
  Eigen::VectorXd tbetal(1);
  tbetal << log(betal + 0.5); // optimizer tbetal
  ClassAsyGMM_nospill f(Z, y, Gy, V, W, nIso, n, Ktheta);
  int status(optim_lbfgs(f, tbetal, fopt, maxit, eps_f, eps_g));
  
  // 2. Output. Put element in the structure Estim
  strEstim out;
  out.theta     = f.theta;
  out.eta       = f.eta;
  out.sV        = f.sV;
  out.sVphi     = f.sVphi;
  out.objective = fopt / pow(S, 2);
  out.gradient  = f.Grad;
  out.status    = status;
  
  return out;
}


////////////////////////////////////////////////////////////////////////
/////////////////// Computes the optimal GMM weight ////////////////////
////////////////////////////////////////////////////////////////////////
Eigen::MatrixXd fAsyWopt(const Eigen::VectorXd& theta,     // Parameters
                         const Eigen::MatrixXd& Z,         // instrument matrix (n x k)
                         const Eigen::VectorXd& eta,       // residuals
                         const Eigen::ArrayXi& cumsn,     // cumulative group indices 
                         const Eigen::ArrayXi& Iso,        // Indice for isolated
                         const Eigen::ArrayXi& nIso,       // Indice for nonisolated
                         const int& dfiso,                 // degree of freedom for isolated
                         const int& dfniso,                // degree of freedom for nonisolated
                         const int& HAC,                   // HAC type
                         const int& S) {                   // Number of subnets
  int Kz(Z.cols());
  // 0. unscaled residuals
  Eigen::VectorXd uneta(eta);
  uneta(nIso) *= (1 + theta(0));
  
  // 1. Variance of the moment
  Eigen::MatrixXd Vm(Eigen::MatrixXd::Zero(Kz, Kz));
  if (HAC == 0) { //1.1 iid
    
    double serr(sqrt(uneta.dot(uneta) / (dfiso + dfniso)));
    Eigen::MatrixXd Ze = Z * serr;
    Ze(nIso, Eigen::all) /= (1 + theta(0));
    Vm = Ze.transpose() * Ze / pow(S, 2); 
    
  } else if (HAC == 1) {// 1.2 iid separately
    
    double serriso(sqrt(eta(Iso).dot(eta(Iso)) / dfiso));
    double serrniso(sqrt(eta(nIso).dot(eta(nIso)) / dfniso));
    Eigen::MatrixXd Ze = Z;
    Ze(Iso, Eigen::all)  *= serriso;
    Ze(nIso, Eigen::all) *= serrniso;
    Vm = Ze.transpose() * Ze / pow(S, 2);
    
  } else if (HAC == 2) { // 1.3 heteroskedasticity
    
    Eigen::MatrixXd Ze    = Z.array().colwise() * eta.array();
    Vm = Ze.transpose() * Ze / pow(S, 2);
    
  } else {
    
    Eigen::ArrayXXd Ze    = Z.array().colwise() * eta.array();
    for (int s(0); s < S; ++ s) { // 1.4 clustering
      int n1(cumsn(s)), n2(cumsn(s + 1) - 1); 
      Eigen::RowVectorXd Zes(Ze(Eigen::seq(n1, n2), Eigen::all).colwise().sum());
      Vm += Zes.transpose() * Zes;
    }
    Vm /= pow(S, 2);
    
  }
  
  // 2. Optimal weighting
  return Vm.inverse();
}



////////////////////////////////////////////////////////////////////////
//////////////// Computes the structural parameters ////////////////////
/////////////////////// and Asymmetric Variance ////////////////////////
////////////////////////////////////////////////////////////////////////
strCov fAsyparms(const Eigen::VectorXd& theta,     // reduced-form parameters
                 const Eigen::VectorXd& eta,       // residuals
                 const Eigen::MatrixXd& sV,        // Scale V (for non-isolated)
                 const Eigen::VectorXd& sVphi,     // sV * phi
                 const Eigen::MatrixXd& Z,         // instrument matrix (n x k)
                 const Eigen::VectorXd& y,         // dependent variable (n x 1)
                 const Eigen::MatrixXd& endo,      // Endogenous variables
                 const Eigen::MatrixXd& X,         // covariates for niso
                 const Eigen::MatrixXd& W,         // weighting matrix
                 const Eigen::ArrayXi& Iso,        // Indice for isolated
                 const Eigen::ArrayXi& nIso,       // Indice for nonisolated
                 const Eigen::ArrayXi& cumsn,     // cumulative group indices 
                 const Eigen::ArrayXi& nc_gamma,   // index of not common columns in Xiso 
                 const int& dfiso,                 // degree of freedom for isolated
                 const int& dfniso,                // degree of freedom for nonisolated
                 const int& HAC,                   // HAC type
                 const int& S) {                   // Number of subnets
  
  int Kx(X.cols()), Kendo(endo.cols()), Knc(nc_gamma.size()),
  Ktheta(1 + Kendo + Kx + Knc), Kz(Z.cols());
  
  // 0. unscaled residuals
  Eigen::VectorXd uneta(eta);
  uneta(nIso) *= (1 + theta(0));
  
  // 1. Variance of S^0.5 * moment
  Eigen::MatrixXd Vm(Eigen::MatrixXd::Zero(Kz, Kz));
  double serr(std::numeric_limits<double>::quiet_NaN());
  double serriso(std::numeric_limits<double>::quiet_NaN());
  double serrniso(std::numeric_limits<double>::quiet_NaN());
  if (HAC == 0) { //1.1 iid
    
    serr = sqrt(uneta.dot(uneta) / (dfiso + dfniso));
    Eigen::MatrixXd Ze = Z * serr;
    Ze(nIso, Eigen::all) /= (1 + theta(0));
    Vm = Ze.transpose() * Ze / S;
    
  } else if (HAC == 1) {// 1.2 iid separately
    
    serriso  = sqrt(eta(Iso).dot(eta(Iso)) / dfiso);
    serrniso = sqrt(eta(nIso).dot(eta(nIso)) / dfniso);
    Eigen::MatrixXd Ze = Z;
    Ze(Iso, Eigen::all)  *= serriso;
    Ze(nIso, Eigen::all) *= serrniso;
    Vm = Ze.transpose() * Ze / S;
    serrniso *= (1 + theta(0));
    
  } else if (HAC == 2) { // 1.3 heteroskedasticity
    
    Eigen::MatrixXd Ze    = Z.array().colwise() * eta.array();
    Vm = Ze.transpose() * Ze / S;
    
  } else {
    
    Eigen::ArrayXXd Ze    = Z.array().colwise() * eta.array();
    for (int s(0); s < S; ++ s) { // 1.4 clustering
      int n1(cumsn(s)), n2(cumsn(s + 1) - 1); 
      Eigen::RowVectorXd Zes(Ze(Eigen::seq(n1, n2), Eigen::all).colwise().sum());
      Vm += Zes.transpose() * Zes;
    }
    Vm /= S;
    
  }
  
  // 2. Jacobian J = d Z'eta / d theta, where eta= y-Vphi/(1+betal) for niso and eta = y - Vphi for iso
  Eigen::MatrixXd  J(Kz, Ktheta);
  J << Z(nIso, Eigen::all).transpose() * sVphi(nIso) / (S * (1 + theta(0))), 
       -Z.transpose() * sV / S;
  
  // 3 Estimator's variance: Var = (J'WJ)^-1 J'W Vm W'J (J'WJ)'^-1,
  Eigen::MatrixXd JtW(J.transpose() * W); 
  Eigen::MatrixXd bread_inv((JtW * J).inverse());
  Eigen::MatrixXd meat(JtW * Vm * JtW.transpose());
  Eigen::MatrixXd Vred(bread_inv * meat * bread_inv.transpose() / S);
  
  // 4 Variance matrice for the structural parameters (DELTA METHOD)
  // Vstr = D Vred D'
  Eigen::MatrixXd D = Eigen::MatrixXd::Identity(Ktheta, Ktheta);
  if (Kendo == 2) {// conformity is asymmetric
    // betal = theta(0)
    // theta(1) = delta + betal, so delta = theta(1) - theta(0)
    // betah - betal = theta(2), so betah = theta(0) + theta(2)
    
    // dbetal / dtheta(0) = 1, dbetal / dtheta(1) = 0, dbetal / dtheta(2) = 0
    // dbetah / dtheta(0) = 1, dbetah / dtheta(1) = 0, dbetah / dtheta(2) = 1
    // ddelta / dtheta(0) = -1, ddelta / dtheta(1) = 1, ddelta / dtheta(2) = 0
    
    D(1,0) = 1.0;
    D(1,1) = 0;
    D(1,2) = 1.0;
    D(2,0) = -1.0;
    D(2,1) = 1.0;
    D(2,2) = 0;
  } 
  
  Eigen::MatrixXd Vstr(D * Vred * D.transpose());  
  
  // structural parameters [betal, betah, delta, gamma] or [beta, delta, gamma]
  Eigen::VectorXd thetastr = D * theta; 
  
  // 5 J-stat (Sargan overidentification test)
  Eigen::VectorXd Ze(Z.transpose() * eta.matrix());
  double stat(Ze.dot(Vm.colPivHouseholderQr().solve(Ze)) / S); // because Vm is the variance of Vm * S^0.5
  double df(Kz - Ktheta);
  // 6 test for asymmetry
  Eigen::ArrayXd TestAsym(2);
  if (Kendo == 2) {
    TestAsym(0) = thetastr(1) - thetastr(0); // betah - betal
    TestAsym(1) = Vstr(0, 0) + Vstr(1, 1) - 2 * Vstr(0, 1);
  } 
  // 7. Output. Put element in the structure Str
  strCov out;
  out.estimate  = thetastr;
  out.cov       = Vstr;
  out.serr      = serr;
  out.serriso   = serriso;
  out.serrniso  = serrniso;
  out.Jstat     = stat;
  out.Jdf       = df;
  out.TestAsym  = TestAsym;     
  return out;
}

////////////////////////////////////////////////////////////////////////
//////////////// Computes the structural parameters ////////////////////
/////////////////////// and Asymmetric Variance ////////////////////////
///////////////// for the model without spillover //////////////////////
////////////////////////////////////////////////////////////////////////
strCov fAsyparms_nospill(const Eigen::VectorXd& theta,     // reduced-form parameters
                         const Eigen::VectorXd& eta,       // residuals
                         const Eigen::MatrixXd& sV,        // Scale V (for non-isolated)
                         const Eigen::VectorXd& sVphi,     // sV * phi
                         const Eigen::MatrixXd& Z,         // instrument matrix (n x k)
                         const Eigen::VectorXd& y,         // dependent variable (n x 1)
                         const Eigen::MatrixXd& endo,      // Endogenous variables
                         const Eigen::MatrixXd& X,         // covariates for niso
                         const Eigen::MatrixXd& W,         // weighting matrix
                         const Eigen::ArrayXi& Iso,        // Indice for isolated
                         const Eigen::ArrayXi& nIso,       // Indice for nonisolated
                         const Eigen::ArrayXi& cumsn,     // cumulative group indices 
                         const Eigen::ArrayXi& nc_gamma,   // index of not common columns in Xiso 
                         const int& dfiso,                 // degree of freedom for isolated
                         const int& dfniso,                // degree of freedom for nonisolated
                         const int& HAC,                   // HAC type
                         const int& S) {                   // Number of subnets
  
  int Kx(X.cols()), Kendo(endo.cols()), Knc(nc_gamma.size()), Ktheta(Kendo + Kx + Knc), 
  Kz(Z.cols());
  
  // 0. unscaled residuals
  Eigen::VectorXd uneta(eta);
  uneta(nIso) *= (1 + theta(0));
  
  // 1. Variance of S^0.5 * moment
  Eigen::MatrixXd Vm(Eigen::MatrixXd::Zero(Kz, Kz));
  double serr(std::numeric_limits<double>::quiet_NaN());
  double serriso(std::numeric_limits<double>::quiet_NaN());
  double serrniso(std::numeric_limits<double>::quiet_NaN());
  if (HAC == 0) { //1.1 iid
    
    serr = sqrt(uneta.dot(uneta) / (dfiso + dfniso));
    Eigen::MatrixXd Ze = Z * serr;
    Ze(nIso, Eigen::all) /= (1 + theta(0));
    Vm = Ze.transpose() * Ze / S;
    
  } else if (HAC == 1) {// 1.2 iid separately
    
    serriso  = sqrt(eta(Iso).dot(eta(Iso)) / dfiso);
    serrniso = sqrt(eta(nIso).dot(eta(nIso)) / dfniso);
    Eigen::MatrixXd Ze = Z;
    Ze(Iso, Eigen::all)  *= serriso;
    Ze(nIso, Eigen::all) *= serrniso;
    Vm = Ze.transpose() * Ze / S;
    serrniso *= (1 + theta(0));
    
  } else if (HAC == 2) { // 1.3 heteroskedasticity
    
    Eigen::MatrixXd Ze    = Z.array().colwise() * eta.array();
    Vm = Ze.transpose() * Ze / S;
    
  } else {
    
    Eigen::ArrayXXd Ze    = Z.array().colwise() * eta.array();
    for (int s(0); s < S; ++ s) { // 1.4 clustering
      int n1(cumsn(s)), n2(cumsn(s + 1) - 1); 
      Eigen::RowVectorXd Zes(Ze(Eigen::seq(n1, n2), Eigen::all).colwise().sum());
      Vm += Zes.transpose() * Zes;
    }
    Vm /= S;
    
  }
  
  // 2. Jacobian J = d Z'eta / d theta, where eta= y-Vphi/(1+betal) for niso and eta = y - Vphi for iso
  Eigen::MatrixXd  J(Kz, Ktheta);
  J << Z(nIso, Eigen::all).transpose() * (sVphi(nIso) - endo(nIso, 0) / (1 + theta(0))) / (S * (1 + theta(0))), 
       -Z.transpose() * sV / S;
  
  // 3 Estimator's variance: Var = (J'WJ)^-1 J'W Vm W'J (J'WJ)'^-1,
  Eigen::MatrixXd JtW(J.transpose() * W); 
  Eigen::MatrixXd bread_inv((JtW * J).inverse());
  Eigen::MatrixXd meat(JtW * Vm * JtW.transpose());
  Eigen::MatrixXd Vred(bread_inv * meat * bread_inv.transpose() / S);
  
  // 4. Covariance matrice for the structural parameters (DELTA METHOD)
  // Vstr = D Vred D'
  Eigen::MatrixXd D = Eigen::MatrixXd::Identity(Ktheta, Ktheta);
  if (Kendo == 2) {
    // betal = theta(0)
    // betah - betal = theta(1), so betah = theta(0) + theta(1)
    
    // dbetal / dtheta(0) = 1, dbetal / dtheta(1) = 0
    // dbetah / dtheta(0) = 1, dbetah / dtheta(1) = 1
    
    D(1, 0) = 1.0;
  } 
  
  Eigen::MatrixXd Vstr(D * Vred * D.transpose());  
  
  // Structural parameters [betal, betah, gamma] or [beta, gamma]
  Eigen::VectorXd thetastr = D * theta; 
  
  // 5. J-stat (Sargan overidentification test)
  Eigen::VectorXd Ze(Z.transpose() * eta.matrix());
  double stat(Ze.dot(Vm.colPivHouseholderQr().solve(Ze)) / S); // because Vm is the variance of Vm * S^0.5
  double df(Kz - Ktheta);
  
  // 6. test for asymmetry
  Eigen::ArrayXd TestAsym(2);
  if (Kendo == 2) {
    TestAsym(0) = thetastr(1) - thetastr(0); // betah - betal
    TestAsym(1) = Vstr(0, 0) + Vstr(1, 1) - 2 * Vstr(0, 1);
  } 
  
  // 7. Output. Put element in the structure Str
  strCov out;
  out.estimate  = thetastr;
  out.cov       = Vstr;
  out.serr      = serr;
  out.serriso   = serriso;
  out.serrniso  = serrniso;
  out.Jstat     = stat;
  out.Jdf       = df;
  out.TestAsym  = TestAsym;     
  return out;
}

////////////////////////////////////////////////////////////////////////
//////////////////// Main function to call from R //////////////////////
///////////// It computes the estimates and the cov matrice ////////////
////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::List fAsyMain(const double betal0,              // starting value for betal
                    const Eigen::MatrixXd& Z,         // instrument matrix (n x k)
                    const Eigen::VectorXd& y,         // dependent variable (n x 1)
                    const Eigen::MatrixXd& endo,      // Endogenous variables
                    const Eigen::MatrixXd& X,         // covariates for niso
                    Eigen::MatrixXd& W,               // weighting matrix
                    const Eigen::ArrayXi& Iso,        // Indice for isolated
                    const Eigen::ArrayXi& nIso,       // Indice for nonisolated
                    const Eigen::ArrayXi& cumsn,     // cumulative group indices 
                    const Eigen::ArrayXi& nc_gamma,   // index of not common columns in Xiso 
                    const int& dfiso,                 // degree of freedom for isolated
                    const int& dfniso,                // degree of freedom for nonisolated
                    const int& HAC,                   // HAC type
                    const int& weight,                // weight type
                    const int& S,                     // Number of subnets
                    const int& maxit,                 // optimizer controls
                    const double& eps_f,
                    const double& eps_g,
                    const bool& spillover) {
  int n(y.size()), Kx(X.cols()), Kendo(endo.cols()), Knc(nc_gamma.size()), 
  Ktheta(spillover + Kendo + Kx + Knc);

  // 1. V
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, Ktheta - 1);
  if (spillover){
    V(Eigen::all, Eigen::seqN(0, Kendo + Kx)) << endo, X;
  } else {
    if (Kendo == 2) {
      V(nIso, 0) = endo(nIso, 1);
    }
    V(Eigen::all, Eigen::seqN(Kendo - 1, Kx)) = X;
  }
  if (Knc > 0) {
    V(Iso, Kendo - (1 - spillover) + nc_gamma).setConstant(0); // These columns will have different coeff for iso
    V(Iso, Eigen::seqN(Kendo - (1 - spillover) + Kx, Knc)) = X(Iso, nc_gamma); // Add them as additional columns for iso
  }
  
  // First step GMM
  strEstim gmm;
  if (spillover) {
    gmm = fAsyGmmEstim(betal0, Z, y, V, W, nIso, n, S, Ktheta, maxit, eps_f, eps_g);
  } else {
    gmm = fAsyGmmEstim_nospill(betal0, Z, y, endo.col(0), V, W, nIso, n, S, Ktheta, 
                               maxit, eps_f, eps_g);
  }
  
  // Optimal weighting matrix if necessary
  if(weight == 2) { 
    W = fAsyWopt(gmm.theta, Z, gmm.eta, cumsn, Iso, nIso, dfiso, dfniso, HAC, S);
    // Optimal GMM
    if (spillover) {
      gmm = fAsyGmmEstim(betal0, Z, y, V, W, nIso, n, S, Ktheta, maxit, eps_f, eps_g);
    } else {
      gmm = fAsyGmmEstim_nospill(betal0, Z, y, endo.col(0), V, W, nIso, n, S, Ktheta, 
                                 maxit, eps_f, eps_g);
    }
  }
  
  // Covariance
  strCov final;
  if (spillover) {
    final = fAsyparms(gmm.theta, gmm.eta, gmm.sV, gmm.sVphi, Z, y, endo, X, W, Iso, 
                      nIso, cumsn, nc_gamma, dfiso, dfniso, HAC, S);
  } else {
    final = fAsyparms_nospill(gmm.theta, gmm.eta, gmm.sV, gmm.sVphi, Z, y, endo, X, W, Iso, 
                              nIso, cumsn, nc_gamma, dfiso, dfniso, HAC, S);
  }
  
  // unscale residual
  gmm.eta(nIso) *= (1 + gmm.theta(0));
  
  return Rcpp::List::create(Rcpp::_["estimate"]      = final.estimate,
                            Rcpp::_["cov"]           = final.cov,
                            Rcpp::_["serr"]          = final.serr,
                            Rcpp::_["serr_iso"]      = final.serriso,
                            Rcpp::_["serr_niso"]     = final.serrniso,
                            Rcpp::_["Jstat"]         = final.Jstat,
                            Rcpp::_["Jdf"]           = final.Jdf,
                            Rcpp::_["TestAsym"]      = final.TestAsym,
                            Rcpp::_["objective"]     = gmm.objective,  
                            Rcpp::_["gradient"]      = gmm.gradient,  
                            Rcpp::_["status"]        = gmm.status,
                            Rcpp::_["unscale.resid"] = gmm.eta);
}
