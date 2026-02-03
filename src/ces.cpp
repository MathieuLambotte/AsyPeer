//////////////////////////////////
/// ACKNOWLEDGMENT
/// Taken from QuantilePeer and modified according to the MIT License; QuantilePeer is cited.
/// Houndetoungan A (2025). QuantilePeer: Quantile Peer Effect Models. 
/// doi:10.32614/CRAN.package.QuantilePeer, R package version 0.0.1,
//////////////////////////////////

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

using namespace Rcpp;
using namespace Eigen;
using namespace std;

///////////////////////////////////////////////////////////////////////
////////////////////////////// Structures /////////////////////////////
///////////////////////////////////////////////////////////////////////
struct strListEstim { // output for many estimations
  std::vector<Eigen::VectorXd> theta;
  Eigen::ArrayXd objective;
  std::vector<Eigen::VectorXd> gradient;
  Eigen::ArrayXi status;
}; 

struct strEstim { // output for a single estimation
  Eigen::VectorXd theta;
  double objective;
  Eigen::VectorXd gradient;
  int status;
}; 

struct strCov { // parms and cov with 2 inst
  Eigen::VectorXd estimate;
  Eigen::MatrixXd cov;
  double serr;
  double serriso;
  double serrniso; 
  Eigen::VectorXd residuals;
  double testrho;
}; 


///////////////////////////////////////////////////
///////////////// Data Preparation ////////////////
///////////////////////////////////////////////////
// [[Rcpp::export]]
Eigen::ArrayXXd fCESdata(const Eigen::ArrayXXd& X, 
                         const Eigen::ArrayXd& y, 
                         const Eigen::ArrayXd& z, 
                         const std::vector<Eigen::ArrayXXd>& G,
                         const std::vector<Eigen::ArrayXi>& friendindex,
                         const Eigen::ArrayXi& cumsn,
                         const Eigen::ArrayXi& frzeroy, 
                         const Eigen::ArrayXi& frzeroz, 
                         const std::vector<Eigen::ArrayXi>& lIso, //in selection 
                         const std::vector<Eigen::ArrayXi>& lnIso,//in selection 
                         const Eigen::ArrayXi& nvec, 
                         const Eigen::ArrayXXd& yFMiMa,
                         const Eigen::ArrayXXd& zFMiMa,
                         const int& n,
                         const int& Kx, 
                         const int& S,
                         const double& rho, 
                         const bool& FE,
                         const bool& deriv,
                         const unsigned int& nthread) {
  Eigen::ArrayXd Gy(n), dGy(n), Gz(n), dGz(n), ddGz(n);
  Gy.setConstant(0);
  dGy.setConstant(0);
  Gz.setConstant(0);
  dGz.setConstant(0);
  ddGz.setConstant(0);
  
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int s = 0; s < S; s++) {
    // Extract data for group s
    int nm(nvec(s)), n1(cumsn(s));
    Eigen::ArrayXXd Xm(X(Eigen::seqN(n1, nm), Eigen::all));
    Eigen::ArrayXd ym(y.segment(n1, nm)), zm(z.segment(n1, nm));
    
    for (int i = 0; i < nm; i++) {
      Eigen::ArrayXi fri = friendindex[n1 + i];
      if (fri.size() > 0) {
        if (rho == R_PosInf) {
          Gy[n1 + i] = ym(fri).maxCoeff();
          Gz[n1 + i] = zm(fri).maxCoeff(); 
        } else if (rho == R_NegInf) {
          Gy[n1 + i] = ym(fri).minCoeff();
          Gz[n1 + i] = zm(fri).minCoeff(); 
        } else if (rho < 0) {
          Eigen::ArrayXd Gri(G[s].row(i).transpose());
          
          if (frzeroy[n1 + i] == 0) {
            Eigen::ArrayXd ymf(ym(fri) / yFMiMa(n1 + i, 0)); // y of friends but scaled
            Eigen::ArrayXd tpy(Gri(fri) * pow(ymf, rho));
            double stpy(tpy.sum());
            Gy(n1 + i)    = yFMiMa(n1 + i, 0) * pow(stpy, 1.0 / rho);
            if (deriv) {
              dGy(n1 + i) = Gy(n1 + i) * ((tpy * ymf.log()).sum() / (rho * stpy) - log(stpy) / pow(rho, 2));// Derivative of Gy for gradient computation
            }
          }
          
          if (frzeroz[n1 + i] == 0) {
            Eigen::ArrayXd zmf(zm(fri) / zFMiMa(n1 + i, 0)); // z of friends but scaled
            Eigen::ArrayXd tpz(Gri(fri) * pow(zmf, rho));
            double stpz(tpz.sum());
            Gz(n1 + i)     = zFMiMa(n1 + i, 0) * pow(stpz, 1.0 / rho);
            double tp((tpz * zmf.log()).sum());
            dGz(n1 + i)    = Gz(n1 + i) * (tp / (rho * stpz) - log(stpz) / pow(rho, 2));
            if (deriv) {
              ddGz(n1 + i) = Gz(n1 + i) * (2 * log(stpz) / pow(rho, 3) - 2 * tp/(stpz * pow(rho, 2)) + 
                ((tpz * pow((zmf).log(), 2)).sum() * stpz - pow(tp, 2)) / (rho * pow(stpz, 2)) + 
                pow(tp / (rho * stpz) - log(stpz) / pow(rho, 2), 2));
            }
          }
        } else {
          Eigen::ArrayXd Gri(G[s].row(i).transpose());
          Eigen::ArrayXd ymf(ym(fri) / yFMiMa(n1 + i, 1)); // y of friends but scaled
          Eigen::ArrayXd tpy(Gri(fri) * pow(ymf, rho));
          double stpy(tpy.sum());
          Gy(n1 + i)    = yFMiMa(n1 + i, 1) * pow(stpy, 1.0 / rho);
          if (deriv) {
            dGy(n1 + i) = Gy(n1 + i) * ((tpy * ymf.log()).sum() / (rho * stpy) - log(stpy) / pow(rho, 2));// Derivative of Gy for gradient computation
          }
          
          Eigen::ArrayXd zmf(zm(fri) / zFMiMa(n1 + i, 1)); // z of friends but scaled
          Eigen::ArrayXd tpz(Gri(fri) * pow(zmf, rho));
          double stpz((tpz).sum());
          Gz(n1 + i)    = zFMiMa(n1 + i, 1) * pow(stpz, 1.0 / rho);
          if (frzeroz[n1 + i] == 0) {
            double tp((tpz * zmf.log()).sum());
            dGz(n1 + i) = Gz(n1 + i) * (tp/(rho * stpz) - log(stpz)/pow(rho, 2));
            if (deriv) {
              ddGz(n1 + i) = Gz(n1 + i) * (2 * log(stpz) / pow(rho, 3) - 2 * tp/(stpz * pow(rho, 2)) + 
                ((tpz * pow(zmf.log(), 2)).sum() * stpz - pow(tp, 2)) / (rho * pow(stpz, 2)) + 
                pow(tp / (rho * stpz) - log(stpz) / pow(rho, 2), 2));
            }
          }
        }
      }
    }
  }
#else
  for (int s = 0; s < S; s++) {
    // Extract data for group s
    int nm(nvec(s)), n1(cumsn(s));
    Eigen::ArrayXXd Xm(X(Eigen::seqN(n1, nm), Eigen::all));
    Eigen::ArrayXd ym(y.segment(n1, nm)), zm(z.segment(n1, nm));
    
    for (int i = 0; i < nm; i++) {
      Eigen::ArrayXi fri = friendindex[n1 + i];
      if (fri.size() > 0) {
        if (rho == R_PosInf) {
          Gy[n1 + i] = ym(fri).maxCoeff();
          Gz[n1 + i] = zm(fri).maxCoeff(); 
        } else if (rho == R_NegInf) {
          Gy[n1 + i] = ym(fri).minCoeff();
          Gz[n1 + i] = zm(fri).minCoeff(); 
        } else if (rho < 0) {
          Eigen::ArrayXd Gri(G[s].row(i).transpose());
          
          if (frzeroy[n1 + i] == 0) {
            Eigen::ArrayXd ymf(ym(fri) / yFMiMa(n1 + i, 0)); // y of friends but scaled
            Eigen::ArrayXd tpy(Gri(fri) * pow(ymf, rho));
            double stpy(tpy.sum());
            Gy(n1 + i)    = yFMiMa(n1 + i, 0) * pow(stpy, 1.0 / rho);
            if (deriv) {
              dGy(n1 + i) = Gy(n1 + i) * ((tpy * ymf.log()).sum() / (rho * stpy) - log(stpy) / pow(rho, 2));// Derivative of Gy for gradient computation
            }
          }
          
          if (frzeroz[n1 + i] == 0) {
            Eigen::ArrayXd zmf(zm(fri) / zFMiMa(n1 + i, 0)); // z of friends but scaled
            Eigen::ArrayXd tpz(Gri(fri) * pow(zmf, rho));
            double stpz(tpz.sum());
            Gz(n1 + i)     = zFMiMa(n1 + i, 0) * pow(stpz, 1.0 / rho);
            double tp((tpz * zmf.log()).sum());
            dGz(n1 + i)    = Gz(n1 + i) * (tp / (rho * stpz) - log(stpz) / pow(rho, 2));
            if (deriv) {
              ddGz(n1 + i) = Gz(n1 + i) * (2 * log(stpz) / pow(rho, 3) - 2 * tp/(stpz * pow(rho, 2)) + 
                ((tpz * pow((zmf).log(), 2)).sum() * stpz - pow(tp, 2)) / (rho * pow(stpz, 2)) + 
                pow(tp / (rho * stpz) - log(stpz) / pow(rho, 2), 2));
            }
          }
        } else {
          Eigen::ArrayXd Gri(G[s].row(i).transpose());
          Eigen::ArrayXd ymf(ym(fri) / yFMiMa(n1 + i, 1)); // y of friends but scaled
          Eigen::ArrayXd tpy(Gri(fri) * pow(ymf, rho));
          double stpy(tpy.sum());
          Gy(n1 + i)    = yFMiMa(n1 + i, 1) * pow(stpy, 1.0 / rho);
          if (deriv) {
            dGy(n1 + i) = Gy(n1 + i) * ((tpy * ymf.log()).sum() / (rho * stpy) - log(stpy) / pow(rho, 2));// Derivative of Gy for gradient computation
          }
          
          Eigen::ArrayXd zmf(zm(fri) / zFMiMa(n1 + i, 1)); // z of friends but scaled
          Eigen::ArrayXd tpz(Gri(fri) * pow(zmf, rho));
          double stpz((tpz).sum());
          Gz(n1 + i)    = zFMiMa(n1 + i, 1) * pow(stpz, 1.0 / rho);
          if (frzeroz[n1 + i] == 0) {
            double tp((tpz * zmf.log()).sum());
            dGz(n1 + i) = Gz(n1 + i) * (tp/(rho * stpz) - log(stpz)/pow(rho, 2));
            if (deriv) {
              ddGz(n1 + i) = Gz(n1 + i) * (2 * log(stpz) / pow(rho, 3) - 2 * tp/(stpz * pow(rho, 2)) + 
                ((tpz * pow(zmf.log(), 2)).sum() * stpz - pow(tp, 2)) / (rho * pow(stpz, 2)) + 
                pow(tp / (rho * stpz) - log(stpz) / pow(rho, 2), 2));
            }
          }
        }
      }
    }
  }
#endif
  
  Eigen::ArrayXXd data(n, Kx + 6);
  data.block(0, 0, n, Kx) = X;
  data.col(Kx)            = y;
  data.col(Kx + 1)        = Gy;
  data.col(Kx + 2)        = Gz;
  data.col(Kx + 3)        = dGz;
  data.col(Kx + 4)        = dGy;
  data.col(Kx + 5)        = ddGz;
  if (FE) {
#ifdef _OPENMP
    omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
    for (int s = 0; s < S; ++ s) {
      // For isolated
      if (lIso[s].size() > 0) {
        data(lIso[s], Eigen::all).rowwise() -= data(lIso[s], Eigen::all).colwise().mean();
      }
      // For non-isolated
      if (lnIso[s].size() > 0) {
        data(lnIso[s], Eigen::all).rowwise() -= data(lnIso[s], Eigen::all).colwise().mean();
      }
    }
#else
    for (int s = 0; s < S; ++ s) {
      // For isolated
      if (lIso[s].size() > 0) {
        data(lIso[s], Eigen::all).rowwise() -= data(lIso[s], Eigen::all).colwise().mean();
      }
      // For non-isolated
      if (lnIso[s].size() > 0) {
        data(lnIso[s], Eigen::all).rowwise() -= data(lnIso[s], Eigen::all).colwise().mean();
      }
    }
#endif
  }
  return data;
}


/////////////////////////////////////////////////////////////////////////
///////////// A class for the estimation with a fixed rho  //////////////
/////////////////////////////////////////////////////////////////////////
class ClassCESGMM_rho: public Numer::MFuncGrad {
private:
  const double& rho;
  const Eigen::VectorXd& y;
  const Eigen::VectorXd& Gy;
  const Eigen::MatrixXd& X;
  const Eigen::MatrixXd& Z;
  const Eigen::MatrixXd& W;
  const Eigen::ArrayXi& nIso; 
  const Eigen::ArrayXi& sel;
  const int& n;               // Sample size
  const int& Ktheta;          // Number of parameters 
public:
  ClassCESGMM_rho(const double& rho_,
                  const Eigen::VectorXd& y_,
                  const Eigen::VectorXd& Gy_,
                  const Eigen::MatrixXd& X_,
                  const Eigen::MatrixXd& Z_,
                  const Eigen::MatrixXd& W_,
                  const Eigen::ArrayXi& nIso_,  
                  const Eigen::ArrayXi& sel_,
                  const int& n_,
                  const int& Ktheta_) :
  rho(rho_),
  y(y_),              
  Gy(Gy_),             
  X(X_),
  Z(Z_),
  W(W_),
  nIso(nIso_), 
  sel(sel_),
  n(n_),
  Ktheta(Ktheta_){}
  
  Eigen::VectorXd theta; // the full reduced-form parameter
  Eigen::VectorXd Grad;  // Gradient to be exported
  
  // function computing the objective function and gradient
  double f_grad(Numer::Constvec& tbeta, Numer::Refvec grad) {
    // 0. beta and theta0 (coefficient of Gy)
    double beta   = tbeta(0);
    double theta0 = beta / (1.0 + beta);
    if (beta <= -0.5) {
      grad << 1e300;
      Grad = grad; 
      return 1e300;
    }
    
    // 1. Scale X for non-isolated
    Eigen::MatrixXd sX = X;
    sX(nIso, Eigen::all) /= (1 + beta);
    
    // 2. Closed-form GMM parameters: 
    Eigen::VectorXd u = y - theta0 * Gy;
    Eigen::MatrixXd ZtsX(Z(sel, Eigen::all).transpose() * sX(sel, Eigen::all)); // kz x ksX: Z'sX
    Eigen::VectorXd Ztu(Z(sel, Eigen::all).transpose() * u(sel));  // kz x 1: Z'u
    Eigen::MatrixXd sXtZW(ZtsX.transpose() * W);   //  ksX x kz: sX'Z W
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Adec(sXtZW * ZtsX); // cholesky decomposition for A = sX'Z W Z' sX
    Eigen::VectorXd b(sXtZW * Ztu); // ksX x 1: sX'Z W Z'u
    Eigen::VectorXd phi(Adec.solve(b));
    // Full parameter
    theta.resize(Ktheta);
    theta << rho, beta, phi; 
    
    // 3. residual (but scaled for nonisolated)
    Eigen::VectorXd eta = u - sX * phi;
    
    // 4. Objective function
    Eigen::VectorXd mom(Z(sel, Eigen::all).transpose() * eta(sel)); // Not divised by S to avoid unreached precision
    Eigen::VectorXd Wmom = W * mom;
    double f = mom.dot(Wmom);
    
    // 5. Gradient with respect to beta
    // 5.1 du: derivative of u
    Eigen::VectorXd du(-Gy(sel) / pow(1.0 + beta, 2));
    
    // 5.2 dX: derivative of X
    Eigen::MatrixXd dsX(Eigen::MatrixXd::Zero(n, Ktheta - 2)); 
    dsX(nIso, Eigen::all) = -sX(nIso, Eigen::all) / (1 + beta); // dsX/dbeta
    dsX = dsX(sel, Eigen::all);
    
    // 5.3 dphi: derivative of phi
    Eigen::MatrixXd ZtdsX(Z(sel, Eigen::all).transpose() * dsX); // kz x ksX: Z'dsX
    Eigen::VectorXd ZtdudsX((Z(sel, Eigen::all).transpose() * du) - (ZtdsX * phi)); // Z' (du - dsX * phi)
    Eigen::VectorXd WZteta(W * (Z(sel, Eigen::all).transpose() * eta(sel))); // kz x 1: W Z' eta
    Eigen::VectorXd dphi = Adec.solve(ZtdsX.transpose() * WZteta + sXtZW * ZtdudsX);
    
    grad = 2 * Wmom.transpose() * (ZtdudsX - ZtsX * dphi);
    Grad = grad; // exported
    
    return f;
  }
};

/////////////////////////////////////////////////////////////////////////
//////////////////// Estimation for beta and gamma //////////////////////
////////// given rho  and using two instruments: Gz and dGz /////////////
/////////////////////////////////////////////////////////////////////////
strListEstim fCESestim_rho_2ins(const Eigen::ArrayXd& gridrho,
                                const Eigen::MatrixXd& X,
                                const Eigen::VectorXd& y,
                                const Eigen::VectorXd& z,
                                const std::vector<Eigen::ArrayXXd>& G,
                                const std::vector<Eigen::ArrayXi>& friendindex,
                                const Eigen::ArrayXi& cumsn,
                                const Eigen::ArrayXi& frzeroy,
                                const Eigen::ArrayXi& frzeroz,
                                const std::vector<Eigen::ArrayXi>& lIso, //in selection
                                const std::vector<Eigen::ArrayXi>& lnIso,//in selection
                                const Eigen::ArrayXi& Iso,               //lIso in vector
                                const Eigen::ArrayXi& nIso,              //lnIso in vector
                                const Eigen::ArrayXi& nvec,
                                const Eigen::ArrayXXd& yFMiMa,
                                const Eigen::ArrayXXd& zFMiMa,
                                const Eigen::MatrixXd& W,
                                const Eigen::ArrayXi& idXiso, // Index for Xiso
                                const Eigen::ArrayXi& idXniso,// Index for Xniso
                                const Eigen::ArrayXi& sel,
                                const bool& FE,              
                                const int& maxit,               
                                const double& eps_f,
                                const double& eps_g,
                                const unsigned int& nthread) {
  int Kiso(idXiso.size()), Kniso(idXniso.size()), nrho(gridrho.size()), Kx(X.cols()), 
  Ktheta(2 + Kx), n(X.rows()), S(G.size());
  bool deriv(false);
  
  strListEstim out;
  out.theta     = std::vector<Eigen::VectorXd>(nrho);
  out.objective = Eigen::ArrayXd(nrho);
  out.gradient  = std::vector<Eigen::VectorXd>(nrho);
  out.status    = Eigen::ArrayXi(nrho);
  
#ifdef _OPENMP
  omp_set_num_threads(1); //Not in parallel for now
#pragma omp parallel for schedule(static)
  for (int k = 0; k < nrho; ++ k) {
    // 0. data
    double rho(gridrho(k));
    Eigen::MatrixXd data = fCESdata(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz,
                                    lIso, lnIso, nvec, yFMiMa, zFMiMa, n, Kx, S,
                                    rho, FE, deriv, nthread);
    // this is how data are organized in data
    // X: 0 to Kx - 1
    // y: Kx
    // Gy: Kx + 1
    // Gz: Kx + 2
    // dGz: Kx + 3
    
    // 1. y, Gy(CES of y), X
    Eigen::VectorXd y_(data.col(Kx));
    Eigen::VectorXd Gy_(data.col(Kx + 1));
    Eigen::MatrixXd X_(data(Eigen::all, Eigen::seqN(0, Kx)));
    
    // 2. Instrument
    Eigen::MatrixXd Z(Eigen::MatrixXd::Zero(n, 2 + Kiso + Kniso));
    Z(nIso, 0)           = data(nIso, Kx + 2);
    Z(nIso, 1)           = data(nIso, Kx + 3);
    if (Kiso > 0) { // if there are isolated
      Z(Iso, Eigen::seqN(2, Kiso)) = data(Iso, idXiso);
    }
    Z(nIso, Eigen::seqN(2 + Kiso, Kniso)) = data(nIso, idXniso);
    
    // 3. Optimization
    double obj; // objective
    
    Eigen::VectorXd tbeta(Eigen::VectorXd::Zero(1)); // optimizer tbetal
    ClassCESGMM_rho f(rho, y_, Gy_, X_, Z, W, nIso, sel, n, Ktheta);
    
    out.status(k)    = optim_lbfgs(f, tbeta, obj, maxit, eps_f, eps_g);

    out.theta[k]     = f.theta;
    out.objective(k) = obj / pow(S, 2);
    out.gradient[k]  = f.Grad / pow(S, 2);
  }
#else
  for (int k = 0; k < nrho; ++ k) {
    // 0. data
    double rho(gridrho(k));
    Eigen::MatrixXd data = fCESdata(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz,
                                    lIso, lnIso, nvec, yFMiMa, zFMiMa, n, Kx, S,
                                    rho, FE, deriv, 1);
    // this is how data are organized in data
    // X: 0 to Kx - 1
    // y: Kx
    // Gy: Kx + 1
    // Gz: Kx + 2
    // dGz: Kx + 3
    
    // 1. y, Gy(CES of y), X
    Eigen::VectorXd y_(data.col(Kx));
    Eigen::VectorXd Gy_(data.col(Kx + 1));
    Eigen::MatrixXd X_(data(Eigen::all, Eigen::seqN(0, Kx)));
    
    // 2. Instrument
    Eigen::MatrixXd Z(Eigen::MatrixXd::Zero(n, 2 + Kiso + Kniso));
    Z(nIso, 0)           = data(nIso, Kx + 2);
    Z(nIso, 1)           = data(nIso, Kx + 3);
    if (Kiso > 0) { // if there are isolated
      Z(Iso, Eigen::seqN(2, Kiso)) = data(Iso, idXiso);
    }
    Z(nIso, Eigen::seqN(2 + Kiso, Kniso)) = data(nIso, idXniso);
    
    // 3. Optimization
    double obj; // objective
    Eigen::VectorXd tbeta(Eigen::VectorXd::Zero(1)); // optimizer tbetal
    ClassCESGMM_rho f(rho, y_, Gy_, X_, Z, W, nIso, sel, n, Ktheta);
    out.status(k)    = optim_lbfgs(f, tbeta, obj, maxit, eps_f, eps_g);
    out.theta[k]     = f.theta;
    out.objective(k) = obj / pow(S, 2);
    out.gradient[k]  = f.Grad / pow(S, 2);
  }
#endif
  
  return out;
}

/////////////////////////////////////////////////////////////////////////
////////////// Export to R Estimation for beta and gamma ////////////////
////////// given rho  and using two instruments: Gz and dGz /////////////
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List fCESplotdata(const Eigen::ArrayXd& gridrho,
                   const Eigen::MatrixXd& X,
                   const Eigen::VectorXd& y,
                   const Eigen::VectorXd& z,
                   const std::vector<Eigen::ArrayXXd>& G,
                   const std::vector<Eigen::ArrayXi>& friendindex,
                   const Eigen::ArrayXi& cumsn,
                   const Eigen::ArrayXi& frzeroy,
                   const Eigen::ArrayXi& frzeroz,
                   const std::vector<Eigen::ArrayXi>& lIso, //in selection
                   const std::vector<Eigen::ArrayXi>& lnIso,//in selection
                   const Eigen::ArrayXi& Iso,               //lIso in vector
                   const Eigen::ArrayXi& nIso,              //lnIso in vector
                   const Eigen::ArrayXi& nvec,
                   const Eigen::ArrayXXd& yFMiMa,
                   const Eigen::ArrayXXd& zFMiMa,
                   const Eigen::ArrayXi& idXiso, // Index for Xiso
                   const Eigen::ArrayXi& idXniso,// Index for Xniso
                   const Eigen::ArrayXi& sel,
                   const bool& FE,              
                   const int& maxit,               
                   const double& eps_f,
                   const double& eps_g,
                   const unsigned int& nthread) {
  
  int Kiso(idXiso.size()), Kniso(idXniso.size()), Kz(2 + Kiso + Kniso), 
  Kx(X.cols()), nrho(gridrho.size());
  Eigen::MatrixXd W = Eigen::MatrixXd::Identity(Kz, Kz);
  
  // 1. Estimation for many rho
  strListEstim estim = fCESestim_rho_2ins(gridrho, X, y, z, G, friendindex, cumsn,
                                          frzeroy, frzeroz, lIso, lnIso, Iso, nIso,
                                          nvec, yFMiMa, zFMiMa, W, idXiso, idXniso,
                                          sel, FE, maxit, eps_f, eps_g, nthread);
  
  // 2. Put estimate in a matrix and gradient in a vector
  Eigen::MatrixXd theta(nrho, 2 + Kx);
  Eigen::VectorXd gradient(nrho);
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int k = 0; k < nrho; ++ k) {
    theta.row(k) = estim.theta[k].transpose();
    gradient(k)  = estim.gradient[k](0);
  }
#else
  for (int k = 0; k < nrho; ++ k) {
    theta.row(k) = estim.theta[k].transpose();
    gradient(k)  = estim.gradient[k](0);
  }
#endif

  return Rcpp::List::create(Rcpp::_["estimate"]  = theta,
                            Rcpp::_["objective"] = estim.objective,
                            Rcpp::_["gradient"]  = gradient,
                            Rcpp::_["status"]    = estim.status);
}

/////////////////////////////////////////////////////////////////////////
//////////////////// Estimation for beta and gamma //////////////////////
/////////// given data and using one instrument: Gz and dGz ////////////
/////////////////////////////////////////////////////////////////////////
strEstim fCESestim_rho_1ins(const double& rho,
                            const Eigen::VectorXd& y,
                            const Eigen::VectorXd& Gy,
                            const Eigen::MatrixXd& X,
                            const Eigen::MatrixXd& Z,
                            const Eigen::ArrayXi& nIso,              //lnIso in vector
                            const Eigen::MatrixXd& W,
                            const Eigen::ArrayXi& sel,
                            const int& n,
                            const int& S,
                            const int& Ktheta,
                            const int& maxit,               
                            const double& eps_f,
                            const double& eps_g) {
  // 1. Optimization
  double obj; // objective
  Eigen::VectorXd tbeta(Eigen::VectorXd::Zero(1)); // optimizer tbetal
  ClassCESGMM_rho f(rho, y, Gy, X, Z, W, nIso, sel, n, Ktheta);
  int status = optim_lbfgs(f, tbeta, obj, maxit, eps_f, eps_g);
  
  // output
  strEstim out;
  out.theta     = f.theta;
  out.objective = obj / pow(S, 2);
  out.gradient  = f.Grad / pow(S, 2);
  out.status    = status;
  
  return out;
}


/////////////////////////////////////////////////////////////////////////
//////////////////////////// GMM Weight /////////////////////////////////
///////////////// using two instruments: Gz and dGz /////////////////////
/////////////////////////////////////////////////////////////////////////
Eigen::MatrixXd fCESWeight_2ins(const Eigen::VectorXd& theta,
                                const Eigen::MatrixXd& X,
                                const Eigen::VectorXd& y,
                                const Eigen::VectorXd& z,
                                const std::vector<Eigen::ArrayXXd>& G,
                                const std::vector<Eigen::ArrayXi>& friendindex,
                                const Eigen::ArrayXi& cumsn,
                                const Eigen::ArrayXi& frzeroy,
                                const Eigen::ArrayXi& frzeroz,
                                const std::vector<Eigen::ArrayXi>& lIso, //in selection
                                const std::vector<Eigen::ArrayXi>& lnIso,//in selection
                                const Eigen::ArrayXi& Iso,               //lIso in vector
                                const Eigen::ArrayXi& nIso,              //lnIso in vector
                                const Eigen::ArrayXi& nvec,
                                const Eigen::ArrayXXd& yFMiMa,
                                const Eigen::ArrayXXd& zFMiMa,
                                const Eigen::MatrixXd& W,
                                const Eigen::ArrayXi& idXiso, // Index for Xiso
                                const Eigen::ArrayXi& idXniso,// Index for Xniso
                                const Eigen::ArrayXi& sel,
                                const bool& FE,
                                const int& HACn,
                                const int& dfiso,
                                const int& dfniso,              
                                const unsigned int& nthread) {
  int Kiso(idXiso.size()), Kniso(idXniso.size()), Kx(X.cols()), n(X.rows()), S(G.size());
  bool deriv(false);
  
  // 0. data
  Eigen::MatrixXd data = fCESdata(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz,
                                  lIso, lnIso, nvec, yFMiMa, zFMiMa, n, Kx, S,
                                  theta(0), FE, deriv, nthread);
  // this is how data are organized in data
  // X: 0 to Kx - 1
  // y: Kx
  // Gy: Kx + 1
  // Gz: Kx + 2
  // dGz: Kx + 3
  
  // 1. y, Gy(CES of y), X, and Scale X for non-isolated
  Eigen::VectorXd y_(data.col(Kx));
  Eigen::VectorXd Gy_(data.col(Kx + 1));
  Eigen::MatrixXd X_(data(Eigen::all, Eigen::seqN(0, Kx)));
  Eigen::MatrixXd sX = X_;
  sX(nIso, Eigen::all) /= (1 + theta(1));
  
  // 2. scaled residuals
  Eigen::VectorXd eta = y_ - (theta(1) / (1 + theta(1))) * Gy_ - sX * theta.segment(2, Kx);
  // unscaled residuals
  Eigen::VectorXd uneta = eta;
  uneta(nIso) *= (1 + theta(1));
  
  // 3. Instrument
  Eigen::MatrixXd Z(Eigen::MatrixXd::Zero(n, 2 + Kiso + Kniso));
  Z(nIso, 0)           = data(nIso, Kx + 2);
  Z(nIso, 1)           = data(nIso, Kx + 3);
  if (Kiso > 0) { // if there are isolated
    Z(Iso, Eigen::seqN(2, Kiso)) = data(Iso, idXiso);
  }
  Z(nIso, Eigen::seqN(2 + Kiso, Kniso)) = data(nIso, idXniso);
  
  // 4. Variance of moments
  Eigen::MatrixXd Vm(2 + Kiso + Kniso, 2 + Kiso + Kniso);
  Eigen::MatrixXd Ze = Z;
  if (HACn == 0) {
    
    double serr(sqrt(uneta(sel).dot(uneta(sel)) / (dfiso + dfniso)));
    Ze *= serr;
    Ze(nIso, Eigen::all) /= (1 + theta(1));
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / pow(S, 2);
    
  } else if (HACn == 1) {
    
    double serriso(sqrt(eta(Iso).dot(eta(Iso)) / dfiso));
    double serrniso(sqrt(eta(nIso).dot(eta(nIso)) / dfniso));
    Ze(Iso, Eigen::all)  *= serriso;
    Ze(nIso, Eigen::all) *= serrniso;
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / pow(S, 2);
    
  } else if(HACn == 2) {
    
    Ze.array().colwise() *= eta.array();
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / pow(S, 2);
    
  } else {
    
    for (int s = 0; s < S; ++s) {
      int nisos(lIso[s].size()), nnisos(lnIso[s].size()), ns(nisos + nnisos);
      Ze.array().colwise() *= eta.array();
      if (ns > 0) {
        Eigen::ArrayXi sels(ns);
        sels << lIso[s], lnIso[s];
        Eigen::RowVectorXd Zes(Ze(sels, Eigen::all).array().colwise().sum());
        Vm += Zes.transpose() * Zes;
      }
    }
    Vm /= pow(S, 2);
    
  }
  
  // 5. W
  return Vm.inverse();
}


/////////////////////////////////////////////////////////////////////////
//////////////////////////// GMM Weight /////////////////////////////////
//////////////////// using one instrument: Gz ///////////////////////////
/////////////////////////////////////////////////////////////////////////
Eigen::MatrixXd fCESWeight_1ins(const Eigen::VectorXd& theta,
                                const Eigen::VectorXd& y,
                                const Eigen::VectorXd& Gy,
                                const Eigen::MatrixXd& X,
                                const Eigen::MatrixXd& Z,
                                const Eigen::ArrayXi& Iso,               //lIso in vector
                                const Eigen::ArrayXi& nIso,              //lnIso in vector
                                const std::vector<Eigen::ArrayXi>& lIso, //in selection
                                const std::vector<Eigen::ArrayXi>& lnIso,//in selection
                                const Eigen::ArrayXi& sel,
                                const int& n,
                                const int& S,
                                const int& Kx,
                                const int& Kz,
                                const int& HACn,
                                const int& dfiso,
                                const int& dfniso) {
  // 1. Scale X for non-isolated
  Eigen::MatrixXd sX = X;
  sX(nIso, Eigen::all) /= (1 + theta(1));
  
  // 2. scaled residuals
  Eigen::VectorXd eta = y - (theta(1) / (1 + theta(1))) * Gy - sX * theta.segment(2, Kx);
  // unscaled residuals
  Eigen::VectorXd uneta = eta;
  uneta(nIso) *= (1 + theta(1));
  
  // 4. Variance of moments
  Eigen::MatrixXd Vm(Kz, Kz);
  Eigen::MatrixXd Ze = Z;
  if (HACn == 0) {
    
    double serr(sqrt(uneta(sel).dot(uneta(sel)) / (dfiso + dfniso)));
    Ze *= serr;
    Ze(nIso, Eigen::all) /= (1 + theta(1));
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / pow(S, 2);
    
  } else if (HACn == 1) {
    
    double serriso(sqrt(eta(Iso).dot(eta(Iso)) / dfiso));
    double serrniso(sqrt(eta(nIso).dot(eta(nIso)) / dfniso));
    Ze(Iso, Eigen::all)  *= serriso;
    Ze(nIso, Eigen::all) *= serrniso;
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / pow(S, 2);
    
  } else if(HACn == 2) {
    
    Ze.array().colwise() *= eta.array();
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / pow(S, 2);
    
  } else {
    
    for (int s = 0; s < S; ++s) {
      int nisos(lIso[s].size()), nnisos(lnIso[s].size()), ns(nisos + nnisos);
      Ze.array().colwise() *= eta.array();
      if (ns > 0) {
        Eigen::ArrayXi sels(ns);
        sels << lIso[s], lnIso[s];
        Eigen::RowVectorXd Zes(Ze(sels, Eigen::all).array().colwise().sum());
        Vm += Zes.transpose() * Zes;
      }
    }
    Vm /= pow(S, 2);
    
  }
  
  // 2. W
  return Vm.inverse();
}


/////////////////////////////////////////////////////////////////////////
/////////// A class for the estimation with a variable rho  /////////////
/////////////////////////////////////////////////////////////////////////
class ClassCESGMM: public Numer::MFuncGrad {
private:
  const Eigen::MatrixXd& X;
  const Eigen::VectorXd& y;
  const Eigen::VectorXd& z;
  const std::vector<Eigen::ArrayXXd>& G;
  const std::vector<Eigen::ArrayXi>& friendindex;
  const Eigen::ArrayXi& cumsn;
  const Eigen::ArrayXi& frzeroy;
  const Eigen::ArrayXi& frzeroz;
  const std::vector<Eigen::ArrayXi>& lIso; //in selection
  const std::vector<Eigen::ArrayXi>& lnIso;//in selection
  const Eigen::ArrayXi& Iso;               //lIso in vector
  const Eigen::ArrayXi& nIso;              //lnIso in vector
  const Eigen::ArrayXi& nvec;
  const Eigen::ArrayXXd& yFMiMa;
  const Eigen::ArrayXXd& zFMiMa;
  const Eigen::MatrixXd& W;
  const Eigen::ArrayXi& idXiso; // Index for Xiso
  const Eigen::ArrayXi& idXniso;// Index for Xniso
  const Eigen::ArrayXi& sel;
  const bool& FE;               // Number of networks
  const int& Ktheta;            // Number of parameters 
  const int& n;
  const int& Kx;
  const int& Kz;
  const int& S;
  const int& Kiso;
  const int& Kniso;
  const double& rhomin;
  const double& rhomax;
public:
  ClassCESGMM(const Eigen::MatrixXd& X_,
              const Eigen::VectorXd& y_,
              const Eigen::VectorXd& z_,
              const std::vector<Eigen::ArrayXXd>& G_,
              const std::vector<Eigen::ArrayXi>& friendindex_,
              const Eigen::ArrayXi& cumsn_,
              const Eigen::ArrayXi& frzeroy_,
              const Eigen::ArrayXi& frzeroz_,
              const std::vector<Eigen::ArrayXi>& lIso_, //in selection
              const std::vector<Eigen::ArrayXi>& lnIso_,//in selection
              const Eigen::ArrayXi& Iso_,               //lIso in vector
              const Eigen::ArrayXi& nIso_,              //lnIso in vector
              const Eigen::ArrayXi& nvec_,
              const Eigen::ArrayXXd& yFMiMa_,
              const Eigen::ArrayXXd& zFMiMa_,
              const Eigen::MatrixXd& W_,
              const Eigen::ArrayXi& idXiso_, // Index for Xiso
              const Eigen::ArrayXi& idXniso_,// Index for Xniso
              const Eigen::ArrayXi& sel_,
              const bool& FE_,               // Number of networks
              const int& Ktheta_,
              const int& n_,
              const int& Kx_,
              const int& Kz_,
              const int& S_,
              const int& Kiso_,
              const int& Kniso_,
              const double& rhomin_,
              const double& rhomax_) :
  
  X(X_),
  y(y_),
  z(z_),
  G(G_),
  friendindex(friendindex_),
  cumsn(cumsn_),
  frzeroy(frzeroy_),
  frzeroz(frzeroz_),
  lIso(lIso_),
  lnIso(lnIso_),
  Iso(Iso_),
  nIso(nIso_),
  nvec(nvec_),
  yFMiMa(yFMiMa_),
  zFMiMa(zFMiMa_),
  W(W_),
  idXiso(idXiso_),
  idXniso(idXniso_),
  sel(sel_),
  FE(FE_),
  Ktheta(Ktheta_),
  n(n_),
  Kx(Kx_),
  Kz(Kz_),
  S(S_),
  Kiso(Kiso_),
  Kniso(Kniso_),
  rhomin(rhomin_),
  rhomax(rhomax_){}
  
  Eigen::VectorXd theta; // the full reduced-form parameter
  Eigen::VectorXd Grad;  // Gradient to be exported
  
  // function computing the objective function and gradient
  double f_grad(Numer::Constvec& ttheta01, Numer::Refvec grad) {
    // 0.1 rho and beta 
    double rho(ttheta01(0)), beta(ttheta01(1));
    if ((beta <= -0.5) || (rho <= rhomin) || (rho >= rhomax)) {
      grad << 1e300, 1e300;
      Grad = grad; 
      return 1e300;
    }
    double theta1(beta / (1.0 + beta)); // Coefficient of Gy

    // 0.3 data
    bool deriv(true);
    Eigen::MatrixXd data = fCESdata(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz,
                                    lIso, lnIso, nvec, yFMiMa, zFMiMa, n, Kx, S,
                                    rho, FE, deriv, 1);
    // X: 0 to Kx - 1
    // y: Kx
    // Gy: Kx + 1
    // Gz: Kx + 2
    // dGz: Kx + 3
    // dGy: Kx + 4
    // ddGz: Kx + 5

    // 1. Instruments and Scale X for non-isolated
    // Instruments
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(n, Kz);
    Z(nIso, 0) = data(nIso, Kx + 2); //Gz
    Z(nIso, 1) = data(nIso, Kx + 3); //dGz
    if (Kiso > 0) { // if there are isolated
      Z(Iso, Eigen::seqN(2, Kiso)) = data(Iso, idXiso);
    }
    Z(nIso, Eigen::seqN(2 + Kiso, Kniso)) = data(nIso, idXniso);

    // Scale X for non-isolated
    Eigen::MatrixXd sX = data.block(0, 0, n, Kx);
    sX(nIso, Eigen::all) /= (1 + beta);
    
    // 2. Closed-form GMM parameters: 
    Eigen::VectorXd u = data.col(Kx) - theta1 * data.col(Kx + 1); // u = y - theta1 * Gy
    Eigen::MatrixXd ZtsX(Z(sel, Eigen::all).transpose() * sX(sel, Eigen::all)); // kz x ksX: Z'sX
    Eigen::VectorXd Ztu(Z(sel, Eigen::all).transpose() * u(sel));  // kz x 1: Z'u
    Eigen::MatrixXd sXtZW(ZtsX.transpose() * W);   //  ksX x kz: sX'Z W
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Adec(sXtZW * ZtsX); // cholesky decomposition for A = sX'Z W Z' sX
    Eigen::VectorXd b(sXtZW * Ztu); // ksX x 1: sX'Z W Z'u
    Eigen::VectorXd phi(Adec.solve(b));
    // Full parameter
    theta.resize(Ktheta);
    theta << rho, beta, phi; 

    // 3. residual (but scaled for nonisolated)
    Eigen::VectorXd eta = u - sX * phi;
    
    // 4. Objective function
    Eigen::VectorXd mom(Z(sel, Eigen::all).transpose() * eta(sel)); // Not divised by S to avoid unreached precision
    Eigen::VectorXd Wmom = W * mom;
    double f = mom.dot(Wmom);
    
    // 5. Gradient 
    // 5.1 derivative of u where  u = y - theta1 * Gy
    /// with respect to rho
    Eigen::VectorXd du_r(-theta1 * data(sel, Kx + 4));   // - theta1 * dGy
    /// with respect to beta
    Eigen::VectorXd du_b(-data(sel, Kx + 1) / pow(1.0 + beta, 2)); // - Gy / (1 + beta)^2
    
    // 5.2 dsX: derivative of sX with respect to beta
    Eigen::MatrixXd dsX_b(Eigen::MatrixXd::Zero(n, Kx)); 
    dsX_b(nIso, Eigen::all) = -sX(nIso, Eigen::all) / (1 + beta); // dsX/dbeta
    dsX_b = dsX_b(sel, Eigen::all);
    
    // 5.3  dZ: derivative of Z with respect to rho
    // column 0 of Z is Gz and column 1 is dGz.
    Eigen::MatrixXd dZ_r(Eigen::MatrixXd::Zero(n, Kz)); 
    dZ_r(nIso, 0)        = data(nIso, Kx + 3); //dGz
    dZ_r(nIso, 1)        = data(nIso, Kx + 5); //ddGz
    dZ_r = dZ_r(sel, Eigen::all);
    
    // 5.4 dphi: derivative of phi
    //  dphi_r
    Eigen::VectorXd WZteta(W * Z(sel, Eigen::all).transpose() * eta(sel)); 
    Eigen::VectorXd dZetaZdu_r = dZ_r.transpose() * eta(sel) + Z(sel, Eigen::all).transpose() * du_r;
    Eigen::VectorXd dphi_r = Adec.solve(sX(sel, Eigen::all).transpose() * dZ_r * WZteta + sXtZW * dZetaZdu_r);
    
    // dphi_b
    Eigen::MatrixXd ZtdsX_b(Z(sel, Eigen::all).transpose() * dsX_b); 
    Eigen::VectorXd ZtdudsX_b(Z(sel, Eigen::all).transpose() * du_b - ZtdsX_b * phi); 
    Eigen::VectorXd dphi_b = Adec.solve(ZtdsX_b.transpose() * WZteta + sXtZW * ZtdudsX_b);
    
    // gradient
    grad(0) = 2 * Wmom.dot(dZetaZdu_r - ZtsX * dphi_r);
    grad(1) = 2 * Wmom.dot(ZtdudsX_b - ZtsX * dphi_b);
    Grad  = grad; // exported

    return f;
  }
};

/////////////////////////////////////////////////////////////////////////
///////////////// Estimation for rho beta and gamma /////////////////////
/////////////////////////////////////////////////////////////////////////
strEstim fCESestim(const Eigen::MatrixXd& X,
                   const Eigen::VectorXd& y,
                   const Eigen::VectorXd& z,
                   const std::vector<Eigen::ArrayXXd>& G,
                   const std::vector<Eigen::ArrayXi>& friendindex,
                   const Eigen::ArrayXi& cumsn,
                   const Eigen::ArrayXi& frzeroy,
                   const Eigen::ArrayXi& frzeroz,
                   const std::vector<Eigen::ArrayXi>& lIso, //in selection
                   const std::vector<Eigen::ArrayXi>& lnIso,//in selection
                   const Eigen::ArrayXi& Iso,               //lIso in vector
                   const Eigen::ArrayXi& nIso,              //lnIso in vector
                   const Eigen::ArrayXi& nvec,
                   const Eigen::ArrayXXd& yFMiMa,
                   const Eigen::ArrayXXd& zFMiMa,
                   const Eigen::MatrixXd& W,
                   const Eigen::ArrayXi& idXiso, // Index for Xiso
                   const Eigen::ArrayXi& idXniso,// Index for Xniso
                   const Eigen::ArrayXi& sel,
                   const double& rhomin,
                   const double& rhomax,
                   const bool& FE,      
                   const int& nthread,
                   const int& maxit,               
                   const double& eps_f,
                   const double& eps_g) {
  int Kiso(idXiso.size()), Kniso(idXniso.size()), Kx(X.cols()), Ktheta(2 + Kx), 
  n(X.rows()), S(G.size()), Kz(2 + Kiso + Kniso);
  
  strEstim out;
  
  if (rhomax * rhomin < 0) { // estimation for positive rho and for negative rho
    std::vector<strEstim> lout(2);
    Eigen::VectorXd lrhomin(2), lrhomax(2);
    lrhomin << rhomin, 0;
    lrhomax << 0, rhomax;
    
    for (int k = 0; k < 2; ++k) {
      ClassCESGMM f(X, y, z, G, friendindex, 
                    cumsn, frzeroy, frzeroz, lIso, lnIso, 
                    Iso, nIso, nvec, yFMiMa, zFMiMa, W, idXiso, idXniso, sel, FE, 
                    Ktheta, n, Kx, Kz, S, Kiso, Kniso, lrhomin(k), lrhomax(k));
      
      Eigen::VectorXd ttheta01(Eigen::VectorXd::Zero(2)); // optimizer ttheta01
      ttheta01(0) = 0.5 * (lrhomin(k) + lrhomax(k));
      double obj; 
      lout[k].status = optim_lbfgs(f, ttheta01, obj, maxit, eps_f, eps_g);
      lout[k].theta     = f.theta;
      lout[k].gradient  = f.Grad / pow(S, 2);
      lout[k].objective = obj / pow(S, 2);
    }
    if (lout[0].objective < lout[1].objective) {
      out = lout[0];
    } else{
      out = lout[1];
    }
  } else{
    ClassCESGMM f(X, y, z, G, friendindex, 
                  cumsn, frzeroy, frzeroz, lIso, lnIso, 
                  Iso, nIso, nvec, yFMiMa, zFMiMa, W, idXiso, idXniso, sel, FE, 
                  Ktheta, n, Kx, Kz, S, Kiso, Kniso, rhomin, rhomax);
    
    Eigen::VectorXd ttheta01(Eigen::VectorXd::Zero(2)); 
    ttheta01(0) = 0.5 * (rhomin + rhomax);
    double obj; 
    out.status = optim_lbfgs(f, ttheta01, obj, maxit, eps_f, eps_g);
    out.theta     = f.theta;
    out.gradient  = f.Grad / pow(S, 2);
    out.objective = obj / pow(S, 2);
  }

  return out;
}


/////////////////////////////////////////////////////////////////////////
//////////////////// full parameters and variances //////////////////////
//////////////////////// with a flexible rho ////////////////////////////
/////////////////////////////////////////////////////////////////////////
strCov fCESparms(const Eigen::VectorXd& theta,
                 const Eigen::MatrixXd& X,
                 const Eigen::VectorXd& y,
                 const Eigen::VectorXd& z,
                 const std::vector<Eigen::ArrayXXd>& G,
                 const std::vector<Eigen::ArrayXi>& friendindex,
                 const Eigen::ArrayXi& cumsn,
                 const Eigen::ArrayXi& frzeroy,
                 const Eigen::ArrayXi& frzeroz,
                 const std::vector<Eigen::ArrayXi>& lIso, //in selection
                 const std::vector<Eigen::ArrayXi>& lnIso,//in selection
                 const Eigen::ArrayXi& Iso,               //lIso in vector
                 const Eigen::ArrayXi& nIso,              //lnIso in vector
                 const Eigen::ArrayXi& nvec,
                 const Eigen::ArrayXXd& yFMiMa,
                 const Eigen::ArrayXXd& zFMiMa,
                 const Eigen::MatrixXd& W,
                 const Eigen::ArrayXi& idXiso, // Index for Xiso
                 const Eigen::ArrayXi& idXniso,// Index for Xniso
                 const Eigen::ArrayXi& sel,
                 const bool& FE,
                 const int& HACn,
                 const int& dfiso,
                 const int& dfniso,              
                 const unsigned int& nthread) {
  int Kiso(idXiso.size()), Kniso(idXniso.size()), Kx(X.cols()), n(X.rows()), 
  S(G.size()), Kz(2 + Kiso + Kniso);
  double rho(theta(0)), beta(theta(1));
  bool deriv(true);
  
  // 0. data
  Eigen::MatrixXd data = fCESdata(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz,
                                  lIso, lnIso, nvec, yFMiMa, zFMiMa, n, Kx, S,
                                  rho, FE, deriv, nthread);
  // this is how data are organized in data
  // X: 0 to Kx - 1
  // y: Kx
  // Gy: Kx + 1
  // Gz: Kx + 2
  // dGz: Kx + 3
  // dGy: Kx + 4
  // ddGz: Kx + 5
  
  // 1. y, Gy(CES of y), X, and Scale X for non-isolated
  Eigen::VectorXd y_(data.col(Kx));
  Eigen::VectorXd Gy_(data.col(Kx + 1));
  Eigen::MatrixXd X_(data(Eigen::all, Eigen::seqN(0, Kx)));
  Eigen::MatrixXd sX = X_;
  sX(nIso, Eigen::all) /= (1 + theta(1));
  
  // 2. scaled residuals
  Eigen::VectorXd eta = y_ - (theta(1) / (1 + theta(1))) * Gy_ - sX * theta.segment(2, Kx);
  // unscaled residuals
  Eigen::VectorXd uneta = eta;
  uneta(nIso) *= (1 + theta(1));
  
  // 3. Instrument
  Eigen::MatrixXd Z(Eigen::MatrixXd::Zero(n, Kz));
  Z(nIso, 0)           = data(nIso, Kx + 2);
  Z(nIso, 1)           = data(nIso, Kx + 3);
  if (Kiso > 0) { // if there are isolated
    Z(Iso, Eigen::seqN(2, Kiso)) = data(Iso, idXiso);
  }
  Z(nIso, Eigen::seqN(2 + Kiso, Kniso)) = data(nIso, idXniso);
  
  // 4. Variance of S^0.5 * moment
  Eigen::MatrixXd Vm(Eigen::MatrixXd::Zero(Kz, Kz));
  double serr(std::numeric_limits<double>::quiet_NaN());
  double serriso(std::numeric_limits<double>::quiet_NaN());
  double serrniso(std::numeric_limits<double>::quiet_NaN());
  Eigen::MatrixXd Ze = Z;
  if (HACn == 0) {
    
    serr = sqrt(uneta(sel).dot(uneta(sel)) / (dfiso + dfniso));
    Ze *= serr;
    Ze(nIso, Eigen::all) /= (1 + theta(1));
    Vm  = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / S;
    
  } else if (HACn == 1) {
    
    serriso  = sqrt(eta(Iso).dot(eta(Iso)) / dfiso);
    serrniso = sqrt(eta(nIso).dot(eta(nIso)) / dfniso);
    Ze(Iso, Eigen::all)  *= serriso;
    Ze(nIso, Eigen::all) *= serrniso;
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / S;
    serrniso *= (1 + theta(1));
    
  } else if(HACn == 2) {
    
    Ze.array().colwise() *= eta.array();
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / S;
    
  } else {
    
    for (int s = 0; s < S; ++s) {
      int nisos(lIso[s].size()), nnisos(lnIso[s].size()), ns(nisos + nnisos);
      Ze.array().colwise() *= eta.array();
      if (ns > 0) {
        Eigen::ArrayXi sels(ns);
        sels << lIso[s], lnIso[s];
        Eigen::RowVectorXd Zes(Ze(sels, Eigen::all).array().colwise().sum());
        Vm += Zes.transpose() * Zes;
      }
    }
    Vm /= S;
    
  }
  
  // 5. Jacobian of the moment
  Eigen::MatrixXd J(Kz, 2 + Kx);
  // 5.1 derivative of u where  u = y - (beta / (1 + beta)) * Gy
  /// with respect to rho
  Eigen::VectorXd du_r(-(beta / (1 + beta)) * data(sel, Kx + 4));   
  /// with respect to beta
  Eigen::VectorXd du_b(-data(sel, Kx + 1) / pow(1 + beta, 2)); 
  
  // 5.2 dsX: derivative of sX with respect to beta
  Eigen::MatrixXd dsX_b(Eigen::MatrixXd::Zero(n, Kx)); 
  dsX_b(nIso, Eigen::all) = -sX(nIso, Eigen::all) / (1 + beta); 
  dsX_b = dsX_b(sel, Eigen::all);
  
  // 5.3  dZ: derivative of Z with respect to rho
  // column 0 of Z is Gz and column 1 is dGz.
  Eigen::MatrixXd dZ_r(Eigen::MatrixXd::Zero(n, Kz)); 
  dZ_r(nIso, 0)        = data(nIso, Kx + 3); //dGz
  dZ_r(nIso, 1)        = data(nIso, Kx + 5); //ddGz
  dZ_r = dZ_r(sel, Eigen::all);
  
  // Jacobian
  // J.col(0) = dZ'/drho * eta + Z' * du/drho
  J.col(0) = dZ_r.transpose() * eta(sel) + Z(sel, Eigen::all).transpose() * du_r;
  // J.col(1) = Z' * (- Gy / (1 + beta)^2 - d sX / dbeta * theta.segment(2, Kx))
  J.col(1) = Z(sel, Eigen::all).transpose() * (-Gy_(sel) / pow(1 + beta, 2) - dsX_b * theta.segment(2, Kx));
  // J.col(2 and +) = Z' * sX
  J.block(0, 2, Kz, Kx) = Z(sel, Eigen::all).transpose() * sX(sel, Eigen::all);
  // Normalization 
  J /= S;
  
  // 6 Compute variance: Var = inv(J'WJ) J'WVmW'J inv(J'WJ) / S
  Eigen::MatrixXd JtW(J.transpose() * W);
  Eigen::MatrixXd bread_inv((JtW * J).inverse());
  Eigen::MatrixXd meat(JtW * Vm * JtW.transpose());
  Eigen::MatrixXd Var(bread_inv * meat * bread_inv.transpose() / S);
  
  // 7. Testing whether rho = 1
  double testrho = pow(rho - 1, 2) / Var(0, 0);
  
  // output
  strCov out;
  out.estimate  = theta;
  out.cov       = Var;
  out.serr      = serr;
  out.serriso   = serriso;
  out.serrniso  = serrniso;
  out.residuals = uneta;
  out.testrho   = testrho;
  return out;
}



/////////////////////////////////////////////////////////////////////////
//////////////////// full parameters and variances //////////////////////
/////////////////////////// with a fixed rho ////////////////////////////
/////////////////////////////////////////////////////////////////////////
strCov fCESparms_rho(const Eigen::VectorXd& theta,
                     const Eigen::VectorXd& y,
                     const Eigen::VectorXd& Gy,
                     const Eigen::MatrixXd& X,
                     const Eigen::MatrixXd& Z,
                     const Eigen::ArrayXi& Iso,               //lIso in vector
                     const Eigen::ArrayXi& nIso,              //lnIso in vector
                     const std::vector<Eigen::ArrayXi>& lIso, //in selection
                     const std::vector<Eigen::ArrayXi>& lnIso,//in selection
                     const Eigen::MatrixXd& W,
                     const Eigen::ArrayXi& sel,
                     const int& n,
                     const int& S,
                     const int& Kx,
                     const int& Kz,
                     const int& HACn,
                     const int& dfiso,
                     const int& dfniso) {
  
  // 1. y, Gy(CES of y), X, and Scale X for non-isolated
  Eigen::MatrixXd sX = X;
  sX(nIso, Eigen::all) /= (1 + theta(1));
  
  // 2. scaled residuals
  Eigen::VectorXd eta = y - (theta(1) / (1 + theta(1))) * Gy - sX * theta.segment(2, Kx);
  // unscaled residuals
  Eigen::VectorXd uneta = eta;
  uneta(nIso) *= (1 + theta(1));
  
  // 3. Variance of S^0.5 * moment
  Eigen::MatrixXd Vm(Eigen::MatrixXd::Zero(Kz, Kz));
  double serr(std::numeric_limits<double>::quiet_NaN());
  double serriso(std::numeric_limits<double>::quiet_NaN());
  double serrniso(std::numeric_limits<double>::quiet_NaN());
  Eigen::MatrixXd Ze = Z;
  if (HACn == 0) {
    
    serr = sqrt(uneta(sel).dot(uneta(sel)) / (dfiso + dfniso));
    Ze *= serr;
    Ze(nIso, Eigen::all) /= (1 + theta(1));
    Vm  = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / S;
    
  } else if (HACn == 1) {
    
    serriso  = sqrt(eta(Iso).dot(eta(Iso)) / dfiso);
    serrniso = sqrt(eta(nIso).dot(eta(nIso)) / dfniso);
    Ze(Iso, Eigen::all)  *= serriso;
    Ze(nIso, Eigen::all) *= serrniso;
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / S;
    serrniso *= (1 + theta(1));
    
  } else if(HACn == 2) {
    
    Ze.array().colwise() *= eta.array();
    Vm = Ze(sel, Eigen::all).transpose() * Ze(sel, Eigen::all) / S;
    
  } else {
    
    for (int s = 0; s < S; ++s) {
      int nisos(lIso[s].size()), nnisos(lnIso[s].size()), ns(nisos + nnisos);
      Ze.array().colwise() *= eta.array();
      if (ns > 0) {
        Eigen::ArrayXi sels(ns);
        sels << lIso[s], lnIso[s];
        Eigen::RowVectorXd Zes(Ze(sels, Eigen::all).array().colwise().sum());
        Vm += Zes.transpose() * Zes;
      }
    }
    Vm /= S;
    
  }
  
  // 5. Jacobian of the moment
  Eigen::MatrixXd J(Kz, 1 + Kx);
  // 5.1 dsX: derivative of sX with respect to beta
  Eigen::MatrixXd dsX(Eigen::MatrixXd::Zero(n, Kx)); 
  dsX(nIso, Eigen::all) = -sX(nIso, Eigen::all) / (1 + theta(1)); 
  dsX = dsX(sel, Eigen::all);
  
  // 5.2 Jacobian
  // J.col(0) = Z' * (- Gy / (1 + beta)^2 - d sX / dbeta * theta.segment(2, Kx))
  J.col(0) = Z(sel, Eigen::all).transpose() * (-Gy(sel) / pow(1 + theta(1), 2) - dsX * theta.segment(2, Kx));
  // J.col(1 and +) = Z' * sX
  J.block(0, 1, Kz, Kx) = Z(sel, Eigen::all).transpose() * sX(sel, Eigen::all);
  // Normalization 
  J /= S;
  
  // 6 Compute variance: Var = inv(J'WJ) J'WVmW'J inv(J'WJ) / S
  Eigen::MatrixXd JtW(J.transpose() * W);
  Eigen::MatrixXd bread_inv((JtW * J).inverse());
  Eigen::MatrixXd meat(JtW * Vm * JtW.transpose());
  Eigen::MatrixXd Var2(bread_inv * meat * bread_inv.transpose() / S);
  Eigen::MatrixXd Var(Eigen::MatrixXd::Zero(2 + Kx, 2 + Kx));
  Var.block(1, 1, 1 + Kx, 1 + Kx) = Var2;
  
  // output
  strCov out;
  out.estimate  = theta;
  out.cov       = Var;
  out.serr      = serr;
  out.serriso   = serriso;
  out.serrniso  = serrniso;
  out.residuals = uneta;
  out.testrho   = std::numeric_limits<double>::quiet_NaN();
  return out;
}


////////////////////////////////////////////////////////////////////////
//////////////////// Main function to call from R //////////////////////
///////////// It computes the estimates and the cov matrice ////////////
////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List fCESMain(const Rcpp::Nullable<double> setrho,
                    const Eigen::MatrixXd& X,
                    const Eigen::VectorXd& y,
                    const Eigen::VectorXd& z,
                    const std::vector<Eigen::ArrayXXd>& G,
                    const std::vector<Eigen::ArrayXi>& friendindex,
                    const Eigen::ArrayXi& cumsn,
                    const Eigen::ArrayXi& frzeroy,
                    const Eigen::ArrayXi& frzeroz,
                    const std::vector<Eigen::ArrayXi>& lIso, //in selection
                    const std::vector<Eigen::ArrayXi>& lnIso,//in selection
                    const Eigen::ArrayXi& Iso,               //lIso in vector
                    const Eigen::ArrayXi& nIso,              //lnIso in vector
                    const Eigen::ArrayXi& nvec,
                    const Eigen::ArrayXXd& yFMiMa,
                    const Eigen::ArrayXXd& zFMiMa,
                    const Eigen::ArrayXi& idXiso, // Index for Xiso
                    const Eigen::ArrayXi& idXniso,// Index for Xniso
                    const Eigen::ArrayXi& sel,
                    const int& n,
                    const int& Kx,
                    const int& S,
                    const int& HACn,
                    const int& dfiso,
                    const int& dfniso,
                    const bool& FE,
                    const double& rhomin,
                    const double& rhomax,
                    const int& maxit,
                    const double& eps_f,
                    const double& eps_g,
                    const unsigned int& nthread){
  strListEstim outlest;
  strEstim outest;
  strCov outcov;
  
  if (setrho.isNotNull()) { // for a fixed value for rho
    int Kiso(idXiso.size()), Kniso(idXniso.size()), Kx(X.cols()), n(X.rows()), 
    S(G.size()), Kz(1 + Kiso + Kniso), Ktheta(2 + Kx);
    double rho = Rcpp::as<double>(setrho);
    bool deriv(false);
    
    // 0.1 data
    Eigen::MatrixXd data = fCESdata(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz,
                                    lIso, lnIso, nvec, yFMiMa, zFMiMa, n, Kx, S,
                                    rho, FE, deriv, nthread);
    // this is how data are organized in data
    // X: 0 to Kx - 1
    // y: Kx
    // Gy: Kx + 1
    // Gz: Kx + 2
    // dGz: Kx + 3
    // dGy: Kx + 4
    // ddGz: Kx + 5

    // 0.2 Instruments
    Eigen::MatrixXd Z(Eigen::MatrixXd::Zero(n, Kz));
    Z(nIso, 0)           = data(nIso, Kx + 2);
    if (Kiso > 0) { // if there are isolated
      Z(Iso, Eigen::seqN(1, Kiso)) = data(Iso, idXiso);
    }
    Z(nIso, Eigen::seqN(1 + Kiso, Kniso)) = data(nIso, idXniso);
    
    // 0.3. Weight
    Eigen::MatrixXd W = (Z(sel, Eigen::all).transpose() * Z(sel, Eigen::all) / S).inverse();
    
    // 1. First estimation
    outest = fCESestim_rho_1ins(rho, data.col(Kx), data.col(Kx + 1), data.block(0, 0, n, Kx),
                                Z, nIso, W, sel, n, S, Ktheta, maxit, eps_f, eps_g);
    
    // 2. Optimal weight
    W = fCESWeight_1ins(outest.theta, data.col(Kx), data.col(Kx + 1), data.block(0, 0, n, Kx),
                        Z, Iso, nIso, lIso, lnIso, sel, n, S, Kx, Kz, HACn, dfiso, dfniso);
    
    // 3. Second estimation
    outest = fCESestim_rho_1ins(rho, data.col(Kx), data.col(Kx + 1), data.block(0, 0, n, Kx),
                                Z, nIso, W, sel, n, S, Ktheta, maxit, eps_f, eps_g);
    
    // 4. Covariance
    outcov = fCESparms_rho(outest.theta, data.col(Kx), data.col(Kx + 1), data.block(0, 0, n, Kx),
                           Z, Iso, nIso, lIso, lnIso, W, sel, n, S, Kx, Kz, HACn, dfiso, dfniso);
    
  } else { // rho is variable
    int Kiso(idXiso.size()), Kniso(idXniso.size()), Kz(2 + Kiso + Kniso), S(G.size());
    Eigen::MatrixXd W = S * Eigen::MatrixXd::Identity(Kz, Kz);

    // 1. Fist estimation
    // Warning: set the range of rho properly before this
    outest = fCESestim(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz, lIso,
                       lnIso, Iso, nIso, nvec, yFMiMa, zFMiMa, W, idXiso, idXniso,
                       sel, rhomin, rhomax, FE, nthread, maxit, eps_f, eps_g);
    
    // 2. Optimal weight
    W = fCESWeight_2ins(outest.theta, X, y, z, G, friendindex, cumsn, frzeroy, frzeroz,
                        lIso, lnIso, Iso, nIso, nvec, yFMiMa, zFMiMa, W, idXiso, idXniso,
                        sel, FE, HACn, dfiso, dfniso, nthread);
    
    // 3. Second estimation
    outest = fCESestim(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz, lIso,
                       lnIso, Iso, nIso, nvec, yFMiMa, zFMiMa, W, idXiso, idXniso,
                       sel, rhomin, rhomax, FE, nthread, maxit, eps_f, eps_g);
    
    // 4. Covariance
    outcov = fCESparms(outest.theta, X, y, z, G, friendindex, cumsn, frzeroy, frzeroz, 
                       lIso, lnIso, Iso, nIso, nvec, yFMiMa, zFMiMa, W, idXiso, idXniso,
                       sel, FE, HACn, dfiso, dfniso, nthread);
  }
  
  return Rcpp::List::create(Rcpp::_["estimate"]      = outcov.estimate,
                            Rcpp::_["cov"]           = outcov.cov,
                            Rcpp::_["serr"]          = outcov.serr,
                            Rcpp::_["serr_iso"]      = outcov.serriso,
                            Rcpp::_["serr_niso"]     = outcov.serrniso,
                            Rcpp::_["Testrho"]       = outcov.testrho,
                            Rcpp::_["objective"]     = outest.objective,  
                            Rcpp::_["gradient"]      = outest.gradient,  
                            Rcpp::_["status"]        = outest.status,
                            Rcpp::_["unscale.resid"] = outcov.residuals);
}