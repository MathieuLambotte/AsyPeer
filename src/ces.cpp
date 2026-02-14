//////////////////////////////////
/// ACKNOWLEDGMENT
/// Taken from QuantilePeer and modified according to the MIT License; QuantilePeer is cited.
/// Houndetoungan A (2025). QuantilePeer: Quantile Peer Effect Models. 
/// doi:10.32614/CRAN.package.QuantilePeer, R package version 0.0.1,
//////////////////////////////////

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#ifdef _OPENMP 
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif
// #define NDEBUG
// #include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;

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
////////////////// returns the value of the objective ///////////////////
////////////////// function of beta, given rho fixed //////////////////// 
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double fCESobjrho(const double& beta,
                  const Eigen::VectorXd& y,
                  const Eigen::VectorXd& Gy,
                  const Eigen::MatrixXd& X,
                  const Eigen::MatrixXd& Z,
                  const Eigen::ArrayXi& nIso,              //lnIso in vector
                  const Eigen::MatrixXd& W,
                  const Eigen::ArrayXi& sel,
                  const int& n,
                  const int& S) {
  
  double theta0 = beta / (1.0 + beta);
  // 1. Scale X for non-isolated
  Eigen::MatrixXd sX = X;
  sX(nIso, Eigen::all) /= (1.0 + beta);
  
  // 2. Closed-form GMM parameters: 
  Eigen::VectorXd u = y - theta0 * Gy;
  Eigen::MatrixXd ZtsX(Z(sel, Eigen::all).transpose() * sX(sel, Eigen::all)); // kz x ksX: Z'sX
  Eigen::VectorXd Ztu(Z(sel, Eigen::all).transpose() * u(sel));  // kz x 1: Z'u
  Eigen::MatrixXd sXtZW(ZtsX.transpose() * W);   //  ksX x kz: sX'Z W
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Adec(sXtZW * ZtsX); // cholesky decomposition for A = sX'Z W Z' sX
  Eigen::VectorXd b(sXtZW * Ztu); // ksX x 1: sX'Z W Z'u
  Eigen::VectorXd phi(Adec.solve(b));
  
  // 3. residual (but scaled for nonisolated)
  Eigen::VectorXd eta = u - sX * phi;
  
  // 4. Objective function
  Eigen::VectorXd mom(Z(sel, Eigen::all).transpose() * eta(sel)); // Not divised by S to avoid unreached precision
  Eigen::VectorXd Wmom = W * mom;
  return mom.dot(Wmom);
}


/////////////////////////////////////////////////////////////////////////
///////////////////// full parameter for fixed rho //////////////////////
///////////////////////////  ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Eigen::VectorXd fCESparmrho(const double& beta,
                            const double& rho,
                            const Eigen::VectorXd& y,
                            const Eigen::VectorXd& Gy,
                            const Eigen::MatrixXd& X,
                            const Eigen::MatrixXd& Z,
                            const Eigen::ArrayXi& nIso,              //lnIso in vector
                            const Eigen::MatrixXd& W,
                            const Eigen::ArrayXi& sel,
                            const int& n,
                            const int& S) {
  
  double theta0 = beta / (1.0 + beta);
  // 1. Scale X for non-isolated
  Eigen::MatrixXd sX = X;
  sX(nIso, Eigen::all) /= (1.0 + beta);
  
  // 2. Closed-form GMM parameters: 
  Eigen::VectorXd u = y - theta0 * Gy;
  Eigen::MatrixXd ZtsX(Z(sel, Eigen::all).transpose() * sX(sel, Eigen::all)); // kz x ksX: Z'sX
  Eigen::VectorXd Ztu(Z(sel, Eigen::all).transpose() * u(sel));  // kz x 1: Z'u
  Eigen::MatrixXd sXtZW(ZtsX.transpose() * W);   //  ksX x kz: sX'Z W
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Adec(sXtZW * ZtsX); // cholesky decomposition for A = sX'Z W Z' sX
  Eigen::VectorXd b(sXtZW * Ztu); // ksX x 1: sX'Z W Z'u
  Eigen::VectorXd phi(Adec.solve(b));
  
  // 3. Full parameters
  Eigen::VectorXd theta(2 + phi.size());
  theta << rho, beta, phi;
  
  return theta;
}


/////////////////////////////////////////////////////////////////////////
////////////////// returns the value of the objective ///////////////////
/////////////////////// function of rho beta //////////////////////////// 
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double fCESobj(const Eigen::VectorXd& theta01,
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
               const double& rhomin,
               const double& rhomax) {
  int Kiso(idXiso.size()), Kniso(idXniso.size()), Kx(X.cols()), n(X.rows()), 
  S(G.size()), Kz(2 + Kiso + Kniso);
  bool deriv(true);
  
  double rho(theta01(0)), beta(theta01(1)), 
  parm1 = beta / (1.0 + beta); // Coefficient of Gy
  
  if ((rho <= rhomin) || (rho >= rhomax) || (beta <= -0.5)) {
    return 1e250;
  }
  // 0. data
  Eigen::MatrixXd data = fCESdata(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz,
                                  lIso, lnIso, nvec, yFMiMa, zFMiMa, n, Kx, S,
                                  rho, FE, deriv, 1);
  
  // 1. Instruments
  Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(n, Kz);
  Z(nIso, 0) = data(nIso, Kx + 2); //Gz
  Z(nIso, 1) = data(nIso, Kx + 3); //dGz
  if (Kiso > 0) { // if there are isolated
    Z(Iso, Eigen::seqN(2, Kiso)) = data(Iso, idXiso);
  }
  Z(nIso, Eigen::seqN(2 + Kiso, Kniso)) = data(nIso, idXniso);
  
  // 2. Scale X for non-isolated
  Eigen::MatrixXd sX = data.block(0, 0, n, Kx);
  sX(nIso, Eigen::all) /= (1.0 + beta);
  
  // 3. Closed-form GMM parameters: 
  Eigen::VectorXd u = data.col(Kx) - parm1 * data.col(Kx + 1); // u = y - parm1 * Gy
  Eigen::MatrixXd ZtsX(Z(sel, Eigen::all).transpose() * sX(sel, Eigen::all)); // kz x ksX: Z'sX
  Eigen::VectorXd Ztu(Z(sel, Eigen::all).transpose() * u(sel));  // kz x 1: Z'u
  Eigen::MatrixXd sXtZW(ZtsX.transpose() * W);   //  ksX x kz: sX'Z W
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Adec(sXtZW * ZtsX); // cholesky decomposition for A = sX'Z W Z' sX
  Eigen::VectorXd b(sXtZW * Ztu); // ksX x 1: sX'Z W Z'u
  Eigen::VectorXd phi(Adec.solve(b));
  
  // 4. residual (but scaled for nonisolated)
  Eigen::VectorXd eta = u - sX * phi;
  
  // 4. Objective function
  Eigen::VectorXd mom(Z(sel, Eigen::all).transpose() * eta(sel)); // Not divised by S to avoid unreached precision
  Eigen::VectorXd Wmom = W * mom;
  return mom.dot(Wmom);
}


/////////////////////////////////////////////////////////////////////////
////////////////// full parameter with flexible rho /////////////////////
///////////////////////////  ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Eigen::VectorXd fCESparm(const Eigen::VectorXd& theta01,
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
                         const bool& FE) {
  int Kiso(idXiso.size()), Kniso(idXniso.size()), Kx(X.cols()), n(X.rows()), 
  S(G.size()), Kz(2 + Kiso + Kniso);
  bool deriv(true);
  
  double rho(theta01(0)), beta(theta01(1)), 
  parm1 = beta / (1.0 + beta); // Coefficient of Gy
  
  Eigen::MatrixXd data = fCESdata(X, y, z, G, friendindex, cumsn, frzeroy, frzeroz,
                                  lIso, lnIso, nvec, yFMiMa, zFMiMa, n, Kx, S,
                                  rho, FE, deriv, 1);
  
  // 1. Instruments
  Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(n, Kz);
  Z(nIso, 0) = data(nIso, Kx + 2); //Gz
  Z(nIso, 1) = data(nIso, Kx + 3); //dGz
  if (Kiso > 0) { // if there are isolated
    Z(Iso, Eigen::seqN(2, Kiso)) = data(Iso, idXiso);
  }
  Z(nIso, Eigen::seqN(2 + Kiso, Kniso)) = data(nIso, idXniso);
  
  // 2. Scale X for non-isolated
  Eigen::MatrixXd sX = data.block(0, 0, n, Kx);
  sX(nIso, Eigen::all) /= (1.0 + beta);
  
  // 3. Closed-form GMM parameters: 
  Eigen::VectorXd u = data.col(Kx) - parm1 * data.col(Kx + 1); // u = y - parm1 * Gy
  Eigen::MatrixXd ZtsX(Z(sel, Eigen::all).transpose() * sX(sel, Eigen::all)); // kz x ksX: Z'sX
  Eigen::VectorXd Ztu(Z(sel, Eigen::all).transpose() * u(sel));  // kz x 1: Z'u
  Eigen::MatrixXd sXtZW(ZtsX.transpose() * W);   //  ksX x kz: sX'Z W
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Adec(sXtZW * ZtsX); // cholesky decomposition for A = sX'Z W Z' sX
  Eigen::VectorXd b(sXtZW * Ztu); // ksX x 1: sX'Z W Z'u
  Eigen::VectorXd phi(Adec.solve(b));
  
  // 3. Full parameters
  Eigen::VectorXd theta(2 + phi.size());
  theta << rho, beta, phi;
  
  return theta;
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// GMM Weight /////////////////////////////////
//////////////////// using one instrument: Gz ///////////////////////////
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
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
//////////////////////////// GMM Weight /////////////////////////////////
///////////////// using two instruments: Gz and dGz /////////////////////
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
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
//////////////////// full parameters and variances //////////////////////
/////////////////////////// with a fixed rho ////////////////////////////
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List fCESparmCovrho(const Eigen::VectorXd& theta,
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
  return Rcpp::List::create(Rcpp::_["estimate"]      = theta,
                            Rcpp::_["cov"]           = Var,
                            Rcpp::_["serr"]          = serr,
                            Rcpp::_["serr_iso"]      = serriso,
                            Rcpp::_["serr_niso"]     = serrniso,
                            Rcpp::_["Testrho"]       = std::numeric_limits<double>::quiet_NaN(),
                            Rcpp::_["unscale.resid"] = uneta);
}


/////////////////////////////////////////////////////////////////////////
//////////////////// full parameters and variances //////////////////////
//////////////////////// with a flexible rho ////////////////////////////
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List fCESparmCov(const Eigen::VectorXd& theta,
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
  return Rcpp::List::create(Rcpp::_["estimate"]      = theta,
                            Rcpp::_["cov"]           = Var,
                            Rcpp::_["serr"]          = serr,
                            Rcpp::_["serr_iso"]      = serriso,
                            Rcpp::_["serr_niso"]     = serrniso,
                            Rcpp::_["Testrho"]       = testrho,
                            Rcpp::_["unscale.resid"] = uneta);
}


