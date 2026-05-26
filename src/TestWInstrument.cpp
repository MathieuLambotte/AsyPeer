// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace std;

// This function computes 
// Conditional F-statistic of Sanderson and Windmeijer (2016)
//[[Rcpp::export]]
Eigen::ArrayXXd fsw_cpp(const Eigen::MatrixXd& y,
                        const Eigen::MatrixXd& Z,
                        const int& Kexo) {
  
  int n = y.rows();
  int Kendo = y.cols();
  int Kins  = Z.cols() - Kexo;
  
  Eigen::ArrayXXd out(Kendo, 4);
  
  for (int k = 0; k < Kendo; ++k) {
    
    // ========= Step 1: build variables =========
    Eigen::VectorXd yk = y.col(k);
    
    Eigen::MatrixXd Yother(n, Kendo - 1);
    int c = 0;
    for (int j = 0; j < Kendo; ++j) {
      if (j != k) Yother.col(c++) = y.col(j);
    }
    
    
    Eigen::MatrixXd Xexo = Z.leftCols(Kexo);
    Eigen::MatrixXd Zins = Z.rightCols(Kins);
    
    // ========= Step 2: conditional IV (# EXACT FIX HERE) =========
    Eigen::MatrixXd X(n, Yother.cols() + Xexo.cols());
    X << Yother, Xexo;
    
    Eigen::MatrixXd Zfull(n, Zins.cols() + Xexo.cols());
    Zfull << Zins, Xexo;
    
    // First stage
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrZ(Zfull);
    Eigen::MatrixXd Xhat = Zfull * qrZ.solve(X);
    
    // Second stage
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrX(Xhat);
    Eigen::VectorXd beta = qrX.solve(yk);
    
    // ✅ CORRECT residual (CRITICAL FIX)
    Eigen::VectorXd condres = yk - Xhat * beta;
    
    // ========= Step 3: unrestricted regression =========
    Eigen::MatrixXd X_u(n, Zins.cols() + Xexo.cols());
    X_u << Zins, Xexo;
    
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrU(X_u);
    Eigen::VectorXd beta_u = qrU.solve(condres);
    Eigen::VectorXd res_u  = condres - X_u * beta_u;
    
    double RSS_u = res_u.squaredNorm();
    
    // ========= Step 4: restricted regression =========
    Eigen::MatrixXd X_r = Xexo;
    
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrR(X_r);
    Eigen::VectorXd beta_r = qrR.solve(condres);
    Eigen::VectorXd res_r  = condres - X_r * beta_r;
    
    double RSS_r = res_r.squaredNorm();
    
    // ========= Step 5: F-stat =========
    int df1_raw = Kins;
    int df2     = n - X_u.cols();
    
    double F =
      ((RSS_r - RSS_u) / df1_raw) /
        (RSS_u / df2);
    
    // ========= Step 6: SW correction =========
    int df1 = df1_raw - (Kendo - 1);
    
    double Fsw = (F * df1_raw) / df1;
    
    // ========= Step 7: output =========
    out(k,0) = Fsw;
    out(k,1) = df1;
    out(k,2) = df2;
    out(k,3) = 1.0 - R::pf(Fsw, df1, df2, 1, 0);
  }
  
  return out;
}

