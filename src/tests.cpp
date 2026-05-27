// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace std;

// Computes sqrt of matrices
Eigen::MatrixXd matrixSqrt(const Eigen::MatrixXd& A) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  Eigen::VectorXd sqrt_evals = es.eigenvalues().array().sqrt();
  return es.eigenvectors() * sqrt_evals.asDiagonal() * es.eigenvectors().transpose();
}

// This function estimates F stats and predict endogenous variables
//[[Rcpp::export]]
Eigen::ArrayXXd fFstat(const Eigen::MatrixXd& endo,
                       const Eigen::MatrixXd& Z,
                       const int& K, 
                       const Eigen::ArrayXi& cumsn,
                       const bool& cluster, 
                       const int& nthread) {
  int n      = endo.rows();
  int Kendo  = endo.cols();
  int L      = Z.cols() - K;
  int ngroup = cumsn.size() - 1;
  int df1    = L;
  int df2    = n - L - K;
  
  // Partialling out X
  Eigen::MatrixXd ry, rZ;
  {
    Eigen::MatrixXd X(Z.leftCols(K));
    Eigen::ColPivHouseholderQR <Eigen::MatrixXd> qr(X);
    ry = endo - X * qr.solve(endo);
    rZ = Z.rightCols(L) - X * qr.solve(Z.rightCols(L));
  }
  
  Eigen::MatrixXd ZZ(rZ.transpose() * rZ);
  Eigen::MatrixXd iZZ((rZ.transpose() * rZ).ldlt().solve(Eigen::MatrixXd::Identity(L, L)));
  Eigen::MatrixXd Zy(rZ.transpose() * ry);
  
  // estimator and residuals
  Eigen::MatrixXd Pi(iZZ * Zy);
  Eigen::MatrixXd e(ry - rZ * Pi);
  
  // F and pvalue
  Eigen::ArrayXd F(Kendo);
  Eigen::ArrayXd prob(Kendo);
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int k = 0; k < Kendo; ++ k) {
    Eigen::MatrixXd V(Eigen::MatrixXd::Zero(L, L));
    
    if (cluster) {
      for (int r(0); r < ngroup; ++ r) {
        int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
        Eigen::VectorXd Ze(rZ(Eigen::seq(n1, n2), Eigen::all).transpose() * e(Eigen::seq(n1, n2), k).matrix());
        V += Ze * Ze.transpose();
      }
    } else {
      Eigen::MatrixXd Ze((rZ.array().colwise() * e.col(k).array()));
      V = Ze.transpose() * Ze;
    }
    
    F(k)    = Pi.col(k).dot((iZZ * V * iZZ).colPivHouseholderQr().solve(Pi.col(k))) / df1;
    prob(k) = R::pf(F(k), df1, df2, 0, 0);
  }
#else
  for (int k = 0; k < Ky; ++ k) {
    Eigen::MatrixXd V(Eigen::MatrixXd::Zero(L, L));
    
    if (cluster) {
      for (int r(0); r < ngroup; ++ r) {
        int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
        Eigen::VectorXd Ze(rZ(Eigen::seq(n1, n2), Eigen::all).transpose() * e(Eigen::seq(n1, n2), k).matrix());
        V += Ze * Ze.transpose();
      }
    } else {
      Eigen::MatrixXd Ze((rZ.array().colwise() * e.col(k).array()));
      V = Ze.transpose() * Ze;
    }
    
    F(k)    = Pi.col(k).dot((iZZ * V * iZZ).colPivHouseholderQr().solve(Pi.col(k))) / df1;
    prob(k) = R::pf(F(k), df1, df2, 0, 0);
  }
#endif
  
  Eigen::ArrayXXd out(Kendo, 4);
  out.col(0).setConstant(df1);
  out.col(1).setConstant(df2);
  out.col(2) = F;
  out.col(3) = prob;
  return out; 
}

// This function computes KP stat
//[[Rcpp::export]]
Eigen::ArrayXXd fKPstat(const Eigen::MatrixXd& endo,
                       const Eigen::MatrixXd& Z,
                       const int& K, 
                       const Eigen::ArrayXi& cumsn,
                       const bool& cluster) {
  int n      = endo.rows();
  int Kendo  = endo.cols();
  int L      = Z.cols() - K;
  int ngroup = cumsn.size() - 1;
  
  // Partialling out X
  Eigen::MatrixXd ry, rZ;
  {
    Eigen::MatrixXd X(Z.leftCols(K));
    Eigen::ColPivHouseholderQR <Eigen::MatrixXd> qr(X);
    ry = endo - X * qr.solve(endo);
    rZ = Z.rightCols(L) - X * qr.solve(Z.rightCols(L));
  }
  
  Eigen::MatrixXd ZZ(rZ.transpose() * rZ);
  Eigen::MatrixXd iZZ((rZ.transpose() * rZ).ldlt().solve(Eigen::MatrixXd::Identity(L, L)));
  Eigen::MatrixXd Zy(rZ.transpose() * ry);
  
  // estimator and residuals
  Eigen::MatrixXd Pi(Zy.transpose() * iZZ);
  Eigen::MatrixXd eps(ry - rZ * Pi.transpose());
  
  // vec(Ze)
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(L * Kendo, L * Kendo));
  for (int s1(0); s1 < L; ++ s1) {
    for (int s2(0); s2 < Kendo; ++ s2) {
      R(s1 * Kendo + s2, s2 * L + s1) = 1;
    }
  }
  
  Eigen::MatrixXd vecZe(n, L*Kendo);
  for (int s(0); s < Kendo; ++ s) {
    vecZe.block(0, s*L, n, L) = rZ.array().colwise()*eps.col(s).array();
  }
  
  // Variance of vec(Ze), covyy and covz
  Eigen::MatrixXd VvecZe(Eigen::MatrixXd::Zero(L*Kendo, L*Kendo));
  if (cluster) {
    for (int r(0); r < ngroup; ++ r) {
      int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
      Eigen::VectorXd tp(vecZe(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix().transpose());
      VvecZe += tp * tp.transpose();
    }
  } else {
    VvecZe = vecZe.transpose() * vecZe;
  } 
  VvecZe /= n;
  
  // Variance of pi
  Eigen::MatrixXd H(R * Eigen::kroneckerProduct(Eigen::MatrixXd::Identity(Kendo, Kendo), iZZ));
  Eigen::MatrixXd varpi(n * H * VvecZe * H.transpose()); 
  
  // normalisation
  Eigen::MatrixXd covy = eps.transpose() * eps / n;
  Eigen::LLT<Eigen::MatrixXd> tpF(ZZ / n);// tpG(covy);
  // Eigen::MatrixXd G(tpG1.matrixL().inverse()); 
  Eigen::MatrixXd F(tpF.matrixL());
  Eigen::LLT<Eigen::MatrixXd> tpG1(covy);
  Eigen::MatrixXd tpG2 = tpG1.matrixL();
  Eigen::MatrixXd G = tpG2.triangularView<Eigen::Lower>()
                          .solve(Eigen::MatrixXd::Identity(Kendo, Kendo));
  
  // Theta and its variance
  Eigen::MatrixXd Theta(G * Pi * F.transpose());
  Eigen::VectorXd theta(Theta.reshaped(L*Kendo, 1));
  Eigen::MatrixXd FG(Eigen::kroneckerProduct(F, G));
  Eigen::MatrixXd vartheta(FG * varpi * FG.transpose());
  // cout << F << endl;
  // cout << G << endl;
  // cout << Pi << endl;
  // cout << vartheta << endl;
  // Until this, replicate using bootstrap
  
  // SDV decomposition of Theta
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Theta, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::MatrixXd U = svd.matrixU(); //Kendo * Kendo
  Eigen::VectorXd d = svd.singularValues();
  Eigen::MatrixXd ddiag = d.asDiagonal();
  Eigen::MatrixXd D(Kendo, L);
  D << ddiag, Eigen::MatrixXd::Zero(Kendo, L - Kendo); //L*Kendo
  Eigen::MatrixXd V = svd.matrixV(); //Kendo * Kendo
  
  //U12, U22, V12, V22
  int q(Kendo - 1);
  Eigen::MatrixXd U12(U.block(0, q, q, Kendo - q));
  Eigen::MatrixXd U22(U.block(q, q, Kendo - q, Kendo - q));
  Eigen::MatrixXd V12(V.block(0, q, q, L - q));
  Eigen::MatrixXd V22(V.block(q, q, L - q, L - q));
  
  // Aqper and Bqper
  Eigen::MatrixXd U12U22(Kendo, Kendo - q), V12V22(L, L - q);
  U12U22 << U12, U22;
  V12V22 << V12, V22;
  
  Eigen::MatrixXd Aper(U12U22 * U22.colPivHouseholderQr().solve(matrixSqrt(U22 * U22.transpose())));
  Eigen::MatrixXd Bper((V12V22 * V22.colPivHouseholderQr().solve(matrixSqrt(V22 * V22.transpose()))).transpose());
  
  // lambda and its varianve
  Eigen::MatrixXd BAper(Eigen::kroneckerProduct(Bper, Aper.transpose()));
  Eigen::VectorXd lambda (BAper * theta);
  Eigen::MatrixXd varlambda (BAper * vartheta * BAper.transpose());
  
  // statistic
  double stat = lambda.dot(varlambda.colPivHouseholderQr().solve(lambda));
  
  // Degree of freedom
  int df = (Kendo - q)*(L - q);
  
  // pvalue
  double prob = R::pchisq(stat, df, 0, 0);
  
  // output
  Eigen::ArrayXXd out(1, 4);
  out << df, std::numeric_limits<double>::quiet_NaN(),
         stat, prob;
  
  return out;
}

// This function computes 
// Conditional F-statistic of Sanderson and Windmeijer (2016)
//[[Rcpp::export]]
Eigen::ArrayXXd fswstat(const Eigen::MatrixXd& endo,
                        const Eigen::MatrixXd& Z,
                        const int& K,
                        const int& nthread) {
  
  int n = endo.rows();
  int Kendo = endo.cols();
  int L     = Z.cols() - K;
  
  Eigen::ArrayXXd out(Kendo, 4);
  
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int k = 0; k < Kendo; ++k) {
    
    // ========= Step 1: build variables =========
    Eigen::VectorXd yk = endo.col(k);
    
    Eigen::MatrixXd Yother(n, Kendo - 1);
    int c = 0;
    for (int j = 0; j < Kendo; ++j) {
      if (j != k) Yother.col(c++) = endo.col(j);
    }
    
    
    Eigen::MatrixXd Xexo = Z.leftCols(K);
    
    // ========= Step 2: conditional IV =========
    Eigen::MatrixXd X(n, Kendo - 1 + K);
    X << Yother, Xexo;
    
    // First stage
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrZ(Z);
    Eigen::MatrixXd Xhat = Z * qrZ.solve(X);
    
    // Second stage
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrX(Xhat);
    Eigen::VectorXd beta = qrX.solve(yk);
    
    // residual 
    Eigen::VectorXd condres = yk - X * beta;
    
    // ========= Step 3: unrestricted regression =========
    Eigen::VectorXd beta_u = qrZ.solve(condres);
    Eigen::VectorXd res_u  = condres - Z * beta_u;
    
    double RSS_u = res_u.squaredNorm();
    
    // ========= Step 4: restricted regression =========
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrR(Xexo);
    Eigen::VectorXd beta_r = qrR.solve(condres);
    Eigen::VectorXd res_r  = condres - Xexo * beta_r;
    
    double RSS_r = res_r.squaredNorm();
    
    // ========= Step 5: F-stat =========
    int df1_raw = L;
    int df2     = n - (K + L);
    
    double F    = ((RSS_r - RSS_u) / df1_raw) / (RSS_u / df2);
    
    // ========= Step 6: SW correction =========
    int df1 = df1_raw - (Kendo - 1);
    
    double Fsw = (F * df1_raw) / df1;
    
    // ========= Step 7: output =========
    out(k,0) = df1;
    out(k,1) = df2;
    out(k,2) = Fsw;
    out(k,3) = R::pf(Fsw, df1, df2, 0, 0);
  }
#else
  for (int k = 0; k < Kendo; ++k) {
    
    // ========= Step 1: build variables =========
    Eigen::VectorXd yk = endo.col(k);
    
    Eigen::MatrixXd Yother(n, Kendo - 1);
    int c = 0;
    for (int j = 0; j < Kendo; ++j) {
      if (j != k) Yother.col(c++) = endo.col(j);
    }
    
    
    Eigen::MatrixXd Xexo = Z.leftCols(K);
    
    // ========= Step 2: conditional IV =========
    Eigen::MatrixXd X(n, Kendo - 1 + K);
    X << Yother, Xexo;
    
    // First stage
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrZ(Z);
    Eigen::MatrixXd Xhat = Z * qrZ.solve(X);
    
    // Second stage
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrX(Xhat);
    Eigen::VectorXd beta = qrX.solve(yk);
    
    // residual 
    Eigen::VectorXd condres = yk - X * beta;
    
    // ========= Step 3: unrestricted regression =========
    Eigen::VectorXd beta_u = qrZ.solve(condres);
    Eigen::VectorXd res_u  = condres - Z * beta_u;
    
    double RSS_u = res_u.squaredNorm();
    
    // ========= Step 4: restricted regression =========
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrR(Xexo);
    Eigen::VectorXd beta_r = qrR.solve(condres);
    Eigen::VectorXd res_r  = condres - Xexo * beta_r;
    
    double RSS_r = res_r.squaredNorm();
    
    // ========= Step 5: F-stat =========
    int df1_raw = L;
    int df2     = n - (K + L);
    
    double F    = ((RSS_r - RSS_u) / df1_raw) / (RSS_u / df2);
    
    // ========= Step 6: SW correction =========
    int df1 = df1_raw - (Kendo - 1);
    
    double Fsw = (F * df1_raw) / df1;
    
    // ========= Step 7: output =========
    out(k,0) = df1;
    out(k,1) = df2;
    out(k,2) = Fsw;
    out(k,3) = R::pf(Fsw, df1, df2, 0, 0);
  }
#endif
  
  return out;
}

