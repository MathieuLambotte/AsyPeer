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

using namespace Rcpp;
using namespace Eigen;
using namespace std;

// This function set nthreads
//[[Rcpp::export]]
int fnthreads(const int& nthread) {
#ifdef _OPENMP
  return nthread;
#else
  return 1;
#endif
}


// Compute statistics such as ybarh, ybarl, gil, gih
//[[Rcpp::export]]
List highlowstat1(const Eigen::MatrixXd& X,
                  const std::vector<Eigen::ArrayXXd>& G,
                  const Eigen::ArrayXi& cumsn,
                  const Eigen::ArrayXi& nvec,
                  const int& ngroup,
                  const unsigned int& nthread) {
  int n(X.rows()), K(X.cols());
  Eigen::ArrayXXd Xb(n, K), Xbh(n, K), Xbl(n, K), gl(n, K), gh(n, K);
  Eigen::ArrayXd g(n);
  
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int m = 0; m < ngroup; ++ m) {
    int n1(cumsn(m)); // Where the group starts in X.
    int nm(nvec(m));
    Eigen::ArrayXXd Xm       = X.block(n1, 0, nm, K);
    Eigen::ArrayXd Gmsum(G[m].rowwise().sum()); g.segment(n1, nm)  = Gmsum;
    Eigen::ArrayXXd GmX(G[m].matrix() * Xm.matrix()); Xb.block(n1, 0, nm, K) = GmX;
    for (int i(0); i < nm; ++ i) {
      for(int k(0); k < K; ++ k) {
        Eigen::ArrayXd Gmhi = (Xm.col(k) > Xm(i, k)).select(G[m].row(i).transpose(), 0);
        gh(n1 + i, k)  = Gmhi.sum();
        Xbh(n1 + i, k) = (Gmhi * Xm.col(k)).sum();
      }
    }
  }
#else
  for (int m = 0; m < ngroup; ++ m) {
    int n1(cumsn(m)); // Where the group starts in X.
    int nm(nvec(m));
    Eigen::ArrayXXd Xm       = X.block(n1, 0, nm, K);
    Eigen::ArrayXd Gmsum(G[m].rowwise().sum()); g.segment(n1, nm)  = Gmsum;
    Eigen::ArrayXXd GmX(G[m].matrix() * Xm.matrix()); Xb.block(n1, 0, nm, K) = GmX;
    for (int i(0); i < nm; ++ i) {
      for(int k(0); k < K; ++ k) {
        Eigen::ArrayXd Gmhi = (Xm.col(k) > Xm(i, k)).select(G[m].row(i).transpose(), 0);
        gh(n1 + i, k)  = Gmhi.sum();
        Xbh(n1 + i, k) = (Gmhi * Xm.col(k)).sum();
      }
    }
  }
#endif
  for(int k(0); k < K; ++ k) {
    gl.col(k) = g - gh.col(k);
  }
  Xbl = Xb - Xbh;
  return List::create(_["Xbar"] = Xb,
                      _["Xbl"]  = Xbl,
                      _["Xbh"]  = Xbh,
                      _["g"]    = g,
                      _["gl"]   = gl,
                      _["gh"]   = gh);
}


// same function as highlowstat1 but high and low are determined by another variable
//[[Rcpp::export]]
List highlowstat2(const Eigen::VectorXd& y,
                  const Eigen::MatrixXd& X,
                  const std::vector<Eigen::ArrayXXd>& G,
                  const Eigen::ArrayXi& cumsn,
                  const Eigen::ArrayXi& nvec,
                  const int ngroup,
                  const unsigned int& nthread) {
  int n(X.rows()), K(X.cols());
  Eigen::ArrayXd yb(n), ybh(n), ybl(n), gl(n), gh(n), g(n);
  Eigen::ArrayXXd Xb(n, K), Xbh(n, K), Xbl(n, K);
  
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int m = 0; m < ngroup; ++ m) {
    int n1(cumsn(m)); // Where the group starts in X.
    int nm(nvec(m));
    Eigen::ArrayXd ym        = y.segment(n1, nm);
    Eigen::MatrixXd Xm       = X.block(n1, 0, nm, K);
    Eigen::ArrayXd Gmsum(G[m].rowwise().sum()); g.segment(n1, nm)  = Gmsum;
    Eigen::ArrayXd Gmy(G[m].matrix() * ym.matrix()); yb.segment(n1, nm)      = Gmy;
    Eigen::ArrayXXd GmX(G[m].matrix() * Xm); Xb.block(n1, 0, nm, K) = GmX;
    for (int i(0); i < nm; ++ i) {
      Eigen::ArrayXd Gmhi = (ym > ym(i)).select(G[m].row(i).transpose(), 0);
      gh(n1 + i)      = Gmhi.sum();
      ybh(n1 + i)     = (Gmhi * ym.array()).sum();
      Xbh.row(n1 + i) = Gmhi.matrix().transpose() * Xm;
    }
  }
#else
  for (int m = 0; m < ngroup; ++ m) {
    int n1(cumsn(m)); // Where the group starts in X.
    int nm(nvec(m));
    Eigen::ArrayXd ym        = y.segment(n1, nm);
    Eigen::MatrixXd Xm       = X.block(n1, 0, nm, K);
    Eigen::ArrayXd Gmsum(G[m].rowwise().sum()); g.segment(n1, nm)  = Gmsum;
    Eigen::ArrayXd Gmy(G[m].matrix() * ym.matrix()); yb.segment(n1, nm)      = Gmy;
    Eigen::ArrayXXd GmX(G[m].matrix() * Xm); Xb.block(n1, 0, nm, K) = GmX;
    for (int i(0); i < nm; ++ i) {
      Eigen::ArrayXd Gmhi = (ym > ym(i)).select(G[m].row(i).transpose(), 0);
      gh(n1 + i)      = Gmhi.sum();
      ybh(n1 + i)     = (Gmhi * ym.array()).sum();
      Xbh.row(n1 + i) = Gmhi.matrix().transpose() * Xm;
    }
  }
#endif
  gl  = g - gh;
  ybl = yb - ybh;
  Xbl = Xb - Xbh;
  return List::create(_["ybar"] = yb,
                      _["ybl"]  = ybl,
                      _["ybh"]  = ybh,
                      _["Xbar"] = Xb,
                      _["Xbl"]  = Xbl,
                      _["Xbh"]  = Xbh,
                      _["g"]    = g,
                      _["gl"]   = gl,
                      _["gh"]   = gh);
}

// This computes peer averages, with power
//[[Rcpp::export]]
Eigen::MatrixXd peeravgpower(const std::vector<Eigen::MatrixXd>& G,
                             const Eigen::ArrayXXd& V,
                             const Eigen::ArrayXi& cumsn,
                             const Eigen::ArrayXi& nvec,
                             const int& power,
                             const int& nthread) {
  int kV(V.cols()), n(nvec.sum()), ngroup(nvec.size());
  Eigen::MatrixXd out(n, kV * power);
  out.block(0, 0, n, kV) = V;
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int m = 0; m < ngroup; ++m) {
    for (int k = 1; k < power; ++k) {
      out.block(cumsn(m), kV * k, nvec(m), kV) = G[m] * out.block(cumsn(m), kV * (k - 1), nvec(m), kV);
    }
  }
#else
  for (int m = 0; m < ngroup; ++m) {
    for (int k = 1; k < power; ++k) {
      out.block(cumsn(m), kV * k, nvec(m), kV) = G[m] * out.block(cumsn(m), kV * (k - 1), nvec(m), kV);
    }
  }
#endif
  return out;
}

// This function removes columns to obtain full rank matrices
// Taken from the QuantilePeer package
//[[Rcpp::export]]
Eigen::Array<bool, Eigen::Dynamic, 1> fcheckrankEigen(const Eigen::MatrixXd& X, const double& tol = 1e-10) {
  int n(X.rows());
  Eigen::RowVectorXd m(X.colwise().mean());
  Eigen::RowVectorXd s(((X.rowwise() - m).array().square().colwise().sum() / n).sqrt());
  m = (s.array() < tol).select(0, m);
  s = (s.array() < tol).select(1, s);
  Eigen::MatrixXd U((X.rowwise() - m).array().rowwise() / s.array());
  U = U.transpose()*U/n;
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(U);
  Eigen::MatrixXd R(qr.matrixQR().topRows(U.cols()));
  // std::cout<<R.diagonal().transpose()<<std::endl;
  return R.diagonal().array().abs() > tol;
}

// Assigning folds to groups
//[[Rcpp::export]]
Eigen::ArrayXi fassignfold(const Eigen::ArrayXi& ddgroup,
                           const int& nfold) {
  int ngroup(ddgroup.maxCoeff() + 1);
  // Number of pairs per group
  Eigen::ArrayXi nvec(Eigen::ArrayXi::Zero(ngroup));
  for (int i(0); i < ddgroup.size(); ++i) {
    nvec(ddgroup(i)) += 1;
  }
  // Fold for each group
  Eigen::ArrayXi fold(ngroup);
  Eigen::ArrayXi foldsize(Eigen::ArrayXi::Zero(nfold));
  for (int s(0); s < ngroup; ++s) {
    int minsize(foldsize.minCoeff());
    for (int k(0); k < nfold; ++k) {
      if (foldsize(k) == minsize) {
        fold(s)      = k;
        foldsize(k) += nvec(s);
        break;
      }
    }
  }
  return fold(ddgroup);
}

//[[Rcpp::export]]
Eigen::ArrayXXd Demean(const Eigen::ArrayXXd& X,
                       const Eigen::ArrayXi& cumsn,
                       const std::vector<Eigen::ArrayXi>& lIso,
                       const std::vector<Eigen::ArrayXi>& lnIso,
                       const int& nthread){
  int ngroup(cumsn.size() - 1);
  Eigen::ArrayXXd out(X);
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int s = 0; s < ngroup; ++ s) {
    // For isolated
    if (lIso[s].size() > 0) {
      out(lIso[s], Eigen::all).rowwise() -= out(lIso[s], Eigen::all).colwise().mean();
    }
    // For non-isolated
    if (lnIso[s].size() > 0) {
      out(lnIso[s], Eigen::all).rowwise() -= out(lnIso[s], Eigen::all).colwise().mean();
    }
  }
#else
  for (int s = 0; s < ngroup; ++ s) {
    // For isolated
    if (lIso[s].size() > 0) {
      out(lIso[s], Eigen::all).rowwise() -= out(lIso[s], Eigen::all).colwise().mean();
    }
    // For non-isolated
    if (lnIso[s].size() > 0) {
      out(lnIso[s], Eigen::all).rowwise() -= out(lnIso[s], Eigen::all).colwise().mean();
    }
  }
#endif
  return out;
}

//[[Rcpp::export]]
std::vector<Eigen::ArrayXXd> fGnormalise(std::vector<Eigen::ArrayXXd>& G, 
                                         const int& nthread = 1) {
  int S(G.size());
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for(int s = 0; s < S; ++s) {
    Eigen::ArrayXd rowsum((G[s].rowwise().sum() * 1e7).round() / 1e7);
    rowsum = (rowsum > 0).select(rowsum, 1);
    G[s].colwise() /= rowsum;
  }
#else
  for(int s = 0; s < S; ++s) {
    Eigen::ArrayXd rowsum((G[s].rowwise().sum() * 1e7).round() / 1e7);
    rowsum = (rowsum > 0).select(rowsum, 1);
    G[s].colwise() /= rowsum;
  }
#endif
  return G;
}


// This function estimates F stats and predict endogenous variables
//[[Rcpp::export]]
Rcpp::List fFstat(const Eigen::MatrixXd& y,
                  const Eigen::MatrixXd& X,
                  const Eigen::VectorXi& index,
                  const Eigen::ArrayXi& cumsn,
                  const int& HAC, 
                  const int& nthread) {
  int n(y.rows()), K(X.cols()), df1(index.size()), df2(n - K), 
  Ky(y.cols()), ngroup(cumsn.size() - 1);
  
  Eigen::MatrixXd XX(X.transpose()*X), iXX(XX.inverse());
  Eigen::MatrixXd b(iXX * X.transpose() * y);
  Eigen::ArrayXXd e(y - X * b);
  Eigen::VectorXd F(Ky);
#ifdef _OPENMP
 omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int k = 0; k < Ky; ++ k) {
    Eigen::MatrixXd V(Eigen::MatrixXd::Zero(K, K));
    if (HAC <= 2) {
      Eigen::MatrixXd Xe((X.array().colwise()*e.col(k)));
      V = Xe.transpose()*Xe;
    } else {
      for (int r(0); r < ngroup; ++ r) {
        int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
        Eigen::VectorXd tp(X(Eigen::seq(n1, n2), Eigen::all).transpose() * e(Eigen::seq(n1, n2), k).matrix());
        V += tp * tp.transpose();
      }
    }
    Eigen::MatrixXd tp(iXX*V*iXX);
    V    = tp(index, index);
    Eigen::VectorXd bk(b(index, k));
    F(k) = (bk.array() * (V.colPivHouseholderQr().solve(bk)).array()).sum() / df1;
  }
#else
  for (int k = 0; k < Ky; ++ k) {
    Eigen::MatrixXd V(Eigen::MatrixXd::Zero(K, K));
    if (HAC <= 2) {
      Eigen::MatrixXd Xe((X.array().colwise()*e.col(k)));
      V = Xe.transpose()*Xe;
    } else {
      for (int r(0); r < ngroup; ++ r) {
        int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
        Eigen::VectorXd tp(X(Eigen::seq(n1, n2), Eigen::all).transpose() * e(Eigen::seq(n1, n2), k).matrix());
        V += tp * tp.transpose();
      }
    }
    Eigen::MatrixXd tp(iXX*V*iXX);
    V    = tp(index, index);
    Eigen::VectorXd bk(b(index, k));
    F(k) = (bk.array() * (V.colPivHouseholderQr().solve(bk)).array()).sum() / df1;
  }
#endif
  return Rcpp::List::create(_["F"] = F, _["df1"] = df1, _["df2"] = df2); // removed _["ru"] = e as exoport
}

// Computes sqrt of matrices
Eigen::MatrixXd matrixSqrt(const Eigen::MatrixXd& A) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  Eigen::VectorXd sqrt_evals = es.eigenvalues().array().sqrt();
  return es.eigenvectors() * sqrt_evals.asDiagonal() * es.eigenvectors().transpose();
}


// This function computes KP stat
//[[Rcpp::export]]
Rcpp::List fKPstat(const Eigen::MatrixXd& endo_,
                   const Eigen::MatrixXd& X,
                   const Eigen::MatrixXd& Z_,
                   const Eigen::VectorXi& index,
                   const Eigen::ArrayXd& cumsn,
                   const int& HAC = 0) {
  Eigen::MatrixXd iXX((X.transpose() * X).inverse());
  Eigen::MatrixXd endo(endo_ - X * iXX * X.transpose() * endo_);
  Eigen::MatrixXd Z(Z_(Eigen::all, index) - X * iXX * X.transpose() * Z_(Eigen::all, index));
  // Eigen::MatrixXd endo(endo_);
  // Eigen::MatrixXd Z(Z_);
  int n(endo.rows()), nendo(endo.cols()), l(Z.cols()), ngroup(cumsn.size() - 1);
  Eigen::MatrixXd ZZ(Z.transpose() * Z);
  Eigen::MatrixXd iZZ(ZZ.inverse());
  Eigen::MatrixXd Zendo(Z.transpose() * endo);
  
  // estimator
  Eigen::MatrixXd Pi(Zendo.transpose() * iZZ);
  Eigen::VectorXd pi(Pi.reshaped(l * nendo, 1)); // Eigen::kroneckerProduct(Eigen::MatrixXd::Identity(l, l), Zendo.transpose()) * ZZ.inverse().reshaped(l*l, 1)
  
  // vec(Ze)
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(l * nendo, l * nendo));
  for (int s1(0); s1 < l; ++ s1) {
    for (int s2(0); s2 < nendo; ++ s2) {
      R(s1 * nendo + s2, s2 * l + s1) = 1;
    }
  }
  
  Eigen::MatrixXd eps(endo - Z * Pi.transpose());
  Eigen::MatrixXd vecZe(n, l*nendo);
  for (int s(0); s < nendo; ++ s) {
    vecZe.block(0, s*l, n, l) = (Z.array().colwise()*eps.col(s).array()).matrix();
  }
  
  // Variance of vec(Ze), covendo and covz
  Eigen::MatrixXd VvecZe(Eigen::MatrixXd::Zero(l*nendo, l*nendo)),
  Eee(Eigen::MatrixXd::Zero(nendo, nendo)),
  Ezz(Eigen::MatrixXd::Zero(l, l));
  if (HAC <= 2) {
    VvecZe = vecZe.transpose() * vecZe;
    Eee    = eps.transpose() * eps;
    Ezz    = Z.transpose() * Z;
  } else {
    for (int r(0); r < ngroup; ++ r) {
      int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
      Eigen::VectorXd tp(vecZe(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix());
      VvecZe += tp * tp.transpose();
      tp      = eps(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix();
      Eee    += tp * tp.transpose();
      tp      = Z(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix();
      Ezz    += tp * tp.transpose();
    }
  }
  
  // Variance of pi
  Eigen::MatrixXd H(R * Eigen::kroneckerProduct(Eigen::MatrixXd::Identity(nendo, nendo), iZZ));
  Eigen::MatrixXd varpi(H * VvecZe * H.transpose()); // O(1/n)
  
  // normalisation
  Eigen::LLT<Eigen::MatrixXd> tpF(ZZ * Ezz.colPivHouseholderQr().solve(ZZ)), tpG(Eee.inverse());
  Eigen::MatrixXd F(tpF.matrixL()); // O(sqrt(n))
  Eigen::MatrixXd G(tpG.matrixL().transpose()); // O(1/sqrt(n))
  F *= sqrt(ngroup); // O(n)
  
  // Theta and its variance
  Eigen::MatrixXd Theta(G * Pi * F.transpose());
  Eigen::VectorXd theta(Theta.reshaped(l*nendo, 1));
  Eigen::MatrixXd FG(Eigen::kroneckerProduct(F, G));
  Eigen::MatrixXd vartheta(FG * varpi * FG.transpose());
  // cout << vartheta << endl;
  
  // SDV decomposition of Theta
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Theta, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::MatrixXd U = svd.matrixU(); //nendo * nendo
  Eigen::VectorXd d = svd.singularValues();
  Eigen::MatrixXd ddiag = d.asDiagonal();
  Eigen::MatrixXd D(nendo, l);
  D << ddiag, Eigen::MatrixXd::Zero(nendo, l - nendo); //l*nendo
  Eigen::MatrixXd V = svd.matrixV(); //nendo * nendo
  
  //U12, U22, V12, V22
  int q(nendo - 1);
  Eigen::MatrixXd U12(U.block(0, q, q, nendo - q));
  Eigen::MatrixXd U22(U.block(q, q, nendo - q, nendo - q));
  Eigen::MatrixXd V12(V.block(0, q, q, l - q));
  Eigen::MatrixXd V22(V.block(q, q, l - q, l - q));
  
  // Aqper and Bqper
  Eigen::MatrixXd U12U22(nendo, nendo - q), V12V22(l, l - q);
  U12U22 << U12, U22;
  V12V22 << V12, V22;
  
  Eigen::MatrixXd Aper(U12U22 * U22.colPivHouseholderQr().solve(matrixSqrt(U22 * U22.transpose())));
  Eigen::MatrixXd Bper((V12V22 * V22.colPivHouseholderQr().solve(matrixSqrt(V22 * V22.transpose()))).transpose());
  
  // lambda and its varianve
  Eigen::MatrixXd BAper(Eigen::kroneckerProduct(Bper, Aper.transpose()));
  Eigen::VectorXd lambda (BAper * theta);
  Eigen::MatrixXd varlambda (BAper * vartheta * BAper.transpose());
  
  // statistic
  double stat((lambda.transpose() * varlambda.colPivHouseholderQr().solve(lambda))(0, 0));
  return Rcpp::List::create(_["stat"] = stat, _["df"] = (nendo - q)*(l - q));
}
