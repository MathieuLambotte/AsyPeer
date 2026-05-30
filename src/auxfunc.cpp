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
                             const Eigen::MatrixXd& V,
                             const Eigen::ArrayXi& cumsn,
                             const Eigen::ArrayXi& nvec,
                             const int& power,
                             const int& nthread) {
  int kV(V.cols()), n(nvec.sum()), ngroup(nvec.size());
  Eigen::MatrixXd out(n, kV * (power + 1));
  out.block(0, 0, n, kV) = V;
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int m = 0; m < ngroup; ++m) {
    for (int k = 1; k <= power; ++k) {
      out.block(cumsn(m), kV * k, nvec(m), kV).noalias() = 
        G[m] * out.block(cumsn(m), kV * (k - 1), nvec(m), kV);
    }
  }
#else
  for (int m = 0; m < ngroup; ++m) {
    for (int k = 1; k <= power; ++k) {
      out.block(cumsn(m), kV * k, nvec(m), kV).noalias() = 
        G[m] * out.block(cumsn(m), kV * (k - 1), nvec(m), kV);
    }
  }
#endif
  return out.rightCols(kV * power);
}

// This function removes columns to obtain full rank matrices
// Taken from the QuantilePeer package
//[[Rcpp::export]]
Eigen::Array<bool, Eigen::Dynamic, 1> fcheckrankEigen(const Eigen::MatrixXd& X, const double& tol = 1e-10) {
  int n(X.rows());
  // Eigen::RowVectorXd m(X.colwise().mean());
  // Eigen::RowVectorXd s(((X.rowwise() - m).array().square().colwise().sum() / n).sqrt());
  // m = (s.array() < tol).select(0, m);
  // s = (s.array() < tol).select(1, s);
  // Eigen::MatrixXd U((X.rowwise() - m).array().rowwise() / s.array());
  Eigen::MatrixXd U = X.transpose()*X/n;
  // std::cout<<U<<std::endl;
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(U);
  Eigen::MatrixXd R(qr.matrixQR().topRows(U.cols()));
  // std::cout<<R.diagonal().transpose()<<std::endl;
  return R.diagonal().array().abs() > tol;
}

// Assigning folds to groups
//[[Rcpp::export]]
Eigen::ArrayXi fassignfold(const Eigen::ArrayXi& group,
                           const int& nfold) {
  std::unordered_set<int> unigroup(group.data(), group.data() + group.size());
  int ngroup(unigroup.size());
  if (ngroup < 2) {
    Rcpp::stop("Only one subnet remains for the intensive/extensive model.");
  }
  
  // group should take values from 0 to ngroup - 1 with possible duplication
  if ((group.minCoeff() != 0) || (group.maxCoeff() != (ngroup - 1))) {
    Rcpp::stop("group should take values from 0, 1, 2, ... without jumps.");
  }
  
  // Number of pairs per group
  Eigen::ArrayXi nvec(Eigen::ArrayXi::Zero(ngroup));
  for (int i(0); i < group.size(); ++i) {
    nvec(group(i)) += 1;
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
  return fold(group);
}

//[[Rcpp::export]]
Eigen::ArrayXXd Demean_separate(const Eigen::ArrayXXd& X,
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
Eigen::ArrayXXd Demean_common(const Eigen::ArrayXXd& X,
                              const Eigen::ArrayXi& cumsn,
                              const int& nthread){
  int ngroup(cumsn.size() - 1);
  Eigen::ArrayXXd out(X);
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int s = 0; s < ngroup; ++ s) {
    int n1(cumsn(s)), ns(cumsn(s + 1) - cumsn(s));
    out(Eigen::seq(n1, ns), Eigen::all).rowwise() -= out(Eigen::seq(n1, ns), Eigen::all).colwise().mean();
  }
#else
  for (int s = 0; s < ngroup; ++ s) {
    int n1(cumsn(s)), ns(cumsn(s + 1) - cumsn(s));
    out(Eigen::seq(n1, ns), Eigen::all).rowwise() -= out(Eigen::seq(n1, ns), Eigen::all).colwise().mean();
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


