// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#if defined(_OPENMP)
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

// Create data Li, Lj for random forest
//[[Rcpp::export]]
Rcpp::List fdataML(const Eigen::ArrayXd& y,
                   const Eigen::ArrayXXd& X,
                   const Eigen::ArrayXi& group, //indicate the group for each i
                   const Eigen::ArrayXi& IDi, // indicate the ID for each i in their group: 0, 1, ..
                   const std::vector<Eigen::ArrayXd>& gij, // gij for each i
                   const std::vector<Eigen::ArrayXi>& idpeer, //ID for each friends in each group
                   const Eigen::ArrayXi& ddni, // Number of friends for each i
                   const Eigen::ArrayXi& ddncs,// cumsum of ddni
                   const Eigen::ArrayXi& ncs,// cumsum of nvec
                   const int& nthread) {
  int n(group.size()), ddn(ddni.sum()), kX(X.cols());
  Eigen::ArrayXXd ddy(ddn, 7); //columns: group, IDi, IDj, gij, yi, yj, Indicator
  Eigen::ArrayXXd ddXi(ddn, kX), ddXj(ddn, kX);
  
#if defined(_OPENMP)
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) {
    if (ddni(i) > 0) {
      ddy.block(ddncs(i), 0, ddni(i), 1).setConstant(group(i));
      ddy.block(ddncs(i), 1, ddni(i), 1).setConstant(IDi(i));
      ddy.block(ddncs(i), 2, ddni(i), 1) = idpeer[i].cast<double>();
      ddy.block(ddncs(i), 3, ddni(i), 1) = gij[i];
      ddy.block(ddncs(i), 4, ddni(i), 1).setConstant(y(i));
      ddy.block(ddncs(i), 5, ddni(i), 1) = y(ncs(group(i)) + idpeer[i]);
      
      ddXi(Eigen::seqN(ddncs(i), ddni(i)), Eigen::all).rowwise() = 
        X.row(i);
      
      ddXj(Eigen::seqN(ddncs(i), ddni(i)), Eigen::all) = 
        X(ncs(group(i)) + idpeer[i], Eigen::all);
    }
  }
#else
  for (int i = 0; i < n; ++i) {
    if (ddni(i) > 0) {
      ddy.block(ddncs(i), 0, ddni(i), 1).setConstant(group(i));
      ddy.block(ddncs(i), 1, ddni(i), 1).setConstant(IDi(i));
      ddy.block(ddncs(i), 2, ddni(i), 1) = idpeer[i].cast<double>();
      ddy.block(ddncs(i), 3, ddni(i), 1) = gij[i];
      ddy.block(ddncs(i), 4, ddni(i), 1).setConstant(y(i));
      ddy.block(ddncs(i), 5, ddni(i), 1) = y(ncs(group(i)) + idpeer[i]);
      
      ddXi(Eigen::seqN(ddncs(i), ddni(i)), Eigen::all).rowwise() = 
        X.row(i);
      
      ddXj(Eigen::seqN(ddncs(i), ddni(i)), Eigen::all) = 
        X(ncs(group(i)) + idpeer[i], Eigen::all);
    }
  }
#endif
  
  ddy.col(6) = (ddy.col(5) > ddy.col(4)).cast<double>();
  return Rcpp::List::create(Rcpp::_["ddy"]  = ddy,
                            Rcpp::_["ddXi"] = ddXi,
                            Rcpp::_["ddXj"] = ddXj);
}

// This transform dyadic prediction to the instrument for check y
//[[Rcpp::export]]
Eigen::ArrayXXd fInstChecky(const Eigen::ArrayXXd& rhoddX,
                            const Eigen::ArrayXi& ddni,
                            const int& nthread) {
  int n(ddni.size()), K(rhoddX.cols());
  Eigen::ArrayXi ddncs(n + 1); ddncs(0) = 0;
  for (int i(0); i < n; ++i) {
    ddncs(i + 1) = ddncs(i) + ddni(i);
  }
  Eigen::ArrayXXd out(Eigen::ArrayXXd::Zero(n, K));
  
#if defined(_OPENMP)
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) {
    if (ddni(i) > 0) {
      out.row(i) = rhoddX.block(ddncs(i), 0, ddni(i), K - 1).colwise().sum();
    }
  }
#else
  for (int i = 0; i < n; ++i) {
    if (ddni(i) > 0) {
      out.row(i) = rhoddX.block(ddncs(i), 0, ddni(i), K - 1).colwise().sum();
    }
  }
#endif
  return out;
}
