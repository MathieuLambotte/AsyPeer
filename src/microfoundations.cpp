/* GENERAL NOTATIONS
 * y       : is the vector of the outcome values
 * X       : is the matrix of the explanatory variables. Add the intercept 
 *           if it is included in the model. The intercept will not be added automatically. 
 * G       : is the network matrix List. That is G[r] is the subnetwork of the group r. 
 *           Gs(i,j) = measures the intensity of the outgoing link from i to j. 
 * ngroup  : is the number of groups.
 * theta   : is the vector of parameters ordered as follow: lambda_tau and explanatory variables (intercept is included)
 * n       : The sample size.
 * nvec    : Vector of agents in each subnet
 * cumsn   : is the cumulative sum for nvec, where the 1st element is 0 and the cumulative sum starts from the 2nd element.
 
 * tol     : A tolerance value for the iterative method solving the game. 
 * maxit   : The maximum number of iterations of the iterative method. If this
 *           number is reached, the algorithm stops and the last y is used as the solution. maxit
 *           is important for numerical reasons if tol is too small. 
 * peffects: is a vector containing (beta-, beta+, delta)
 */

// [[Rcpp::depends(RcppEigen)]]
// #include <RcppArmadillo.h>
#include <RcppEigen.h>
// #define NDEBUG
// #include <RcppNumerical.h>
// #include <RcppEigen.h>

// typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
// typedef Eigen::Map<Eigen::VectorXd> MapVect;

using namespace Rcpp;
using namespace Eigen;
using namespace std;

// Sort unique values
Eigen::ArrayXd sort_unique(const Eigen::ArrayXd& v) {
  std::set<double> s(v.data(), v.data() + v.size()); // unique and sorted
  Eigen::ArrayXd out(s.size());
  int i = 0;
  for (double val : s) {
    out(i++) = val;
  }
  return out;
}

// Compute average of friends values for some vector or vector u
Eigen::VectorXd peeravg(const Eigen::VectorXd& u,
                        const std::vector<Eigen::MatrixXd>& G,
                        const Eigen::ArrayXi& cumsn,
                        const Eigen::ArrayXi& nvec,
                        const int& ngroup,
                        const unsigned int& nthread) {
  Eigen::VectorXd out(u.size());
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#endif
#pragma omp parallel for schedule(static)
  for (int m = 0; m < ngroup; ++ m) {
    out.segment(cumsn(m), nvec(m)) =  G[m] * u.segment(cumsn(m), nvec(m));
  }
  return out;
}


// Best response function
Eigen::ArrayXd BR(const Eigen::ArrayXd& alpha,
                  const Eigen::ArrayXd& y,
                  const std::vector<Eigen::MatrixXd>& G,
                  const Eigen::ArrayXd& peffects,
                  const Eigen::ArrayXi& cumsn,
                  const Eigen::ArrayXi& nvec,
                  const std::vector<Eigen::ArrayXi>& idpeer,
                  const Eigen::ArrayXi& d,
                  const int& ngroup,
                  const unsigned int& nthread){
  int n(alpha.size());
  // parameter
  double bl(peffects(0)), bh(peffects(1)), delta(peffects(2));
  
  // Compute ybar
  Eigen::ArrayXd ybar = peeravg(y, G, cumsn, nvec, ngroup, nthread);
  // Compute new y
  Eigen::ArrayXd ynew(alpha);
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#endif
#pragma omp parallel for schedule(static)
  for (int m = 0; m < ngroup; ++ m) {
    int l(cumsn(m)); // Where the group starts in y.
    Eigen::ArrayXd ym        = y.segment(l, nvec(m));
    Eigen::MatrixXd Gmt      = G[m].transpose();
    for (int i(0); i < nvec(m); ++ i) {
      if (d(l) > 0) {
        Eigen::ArrayXd uypeer(sort_unique(ym(idpeer[l]))); // unique values of peer outcomes
        Eigen::ArrayXd Ai(uypeer.size() + 2); Ai << R_NegInf, uypeer, R_PosInf; // Ai: -inf, yi(1), ..., yi(ni), +inf
        int ell(1); // position of ai in Ai
        bool cont(true); // says if ell should be incremented
        while (cont) { // continue and remaining some ai
          // compute the marginal utility at Ai(ell)
          Eigen::ArrayXd Gmhi((ym > Ai(ell)).select(Gmt.col(i), 0));
          double yhi((Gmhi*ym).sum());
          double gih(Gmhi.sum());
          double marg((alpha(l) + (delta + bl)*ybar(l) + (bh - bl)*yhi) - (1 + bl * d(l) + (bh - bl) * gih) * Ai(ell));
          if (marg >= 0) {
            ++ ell;
          } else {
            cont = false;
          }
        }
        // We know the upper bound ell + 1 (in the paper), which is ell here.
        Eigen::ArrayXd Gmhi((ym >= Ai(ell)).select(Gmt.col(i), 0));
        double yhi((Gmhi*ym).sum());
        double gih(Gmhi.sum());
        ynew(l) = (alpha(l) + (delta + bl)*ybar(l) + (bh - bl)*yhi) / (1 + bl * d(l) + (bh - bl) * gih);
      }
      ++ l;
    }
  }
  return ynew;
}

// Nash Equilibrium
// y is initial solution
//[[Rcpp::export]]
int fNashE(Eigen::Map<Eigen::VectorXd> y,
           const Eigen::ArrayXd& alpha,
           const std::vector<Eigen::MatrixXd>& G,
           const Eigen::ArrayXd& peffects,
           const Eigen::ArrayXi& cumsn,
           const Eigen::ArrayXi& nvec,
           const std::vector<Eigen::ArrayXi>& idpeer,
           const Eigen::ArrayXi& d,
           const int& ngroup,
           const double& tol,
           const int& maxit,
           const unsigned int& nthread){
  int t(0);
  computeBR: ++t; // Best response dynamics
  
  // New y
  // cout<<y.transpose()<<endl;
  // cout<<"*****"<<endl;
  Eigen::ArrayXd yst = BR(alpha, y, G, peffects, cumsn, nvec, 
                          idpeer, d, ngroup, nthread);
  y.array().maxCoeff();
  
  // check convergence
  double dist = ((yst - y.array())/(y.array() + 1e-50)).abs().maxCoeff();
  y           = yst;
  cout<<"Iteration: "<<t<<" Distance: "<<dist<<endl;
  if (dist > tol && t < maxit) goto computeBR;
  return t;
}



