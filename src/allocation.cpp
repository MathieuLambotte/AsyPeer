// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <random>
#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace std;

////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Auxiliaries Function 
// Choose combination
// [[Rcpp::export]]
std::vector<std::vector<int>> combn(const std::vector<int>& elements, 
                                    unsigned int n,
                                    unsigned int k) {
  std::vector<std::vector<int>> result;
  
  std::vector<bool> mask(n, false);
  std::fill(mask.begin(), mask.begin() + k, true);
  
  do {
    std::vector<int> current_combination;
    for (unsigned int i = 0; i < n; ++i) {
      if (mask[i]) {
        current_combination.push_back(elements[i]);
      }
    }
    result.push_back(current_combination);
  } while (std::prev_permutation(mask.begin(), mask.end()));
  
  return result;
}


// Best response function for a single network
//[[Rcpp::export]]
Eigen::ArrayXd BRSingle(const Eigen::ArrayXd& alpha,
                        const Eigen::ArrayXd& y,
                        const Eigen::MatrixXd& G,
                        const Eigen::ArrayXd& betadelta,
                        const unsigned int& n,
                        const std::vector<Eigen::ArrayXi>& idpeer,
                        const Eigen::ArrayXd& d){
  // int n(alpha.size());
  // parameter
  double bl(betadelta(0)), bh(betadelta(1)), delta(betadelta(2));
  
  // Compute ybar
  Eigen::ArrayXd ybar = G * y.matrix(); 
  
  // Compute new y
  Eigen::ArrayXd ynew(alpha);
  Eigen::MatrixXd Gt      = G.transpose();
  for (unsigned int i(0); i < n; ++ i) {
    if (idpeer[i].size() > 0) {
      
      // unique values of peer outcomes
      Eigen::ArrayXd ypeer = y(idpeer[i]);
      std::set<double> uypeer_vec(ypeer.data(), ypeer.data() + ypeer.size()); // unique and sorted
      Eigen::ArrayXd uypeer(uypeer_vec.size());
      std::copy(uypeer_vec.begin(), uypeer_vec.end(), uypeer.data());
      
      // Ai: -inf, yi(1), ..., yi(ni), +inf
      Eigen::ArrayXd Ai(uypeer.size() + 2); Ai << R_NegInf, uypeer, R_PosInf; 
      int ell(1); // position of ai in Ai
      bool cont(true); // says if ell should be incremented
      while (cont) { // continue and remaining some ai
        // compute the marginal utility at Ai(ell)
        Eigen::ArrayXd Ghi((y > Ai(ell)).select(Gt.col(i), 0));
        double yhi((Ghi*y).sum());
        double gih(Ghi.sum());
        double marg((alpha(i) + (delta + bl)*ybar(i) + (bh - bl)*yhi) - (1 + bl * d(i) + (bh - bl) * gih) * Ai(ell));
        if (marg >= 0) {
          ++ ell;
        } else {
          cont = false;
        }
      }
      // We know the upper bound ell + 1 (in the paper), which is ell here.
      Eigen::ArrayXd Ghi((y >= Ai(ell)).select(Gt.col(i), 0));
      double yhi((Ghi*y).sum());
      double gih(Ghi.sum());
      ynew(i) = (alpha(i) + (delta + bl)*ybar(i) + (bh - bl)*yhi) / (1 + bl * d(i) + (bh - bl) * gih);
    }
  }
  
  return ynew;
}

// Nash Equilibrium for a single network
//[[Rcpp::export]]
Eigen::VectorXd NashSingle(const Eigen::ArrayXd y0,
                           const Eigen::ArrayXd& alpha,
                           const Eigen::MatrixXd& G,
                           const Eigen::ArrayXd& betadelta,
                           const std::vector<Eigen::ArrayXi>& idpeer,
                           const Eigen::ArrayXd& d,
                           const double& tol){
  Eigen::ArrayXd y = y0;
  int n = y0.size();
  
  computeBR: Eigen::ArrayXd yst = BRSingle(alpha, y, G, betadelta, n, idpeer, d);
  
  // check convergence
  double dist = ((yst - y.array())/(y.array() + 1e-50)).abs().maxCoeff();
  y           = yst;
  if (dist > tol) goto computeBR;
  
  return y;
}


// This function compute alpha
//[[Rcpp::export]]
Eigen::ArrayXd falpha(const Eigen::ArrayXd& betadelta,
                      const Eigen::MatrixXd& G,
                      const Eigen::ArrayXd& y,
                      const Eigen::ArrayXi& isolates) {
  
  // Parameters
  unsigned int n = y.size();
  double betal(betadelta(0));
  double betah(betadelta(1));
  double delta(betadelta(2));
  double theta1 = (betal + delta) / (1.0 + betal);
  double theta2 = (betah - betal) / (1.0 + betal);
  
  // gh, and ybarhigh
  Eigen::ArrayXd gh(n), ybh(n);
  for (unsigned int i(0); i < n; ++ i) {
    Eigen::ArrayXd Ghi = (y > y(i)).select(G.row(i).transpose(), 0);
    gh(i)  = Ghi.sum();
    ybh(i) = (Ghi * y).sum();
  }
  
  // ybar and ycheck
  Eigen::ArrayXd ybar   = G * y.matrix();
  Eigen::ArrayXd ycheck = ybh - gh * y; 
  
  //alpha
  Eigen::ArrayXd alpha  = y - ybar * theta1 - ycheck * theta2;
  alpha *= (isolates == 0).select(1 + betal, Eigen::ArrayXd::Constant(n, 1.0));
  
  return alpha;
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////// Rank according to the symmetric model
// [[Rcpp::export]]
Rcpp::List fRankSym(const double& beta, // for the misspecified symmetric model
                    const double& delta, // for the misspecified symmetric model
                    const Eigen::ArrayXd& betadelta, // for the asymmetric model
                    const Eigen::ArrayXd& y,
                    const Eigen::ArrayXd& alpha,
                    const Eigen::MatrixXd& G,
                    const std::vector<Eigen::ArrayXi>& idpeer,
                    const Eigen::ArrayXd& treat,
                    const Eigen::ArrayXd& d,
                    const Eigen::ArrayXi& isolates,
                    const double& tol,
                    const unsigned int& nthread) {
  double theta = (beta + delta) / (1.0 + beta);
  unsigned int n = G.rows();
  
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n, n) - theta * G.transpose();
  Eigen::ArrayXd AinvOne = A.colPivHouseholderQr().solve(Eigen::VectorXd::Ones(n));
  AinvOne *= (isolates == 0).select(1.0/(1 + beta), Eigen::ArrayXd::Constant(n, 1.0));
  
  // Rank in AinvOne
  Eigen::ArrayXi RAinvOne = Eigen::ArrayXi::LinSpaced(n, 0, n - 1);
  std::sort(RAinvOne.data(), RAinvOne.data() + n,
            [&AinvOne](int a, int b) {
              return AinvOne(a) > AinvOne(b);
            });
  
  // For each i we define the set of individuals that should be treated together
  // This is important for ties
  std::vector<Eigen::ArrayXi> setRank(n);
  Eigen::ArrayXi NTreat(n);
  Eigen::ArrayXd sAinvOne = AinvOne(RAinvOne);
  unsigned int i = 0;
  while (i < n) {
    unsigned int j = 1;
    
    while ((i + j) < n &&
           std::abs(sAinvOne(i) - sAinvOne(i + j)) < tol) {
      ++j;
    }
    
    setRank[i] = RAinvOne.segment(i, j);
    for (unsigned int ii = (i + 1); ii < (i + j); ++ ii) {
      setRank[ii] = setRank[i];
    }
    
    NTreat.segment(i, j) = j;
    
    i += j;
  }
  
  // Compute spillover 
  Eigen::ArrayXd Spill(n);
  Eigen::ArrayXd Diffy(n);
  Eigen::ArrayXd streat(n);
  double ysum = y.sum();
  
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for(unsigned int k = 0; k < n; ++ k) {
    // treatment
    Eigen::ArrayXd alphak = alpha;
    double treatsum = 0;
    for (unsigned int l = 0; l <= k ; ++ l) {
      alphak(setRank[l]) += treat(setRank[l]) / NTreat(l);
      treatsum += treat(setRank[l]).sum() / NTreat(l);
    }
    // Nash equilibrium
    Diffy(k)  = NashSingle(y, alphak, G, betadelta, idpeer, d, tol).sum() - ysum;
    Spill(k)  = Diffy(k) / treatsum;
    streat(k) = treatsum;
  }
#else
  for(unsigned int k = 0; k < n; ++ k) {
    // treatment
    Eigen::ArrayXd alphak = alpha;
    double treatsum = 0;
    for (unsigned int l = 0; l <= k ; ++ l) {
      alphak(setRank[l]) += treat(setRank[l]) / NTreat(l);
      treatsum += treat(setRank[l]).sum() / NTreat(l);
    }
    // Nash equilibrium
    Diffy(k)  = NashSingle(y, alphak, G, betadelta, idpeer, d, tol).sum() - ysum;
    Spill(k)  = Diffy(k) / treatsum;
    streat(k) = treatsum;
  }
#endif
  
  return Rcpp::List::create(
    Rcpp::_["sum.treat"] = streat,
    Rcpp::_["setRank"]   = setRank,
    Rcpp::_["Spillover"] = Spill,
    Rcpp::_["Diff.sumy"] = Diffy
  );
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////// Rank according to the asymmetric model 
//////////////////////// but using the sequential algo
// [[Rcpp::export]]
Rcpp::List fRankASym(const Eigen::ArrayXd& y,
                     const Eigen::ArrayXd& alpha,
                     const Eigen::ArrayXd& betadelta,
                     const Eigen::MatrixXd& G,
                     const std::vector<Eigen::ArrayXi>& idpeer,
                     const Eigen::ArrayXd& treat,
                     const Eigen::ArrayXd& d,
                     const double& tol,
                     const unsigned int& nthread,
                     const unsigned int& seed,
                     const bool& print) {
  unsigned int n = G.rows();
  double ysum(y.sum());
  
  // random generator
  std::mt19937 rng(seed);
  
  // initial alpha and outcome for forward
  Eigen::ArrayXd alphaF = alpha;
  Eigen::ArrayXd yF(y);// Already compatible with alphaF
  
  // initial alpha for backward
  Eigen::ArrayXd alphaB = alpha + treat;
  Eigen::ArrayXd yB = NashSingle(y, alphaB, G, betadelta, idpeer, d, tol);
  
  // Rank for forward and backward
  Eigen::ArrayXXi Rk(n, 2);
  
  // Spillover for forward and backward
  Eigen::ArrayXXd Spill(n, 2), Diffy(n, 2), streat(n, 2);
  Diffy(n - 1, 1) = yB.sum() - ysum;
  
  // Remaining indices for forward and backward
  std::vector<int> RF(n), RB(n);
  std::iota(RF.begin(), RF.end(), 0);
  RB = RF;
  
  // Sequential algorithm
  Progress p(n, print);
  for (unsigned int k = 0; k < n; ++ k) {
    if (Progress::check_abort())
      Rcpp::stop("interrupted");
    
    // Save and y
    Eigen::ArrayXXd simyF(n, n - k);
    Eigen::ArrayXXd simyB(n, n - k);
    
#ifdef _OPENMP
    omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < (n - k); ++ i) {
      
      // Forward 
      Eigen::ArrayXd alphak = alphaF;
      alphak(RF[i]) += treat(RF[i]);
      simyF.col(i) = NashSingle(yF, alphak, G, betadelta, idpeer, d, tol);
      
      // Backward
      alphak = alphaB;
      alphak(RB[i]) -= treat(RB[i]);
      simyB.col(i) = NashSingle(yB, alphak, G, betadelta, idpeer, d, tol);
    }
#else
    for (unsigned int i = 0; i < (n - k); ++ i) {
      
      // Forward 
      Eigen::ArrayXd alphak = alphaF;
      alphak(RF[i]) += treat(RF[i]);
      simyF.col(i) = NashSingle(yF, alphak, G, betadelta, idpeer, d, tol);
      
      // Backward
      alphak = alphaB;
      alphak(RB[i]) -= treat(RB[i]);
      simyB.col(i) = NashSingle(yB, alphak, G, betadelta, idpeer, d, tol);
    }
#endif
    
    // Find  optimal index to target (forward) or remove (backword)
    Eigen::Index iF, iB; 
    // Forward
    {
      Eigen::ArrayXd SpilF = simyF.colwise().sum() - ysum;
      double best = SpilF.maxCoeff();
      Diffy(k, 0) = best;
      
      // Randomize rank in case of ties
      std::vector<int> ties;
      for (unsigned int i = 0; i < (n - k); ++ i) {
        if (std::abs(best - SpilF(i)) < tol) {
          ties.push_back(i);
        }
      }
      std::shuffle(ties.begin(), ties.end(), rng);
      iF = ties[0];
      
    }
    // Backward
    {
      Eigen::ArrayXd SpilB = simyB.colwise().sum() - ysum;
      double best = SpilB.maxCoeff();
      if (k < n - 1) {
        Diffy(n - 2 - k, 1) = best;
      }
      
      // Randomize in case of ties
      std::vector<int> ties;
      for (unsigned int i = 0; i < (n - k); ++ i) {
        if (std::abs(best - SpilB(i)) < tol) {
          ties.push_back(i);
        }
      }
      std::shuffle(ties.begin(), ties.end(), rng);
      iB = ties[0];
    }
    
    // update alphaF alphaB
    alphaF(RF[iF]) += treat(RF[iF]);
    alphaB(RB[iB]) -= treat(RB[iB]);
    
    // Update y
    yF = simyF.col(iF);
    yB = simyB.col(iB);
    
    // Save ranks
    Rk(k, 0) = RF[iF];
    Rk(n - k - 1, 1) = RB[iB];
    
    // Remove selected individual from RF and RB
    RF.erase(RF.begin() + iF);
    RB.erase(RB.begin() + iB);
    
    // update progress bar
    p.increment();
  }
  
  // Spillover
  for (unsigned int k = 0; k < n; ++ k) {
    streat(k, 0) = treat(Rk.col(0).head(k + 1)).sum();
    Spill(k, 0)  = Diffy(k, 0) / streat(k, 0);
    streat(k, 1) = treat(Rk.col(1).head(k + 1)).sum();
    Spill(k, 1)  = Diffy(k, 1) / treat(Rk.col(0).head(k + 1)).sum();
  }
  
  return Rcpp::List::create(Rcpp::_["sum.treat"] = streat,
                            Rcpp::_["Rank"]      = Rk,
                            Rcpp::_["Spillover"] = Spill,
                            Rcpp::_["Diff.sumy"] = Diffy);
}

