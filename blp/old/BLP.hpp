#ifndef BLPHEADERDEF 
#define BLPHEADERDEF

#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;


class BLP
{

public:
  BLP(const std::string initguess_f, double min_share_, const double\
      contract_tol_, const unsigned max_iter_contract_, const double\
      penalty_param1_, const unsigned penalty_param2_, const double\
      init_tetra_size1, const double init_tetra_size2, unsigned const\
      max_threads=64);

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & s_obs_wg;
    ar & mkt_id;
    ar & pop_ave;
    ar & X;
    ar & Z;
  }
  unsigned params_nbr;
  unsigned num_threads;
  bool do_NelderMead;
  bool halt_check;
  void allocate();
  void calc_objective(unsigned pt, bool add_penalty=false);
  void updatePs(const double inc);
  void grad_calc(const double inc, const double tol);
  void step(const double step_size, const double max_step, const double\
	    step_factor, unsigned iter_nbr);
  void updatePs_NM();
  void nelder_mead(const double alpha, const double beta, const double gamma);
  void updateP0_NM(unsigned iter_nbr);
  void variance();
  void persist(const std::string persist_file2);
  void persist_ig(const std::string initguess_f);
  
private:
  // Params
  double min_share;
  double contract_tol;
  unsigned max_iter_contract;
  double penalty_param1;
  unsigned penalty_param2;
  double init_tetra_size1;
  double init_tetra_size2;
  // Exogenous vars
  unsigned N;
  ublas::vector<double> s_obs_wg;
  ublas::vector<double> mkt_id;
  ublas::vector<double> pop_ave;
  ublas::matrix<double> X;
  /* X matrix structure: 0 constant
                         1 fare (bin center R$/1000)
                         2 distance (km / 1000)
                         3 distance^2
                         4 dummy TAM
                         5 dummy GLO
                         6 dummy ONE
                         7 dummy other (base AZU)
                         8 time dummy second period
                         9 time dummy third period */
  ublas::matrix<double> Z;
  /* Z matrix structure: 0 constant
                         1 aviation fuel price (querosene, instrument)
                         2 distance (km / 1000)
                         3 distance^2
                         4 dummy TAM
                         5 dummy GLO
                         6 dummy ONE
                         7 dummy other (base AZU)
                         8 time dummy second period
                         9 time dummy third period */

  std::vector<ublas::vector<double>> P;
  /* Params are 23 total:
     0: beta1_0, 1: beta1_1, 2: beta1_2, 3: beta1_3, 4: beta1_4, 5: beta1_5,
     6: beta1_6, 7: beta1_7, 8: beta1_8, 9: beta1_9;
     10: beta2_0, 11: beta2_1, 12: beta2_2, 13: beta2_3, 14: beta2_4,
     15: beta2_5; 16: beta2_6, 17: beta2_7, 18: beta2_8, 19: beta2_9;
     20: gamma, 21: lambda, 22: mu 
     NM procedures takes number of params + 1 = 24 P's, + 
       P^bar (params_nbr+1), P* (params_nbr+2) and P** (params_nbr+3) */

  // Objective function values
  std::vector<double> y;
  
  /// Calc vars
  // unobserved utility
  std::vector<ublas::vector<double>> xi0;
  std::vector<ublas::vector<double>> xi1;
  // other vars
  std::vector<ublas::vector<double>> s_aux1;
  std::vector<ublas::vector<double>> s_aux2;
  std::vector<ublas::vector<double>> D1;
  std::vector<ublas::vector<double>> D2;
  std::vector<ublas::vector<double>> s_calc;
  std::vector<ublas::vector<double>> ln_s_obs;
  /// NR specific
  ublas::vector<double> grad;
};

#endif
