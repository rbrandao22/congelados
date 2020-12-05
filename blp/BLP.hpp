#ifndef BLPHEADERDEF 
#define BLPHEADERDEF

#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;


class BLP
{

public:
  BLP(const std::string persist_file2_, const unsigned num_periods, const\
      unsigned num_bins_renda, const unsigned num_bins_idade, const\
      std::vector<std::vector<unsigned>> areas, unsigned ns_,\
      std::vector<double> theta2_, double contract_tol_, const double num_lim_,\
      const unsigned max_threads=64);
  ~BLP()
  {
    S.clear();
    delta.clear();
    X1.clear();
    X2.clear();
    Z.clear();
    mkt_id.clear();
    area_id.clear();
    theta1.clear();
    theta2.clear();
    exp_delta1.clear();
    exp_delta2.clear();
    phi_inv.clear();
    omega.clear();
    s_ijt.clear();
    Ddelta.clear();
    Ddelta1.clear();
    Ddelta2.clear();
    Ddelta1_aux.clear();
    Ddelta2a_aux.clear();
    Ddelta2b_aux.clear();
    grad.clear();
    grad_aux.clear();
    grad_adj.clear();
    obj_aux.clear();
    obj_aux2.clear();
  }
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & S;
    ar & delta;
    ar & X1;
    ar & X2;
    ar & Z;
    ar & mkt_id;
    ar & area_id;
  }
  void allocate();
  void calc_shares();
  void contraction();
  void calc_phi_inv();
  void calc_theta1();
  void calc_Ddelta();
  void perturb_params();
  void grad_calc();
  void objective_calc();
  void gmm(double nr_tol, double step_size, const unsigned max_iter);
  void persist();
  
private:
  std::string persist_file2;
  unsigned ns;
  unsigned num_threads;
  double contract_tol;
  double num_lim;
  bool ctol_inc;  // following Nevo, contraction tol may initially increase
  double grad_norm;
  double obj_value;
  bool perturb;
  ublas::vector<double> S;
  ublas::vector<double> delta;
  ublas::matrix<double> X1;
  ublas::matrix<double> X2;
  ublas::matrix<double> Z;
  ublas::vector<unsigned> mkt_id;
  ublas::vector<unsigned> area_id;
  // random coeff draws
  boost::multi_array<double, 3> v;
  boost::multi_array<double, 3> D;
  // params
  ublas::vector<double> theta1;
  ublas::vector<double> theta2;
  /// calc objs
  std::vector<ublas::vector<double>> s_aux;
  std::vector<ublas::vector<double>> s_calc;
  ublas::vector<double> exp_delta1;
  ublas::vector<double> exp_delta2;
  ublas::matrix<double> phi_inv;
  ublas::vector<double> omega;
  // grad specific
  ublas::matrix<double> s_ijt;
  ublas::matrix<double> Ddelta;
  std::vector<ublas::matrix<double>> Ddelta1;
  std::vector<ublas::matrix<double>> Ddelta2;
  std::vector<ublas::matrix<double>> Ddelta1_aux;
  std::vector<ublas::matrix<double>> Ddelta2a_aux;
  std::vector<ublas::matrix<double>> Ddelta2b_aux;
  ublas::vector<double> grad;
  ublas::matrix<double> grad_aux;
  ublas::vector<double> grad_adj;
  ublas::vector<double> obj_aux;
  ublas::matrix<double> obj_aux2;
};

#endif
