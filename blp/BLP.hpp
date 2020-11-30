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
  BLP(const unsigned num_periods, const unsigned num_bins_renda, const unsigned\
      num_bins_idade, const std::vector<std::vector<unsigned>> areas, unsigned\
      ns_, std::vector<double> theta2_, double contract_tol_, const unsigned\
      max_threads=64);
  ~BLP()
  {
    S.clear();
    delta.clear();
    X1.clear();
    X2.clear();
    Z.clear();
    mkt_id.clear();
    area_id.clear();
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
  unsigned ns;
  unsigned num_threads;
  double contract_tol;
  void allocate();
  void calc_shares();
  void contraction(bool increase_tol=true);
  void calc_phi();
  void calc_theta1();
  void gmm();
  
private:
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
  // calc objs
  std::vector<ublas::vector<double>> s_aux;
  std::vector<ublas::vector<double>> s_calc;
  ublas::vector<double> exp_delta1;
  ublas::vector<double> exp_delta2;
  ublas::matrix<double> phi;
  ublas::matrix<double> phi_inv;
  ublas::vector<double> omega;
};

#endif
