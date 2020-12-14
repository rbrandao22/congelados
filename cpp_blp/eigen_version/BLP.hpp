#ifndef BLPHEADERDEF 
#define BLPHEADERDEF

#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <Eigen/Dense>

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
    mkt_id.clear();
    area_id.clear();
  }
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & S_;
    ar & delta_;
    ar & X1_;
    ar & X2_;
    ar & Z_;
    ar & mkt_id;
    ar & area_id;
  }
  unsigned ns;
  unsigned jt_size;
  unsigned num_threads;
  double contract_tol;
  void transf_eigen();
  void allocate();
  void calc_shares();
  void contraction(bool increase_tol=true);
  void calc_phis();
  void calc_theta1();
  void calc_Ddelta();
  void gmm();
  
private:
  ublas::vector<double> S_;
  ublas::vector<double> delta_;
  ublas::matrix<double> X1_;
  ublas::matrix<double> X2_;
  ublas::matrix<double> Z_;
  ublas::vector<unsigned> mkt_id;
  ublas::vector<unsigned> area_id;
  // random coeff draws
  boost::multi_array<double, 3> v;
  boost::multi_array<double, 3> D;
  // Eigen transfs
  Eigen::VectorXd S;
  Eigen::VectorXd delta;
  Eigen::MatrixXd X1;
  Eigen::MatrixXd X2;
  Eigen::MatrixXd Z;
  // params
  Eigen::VectorXd theta1;
  Eigen::VectorXd theta2;
  /// calc objs
  std::vector<Eigen::VectorXd> s_aux;
  std::vector<Eigen::VectorXd> s_calc;
  Eigen::VectorXd exp_delta1;
  Eigen::VectorXd exp_delta2;
  Eigen::MatrixXd phi;
  Eigen::MatrixXd phi_inv;
  Eigen::VectorXd omega;
  // grad specific
  Eigen::MatrixXd s_ijt;
  Eigen::MatrixXd Ddelta;
  std::vector<Eigen::MatrixXd> Ddelta1;
  std::vector<Eigen::MatrixXd> Ddelta2;
  std::vector<Eigen::MatrixXd> Ddelta1_aux;
  std::vector<Eigen::MatrixXd> Ddelta2a_aux;
  std::vector<Eigen::MatrixXd> Ddelta2b_aux;
};

#endif
