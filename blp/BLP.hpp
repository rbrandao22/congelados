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
      num_bins_idade, unsigned ns, const std::vector<std::vector<unsigned>>\
      areas, const unsigned max_threads=64);
  ~BLP()
  {
    S.clear();
    delta.clear();
    X1.clear();
    X2.clear();
    Z.clear();
    mkt_id.clear();
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
  void allocate();
  void calc_objective(std::vector<double> theta2_);
  
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
  std::vector<double> theta2;
  // calc objs
  std::vector<ublas::vector<double>> s_aux;
  std::vector<ublas::vector<double>> s_calc;
};

#endif
