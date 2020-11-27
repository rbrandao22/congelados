#ifndef GENARRAYSHEADERDEF
#define GENARRAYSHEADERDEF

#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;


class GenArrays
{
public:
  GenArrays(const unsigned num_periods, const unsigned num_draws, const\
	    std::vector<std::vector<unsigned>> areas, const unsigned\
	    num_bins_renda, const unsigned num_bins_idade);
  ~GenArrays()
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
  }

private:
  double num_mkts;
  double num_prods;
  ublas::vector<double> S;
  ublas::vector<double> delta;
  ublas::matrix<double> X1;
  ublas::matrix<double> X2;
  ublas::matrix<double> Z;
  ublas::vector<unsigned> mkt_id;
  boost::multi_array<double, 3> v;
  boost::multi_array<double, 3> D;
};

#endif
