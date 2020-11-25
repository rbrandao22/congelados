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
  GenArrays(const unsigned num_mkts, const unsigned num_draws, const\
	    std::vector<std::vector<unsigned>> areas, const unsigned\
	    num_bins_renda, const unsigned num_bins_idade);

private:
  ublas::matrix<double> S;
  ublas::matrix<double> X_p;
  boost::multi_array<double, 3> v;
  boost::multi_array<double, 3> D;
};

#endif
