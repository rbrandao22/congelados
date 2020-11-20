#ifndef GENARRAYSHEADERDEF
#define GENARRAYSHEADERDEF

#include <string>
#include <valarray>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

class GenArrays
{
public:
    GenArrays(const std::vector<std::string> dates_, const std::valarray\
            <double> bins_, const double pop_thres);
    ~GenArrays()
    {
        s_obs_wg.clear();
        mkt_id.clear();
        pop_ave.clear();
        X.clear();
	Z.clear();
    }
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & s_obs_wg;
        ar & mkt_id;
        ar & pop_ave;
        ar & X;
        ar & Z;
    }
    void gen_arrays();
    void gen_instruments();
private:
    std::vector<std::string> dates;
    /* products table structure: 0 string origem
                                 1 string destino
                                 2 string empresa
                                 3 double tarifa_lower (bin)
                                 4 double tarifa_upper (bin)
                                 5 unsigned empresa_nbr
                                 6 unsigned period_nbr */
    std::vector<std::tuple<std::string, std::string, std::string, double, \
            double, unsigned, unsigned>> products;
    ublas::vector<double> s_obs_wg;
    ublas::vector<double> mkt_id;
    ublas::vector<double> pop_ave;
    /* X matrix structure: 0 constant
                           1 fare (bin center)
                           2 distance
                           3 distance^2
                           4 carrier dummy
                           5 time dummy
    */
    ublas::matrix<double> X;
    ublas::vector<double> instruments;
    /* Z matrix structure: 0 constant
                           1 aviation fuel price (instrument, querosene)
                           2 distance
                           3 distance^2
                           4 carrier dummy
                           5 time dummy
    */
    ublas::matrix<double> Z;
  
};

#endif
