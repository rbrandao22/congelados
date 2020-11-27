#include <cmath>
#include <pqxx/pqxx>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "GenArrays.hpp"

namespace ublas = boost::numeric::ublas;


GenArrays::GenArrays(const unsigned num_periods, const\
		     std::vector<std::vector<unsigned>> areas)
{
  pqxx::connection C("dbname = congelados user = postgres password = passwd"\
          " hostaddr = 127.0.0.1 port = 5432");
  if (C.is_open()) {
      std::cout << "Opened database successfully: " << C.dbname() << \
          std::endl;
  } else {
      std::cout << "Can't open database" << std::endl;
      throw std::runtime_error("aborting");
  }
  pqxx::nontransaction N(C);

  /// fill observed shares - S, and mkt_id
  std::string lanches_shares_query = "SELECT * FROM lanches_shares;";
  pqxx::result R_ls(N.exec(lanches_shares_query));
  num_mkts = num_periods * areas.size();
  num_prods = R_ls.size();
  S.resize(num_prods * num_mkts);
  mkt_id.resize(num_prods * num_mkts);
  unsigned i = 0;
  unsigned mkt_counter = 0;
  for (unsigned t = 0; t != num_mkts; ++t) {
    for (auto c = R_ls.begin(); c != R_ls.end(); ++c) {
      if (!c[t+1].is_null()) {
	S(i) = c[t+1].as<double>();  // +1 skips 1st column (prod name)
      } else {
	S(i) = NAN;
      }
      mkt_id[i] = mkt_counter;
      ++i;
    }
    ++mkt_counter;
  }

  /// fill delta (initial value log S_jt - log S_0t)
  delta.resize(num_prods * num_mkts);
  for (unsigned t = 0; t != num_mkts; ++t) {
    for (unsigned j = 0; j != num_prods - 1; ++j) {
      if (!std::isnan(S(t * num_prods + j))) {
	delta(t * num_prods + j) = std::log(S(t * num_prods + j)) -\
	  std::log(S((t+1) * num_prods - 1));
      } else {
	delta(t * num_prods + j) = NAN;
      }
    }
    delta((t+1) * num_prods - 1) = 0.;
  }

  /// fill X2 - precos
  X2.resize(num_prods * num_mkts, 1);
  std::string lanches_precos_query = "SELECT * FROM lanches_precos;";
  pqxx::result R_lp(N.exec(lanches_precos_query));
  i = 0;
  for (unsigned t = 0; t != num_mkts; ++t) {
    for (auto c = R_lp.begin(); c != R_lp.end(); ++c) {
      if (!c[t+1].is_null()) {
	X2(i, 0) = c[t+1].as<double>();  // +1 skips 1st column (prod name)
      } else {
	X2(i, 0) = NAN;
      }
      ++i;
    }
    X2(i, 0) = 0.;  // outside good
    ++i;
  }


  /// fill X1 - precos, area dummies e brand dummies
  X1.resize(num_prods * num_mkts, 1 + areas.size() + num_prods - 2);
  i = 0;
  unsigned area;
  for (unsigned t = 0; t != num_mkts; ++t) {
    area = t / num_periods;
    for (unsigned j = 0; j != num_prods; ++j) {
      // precos
      X1(i, 0) = X2(i, 0);
      // area dummies
      for (unsigned col = 1; col != 1 + areas.size(); ++col) {
	if (area == 0) {
	  X1(i, col) = 0;  // this make coefficients relative to Area I
	} else if (area != 0 && col == area) {
	  X1(i, col) = 1;
	} else {
	  X1(i, col) = 0;
	}
      }
      // brand dummies
      for (unsigned col = areas.size(); col != X1.size2(); ++col) {
	if ((j == col - areas.size()) && (j != num_prods - 1)) {
	  X1(i, col) = 1;  // coeffs are relative to outside good
	} else {
	  X1(i, col) = 0;
	}
      }
      ++i;
    }
  }

  /// fill Z - area dummies, brand dummies e instrumentos: precos de outras
  // areas no mesmo periodo
  Z.resize(num_prods * num_mkts, areas.size() + num_prods - 2 + areas.size()\
	   - 1);
  unsigned period;
  unsigned col;
  unsigned aux_row;
  i = 0;
  for (unsigned t = 0; t != num_mkts; ++t) {
    period = t % num_periods;
    for (unsigned j = 0; j != num_prods; ++j) {
      for (col = 1; col != X1.size2(); ++col) {
	Z(i, col-1) = X1(i, col);
      }
      col = X1.size2() - 1;
      for (unsigned k = 0; k != areas.size(); ++k) {
	aux_row = j + period * num_prods + k * (num_prods * num_periods);
	if (i != aux_row) {
	  if (!(std::isnan(X1(aux_row, 0)))) {
	    Z(i, col) = X1(aux_row, 0);
	  } else {
	    Z(i, col) = 0;
	  }
	++col;
	}
      }
      ++i;
    }
  }
}

void GenArrays::elim_nans()
{
  std::vector<unsigned> nan_rows;
  for (unsigned i = 0; i != S.size(); ++i) {
    if (std::isnan(S(i))) {
      nan_rows.push_back(i);
    }
  }
  ublas::vector<double> S_tmp;
  S_tmp.resize(S.size() - nan_rows.size());
  ublas::vector<double> delta_tmp;
  delta_tmp.resize(S.size() - nan_rows.size());
  ublas::matrix<double> X1_tmp;
  X1_tmp.resize(S.size() - nan_rows.size(), X1.size2());
  ublas::matrix<double> X2_tmp;
  X2_tmp.resize(S.size() - nan_rows.size(), X2.size2());
  ublas::matrix<double> Z_tmp;
  Z_tmp.resize(S.size() - nan_rows.size(), Z.size2());
  ublas::vector<unsigned> mkt_id_tmp;
  mkt_id_tmp.resize(S.size() - nan_rows.size());
  unsigned j = 0;
  unsigned k = 0;
  for (unsigned i = 0; i != S.size(); ++i) {
    if (i != nan_rows[k]) {
      S_tmp(j) = S(i);
      delta_tmp(j) = delta(i);
      for (unsigned col = 0; col != X1.size2(); ++col) {
	X1_tmp(j, col) = X1(i, col);
      }
      for (unsigned col = 0; col != X2.size2(); ++col) {
	X2_tmp(j, col) = X2(i, col);
      }
      for (unsigned col = 0; col != Z.size2(); ++col) {
	Z_tmp(j, col) = Z(i, col);
      }
      mkt_id_tmp(j) = mkt_id(i);
    ++j;
    } else {
      ++k;
    }
  }
  S.clear();
  delta.clear();
  X1.clear();
  X2.clear();
  Z.clear();
  mkt_id.clear();
  S = S_tmp;
  delta = delta_tmp;
  X1 = X1_tmp;
  X2 = X2_tmp;
  Z = Z_tmp;
  mkt_id = mkt_id_tmp;
}		     
