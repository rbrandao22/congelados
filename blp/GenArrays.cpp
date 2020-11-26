#include <cmath>
#include <pqxx/pqxx>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "GenArrays.hpp"


GenArrays::GenArrays(const unsigned num_periods, const unsigned num_draws, const\
		     std::vector<std::vector<unsigned>> areas, const unsigned\
		     num_bins_renda, const unsigned num_bins_idade)
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
      for (unsigned col = 1 + areas.size(); col != X1.size2(); ++col) {
	if (j == col - 1 - areas.size() && j != num_prods - 1) {
	  X1(i, col) = 1;  // coeffs are relative to outside good
	} else {
	  X1(i, col) = 0;
	}
      }
      ++i;
    }
  }

  /// draw v
  v.resize(boost::extents[num_draws][1][num_periods]);
  std::default_random_engine generator_n;
  std::normal_distribution<double> normal_dist(0., 1.);
  for (unsigned i = 0; i != num_draws; ++i) {
    for (unsigned j = 0; j != num_periods; ++j) {
      v[i][0][j] = normal_dist(generator_n);
    }
  }

  /// draw D
  
  // allocate and grab from db
  D.resize(boost::extents[num_draws][3][areas.size()]);
  std::string renda_query = "SELECT * FROM renda;";
  pqxx::result R_renda(N.exec(renda_query));
  std::string idade_query = "SELECT * FROM idade;";
  pqxx::result R_idade(N.exec(idade_query));
  std::string dict_idades_query = "SELECT * FROM dict_idades;";
  pqxx::result R_dict_idades(N.exec(dict_idades_query));
  std::default_random_engine generator_u;
  std::uniform_real_distribution<double> uniform_dist(0., 1.);
  
  for (unsigned i = 0; i != areas.size(); ++i) {
    std::vector<double> share_pop;
    std::string pop_query = "SELECT populacao FROM populacao WHERE id = ";
    for (const auto& estado : areas[i]) {
      pqxx::result R_pop(N.exec(pop_query + std::to_string(estado) + ";"));
      auto c = R_pop.begin();
      share_pop.push_back(c[0].as<double>());
    }
    double share_pop_total = std::accumulate(share_pop.begin(), share_pop.end(), 0.);
    std::vector<double> share_pop_acum;
    share_pop_acum.push_back(0.);
    for (unsigned j = 1; j != share_pop.size(); ++j) {
      share_pop_acum.push_back(std::accumulate(share_pop.begin(),\
					       share_pop.begin() + j, 0.) /\
			       share_pop_total);
    }
    share_pop_acum.push_back(1.+1e-6);
    // select estado from a given area or mkt
    for (unsigned draw_counter = 0; draw_counter < num_draws; ++draw_counter) {
      double draw_estado = uniform_dist(generator_u);
      unsigned num_estado;
      for (unsigned j = 0; j != share_pop_acum.size()-1; ++j) {
	if (draw_estado >= share_pop_acum[j] && draw_estado <\
	    share_pop_acum[j+1])
	  num_estado = areas[i][j];
      }
      
      // take draw for renda
      auto c = R_renda.begin() + num_estado - 1;
      double draw_renda = uniform_dist(generator_u);
      double renda;
      for (unsigned j = 0; j < num_bins_renda; ++j) {
	if (j == 0 && draw_renda < c[num_bins_renda+1].as<double>()) {
	  renda = c[j+1].as<double>();
	  j = num_bins_renda;
	} else if (j == num_bins_renda-1 && draw_renda >=\
		   c[num_bins_renda+j].as<double>()) {
	  renda = c[j+1].as<double>();
	} else if (draw_renda >= c[num_bins_renda+j].as<double>() &&\
		   draw_renda < c[num_bins_renda+j+1].as<double>()) {
	  renda = c[j+1].as<double>();
	  j = num_bins_renda;
	}
      }
      
      // take draw for idade
      auto c2 = R_idade.begin() + num_estado - 1;
      double draw_idade = uniform_dist(generator_u);
      double idade;
      for (unsigned j = 0; j < num_bins_idade; ++j) {
	if (j == 0 && draw_idade < c2[j+1].as<double>()) {
	  auto c3 = R_dict_idades.begin() + j;
	  idade = c3[1].as<double>();
	  j = num_bins_idade;
	} else if (j == num_bins_idade-1 && draw_idade >=\
		   c2[j].as<double>()) {
	  auto c3 = R_dict_idades.begin() + j;
	  idade = c3[1].as<double>();
	} else if (draw_idade >= c2[j].as<double>() &&\
		   draw_idade < c2[j+1].as<double>()) {
	  auto c3 = R_dict_idades.begin() + j;
	  idade = c3[1].as<double>();
	  j = num_bins_idade;
	}
      }

      // fill D
      D[draw_counter][0][i] = std::log(renda);
      D[draw_counter][1][i] = std::pow(D[draw_counter][0][i], 2);
      D[draw_counter][2][i] = idade;
      auto debug = false;
      if (debug)
	std::cout << "debug" << std::endl; //DEBUG
    }
  }
  N.commit();
}
