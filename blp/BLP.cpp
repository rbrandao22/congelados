#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <pqxx/pqxx>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "BLP.hpp"

namespace ublas = boost::numeric::ublas;


BLP::BLP(const unsigned num_periods, const unsigned num_bins_renda, const\
	 unsigned num_bins_idade, unsigned ns_, const\
	 std::vector<std::vector<unsigned>> areas, const unsigned max_threads)
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

  ns = ns_; // num of draws
  
  /// draw v & D
  
  // allocate and grab from db
  v.resize(boost::extents[ns][1][num_periods]);
  std::default_random_engine generator_n;
  std::normal_distribution<double> normal_dist(0., 1.);
  D.resize(boost::extents[ns][3][areas.size()]);
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
    for (unsigned draw_counter = 0; draw_counter < ns; ++draw_counter) {
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

      // fill v & D
      v[draw_counter][0][i] = normal_dist(generator_n);
      D[draw_counter][0][i] = std::log(renda);
      D[draw_counter][1][i] = std::pow(D[draw_counter][0][i], 2);
      D[draw_counter][2][i] = idade;
    }
  }
  N.commit();
  
  /// init parallel params
  unsigned hardware_threads = std::thread::hardware_concurrency();
  num_threads = std::min(hardware_threads != 0 ? hardware_threads : 1,\
			 max_threads);
}

void BLP::allocate()
{
  ublas::vector<double> auxV;
  auxV.resize(S.size());
  for (unsigned i = 0; i != num_threads; ++i) {
    s_calc.push_back(auxV);
  }
}

void BLP::calc_objective(std::vector<double> theta2_)
{
  theta2 = theta2_;
  std::vector<std::thread> threads;
  unsigned th, j, k, block_size;
  block_size = ns / num_threads;
  auto s_calc_L = [&] (unsigned th, unsigned begin, unsigned end) {
		     for (unsigned i =  begin; i < end; ++i) {
		       for (unsigned jt = 0; jt < S.size(); ++jt) {
			 s_calc[th][jt] = std::exp(delta(jt) + X2(jt) *\
						   (theta2[0] *\
						    v[i][0][area_id[jt]] +\
						    theta2[1] *\
						    D[i][0][area_id[jt]] +\
						    theta2[2] *\
						    D[i][1][area_id[jt]] +\
						    theta2[3] *\
						    D[i][2][area_id[jt]]));
		       }
		     }
		     // CONTINUE divide by mkt totals and take average from threads
		   };
  j = 0;
  for (th = 0; th < (num_threads - 1); ++th) {
    k = j + block_size;
    threads.push_back(std::thread(s_calc_L, th, j, k));
    j = k;
  }
  threads.push_back(std::thread(s_calc_L, th, j, ns));
  for (auto& thread : threads) {
    thread.join();
  }
}
