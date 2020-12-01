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

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "BLP.hpp"

namespace ublas = boost::numeric::ublas;


BLP::BLP(const unsigned num_periods, const unsigned num_bins_renda, const\
	 unsigned num_bins_idade, const std::vector<std::vector<unsigned>> areas,\
	 unsigned ns_, std::vector<double> theta2_, double contract_tol_,\
	 const unsigned max_threads)
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

  // initialization of vars
  ns = ns_;  // num of draws
  theta2.resize(theta2_.size());  // initial guess for theta2
  for (unsigned i = 0; i != theta2.size(); ++i) {
    theta2[i] = theta2_[i];
  }
  contract_tol = contract_tol_;
  
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
    s_aux.push_back(auxV);
    s_calc.push_back(auxV);
  }
  exp_delta1.resize(S.size());
  exp_delta2.resize(S.size());
}

void BLP::calc_shares()
{
  std::vector<std::thread> threads;
  unsigned th, j, k, block_size;
  block_size = ns / num_threads;
  auto s_calc_L = [&] (unsigned th, unsigned begin, unsigned end) {
		     for (unsigned i =  begin; i < end; ++i) {
		       for (unsigned jt = 0; jt < S.size(); ++jt) {
			 s_aux[th][jt] = std::exp(delta(jt) + X2(jt) *\
						   (theta2[0] *\
						    v[i][0][area_id[jt]] +\
						    theta2[1] *\
						    D[i][0][area_id[jt]] +\
						    theta2[2] *\
						    D[i][1][area_id[jt]] +\
						    theta2[3] *\
						    D[i][2][area_id[jt]]));
		       }
		       double mkt_sum = 0;
		       unsigned aux_mkt_id = 0;
		       unsigned aux_init_jt = 0;
		       for (unsigned jt = 0; jt < S.size(); ++jt) {
		         if (aux_mkt_id == mkt_id[jt]) {
			   mkt_sum += s_aux[th][jt];
		         } else {
			   for (unsigned jt2 = aux_init_jt; jt2 < jt; ++jt2) {
			     s_aux[th][jt2] /= mkt_sum;
			   }
			   mkt_sum = 0;
		       	   aux_init_jt = jt;
		       	   aux_mkt_id = mkt_id[jt];
		         }
		         if (jt == S.size() - 1) {
			   for (unsigned jt2 = aux_init_jt; jt2 <= jt; ++jt2) {
			     s_aux[th][jt2] /= mkt_sum;
			   }
		         }
		       }
		       for (unsigned jt = 0; jt < S.size(); ++jt) {
			 if (i == begin) {
			   s_calc[th][jt] = s_aux[th][jt];
			 } else {
			   s_calc[th][jt] += s_aux[th][jt];
			 }
		       }
		     }
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
  // take average from threads
  for (th = 1; th < num_threads; ++th) {
    s_calc[0] += s_calc[th];
  }
  s_calc[0] /= ns;
}

void BLP::contraction(bool increase_tol)
{
  // increase tol params
  unsigned iters_nbr1 = 50;  // iterations before first increase
  unsigned iters_nbr2 = 20;  // number of iters for further increases
  unsigned tol_factor = 1e1;  // increase factor
  // init threads
  std::vector<std::thread> threads;
  unsigned j, k, block_size;
  block_size = S.size() / num_threads;
  // contraction lambda function
  auto contract_L = [&] (unsigned begin, unsigned end) {
		      for (unsigned i = begin; i < end; ++i) {
			exp_delta2[i] = exp_delta1[i] * S[i] / s_calc[0][i];
			if (std::isnan(exp_delta2[i])) {
			  throw std::runtime_error("NAN in contract_L, check");
			}
			delta[i] = std::log(exp_delta2[i]);
		      }
		    };
  // initialization of exp_delta
  for (unsigned i = 0; i != delta.size(); ++i) {
    exp_delta1[i] = std::exp(delta[i]);
  }
  bool conv_check = false;
  unsigned iter_count = 0;
  while (!conv_check) {
    // calc shares
    this->calc_shares();
    // contraction step
    threads.clear();
    j = 0;
    for (unsigned i = 0; i < num_threads - 1; ++i) {
      k = j + block_size;
      threads.push_back(std::thread(contract_L, j , k));
      j = k;
    }
    threads.push_back(std::thread(contract_L, j, S.size()));
    for (auto& thread : threads) {
      thread.join();
    }
    // check for convergence
    for (unsigned i = 0; i <= exp_delta2.size(); ++i) {
      if (i == exp_delta2.size()) {
	conv_check = true;
	break;
      }
      if (std::abs(exp_delta2[i] - exp_delta1[i]) < contract_tol) {
	continue;
      } else {
	break;
      }
    }
    ++iter_count;
    if (iter_count == iters_nbr1 || (iter_count > iters_nbr1 && (iter_count -\
								 iters_nbr1) %\
				     iters_nbr2 == 0)) {
      contract_tol *= tol_factor;
    }
  }
}

void BLP::calc_phis()
{
  phi = ublas::prod(ublas::trans(Z), Z);
  phi_inv = ublas::identity_matrix<double> (phi.size1());
  ublas::permutation_matrix<size_t> pm(phi.size1());
  ublas::lu_factorize(phi, pm);
  ublas::lu_substitute(phi, pm, phi_inv);
  phi_inv *= -1;
}

void BLP::calc_theta1()
{
  this->calc_phis();
  ublas::matrix<double> aux_mat1;
  ublas::matrix<double> aux_mat1_inv;
  ublas::matrix<double> aux_mat2;
  // aux_mat1 = X1'Z*phi_inv*Z'X1
  aux_mat1 = ublas::prod(ublas::trans(X1), Z);
  aux_mat1 = ublas::prod(aux_mat1, phi_inv);
  aux_mat1 = ublas::prod(aux_mat1, ublas::trans(Z));
  aux_mat1 = ublas::prod(aux_mat1, X1);
  aux_mat1_inv = ublas::identity_matrix<double> (aux_mat1.size1());
  ublas::permutation_matrix<size_t> pm(aux_mat1.size1());
  ublas::lu_factorize(aux_mat1, pm);
  ublas::lu_substitute(aux_mat1, pm, aux_mat1_inv);
  aux_mat1_inv *= -1;
  // aux_mat2 = aux_mat1_inv*X1'Z*phi_inv*Z'
  aux_mat2 = ublas::prod(aux_mat1_inv, ublas::trans(X1));
  aux_mat2 = ublas::prod(aux_mat2, Z);
  aux_mat2 = ublas::prod(aux_mat2, phi_inv);
  aux_mat2 = ublas::prod(aux_mat2, ublas::trans(Z));
  // theta1 = aux_mat2*delta
  theta1 = ublas::prod(aux_mat2, delta);
}

void BLP::gmm()
{
  this->contraction();
  this->calc_theta1();
  // auto error_calc
  omega = delta - ublas::prod(X1, theta1);
  unsigned x = 0; //DEBUG
}
