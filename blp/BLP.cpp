#include <cmath>
#include <fstream>
#include <iostream>
#include <pqxx/pqxx>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <Eigen/Dense>

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
  v.resize(boost::extents[ns][1][areas.size()]);
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
  s_ijt.resize(ns, S.size());
  ublas::matrix<double> auxM1;
  auxM1.resize(S.size(), S.size());
  for (unsigned i = 0; i != num_threads; ++i) {
    Ddelta1.push_back(auxM1);
    Ddelta1_aux.push_back(auxM1);
  }
  ublas::matrix<double> auxM2;
  auxM2.resize(S.size(), theta2.size());
  ublas::matrix<double> auxM3;
  auxM3.resize(S.size(), 1);
  for (unsigned i = 0; i != num_threads; ++i) {
    Ddelta2.push_back(auxM2);
    Ddelta2a_aux.push_back(auxM2);
    Ddelta2b_aux.push_back(auxM3);
  }
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
		       double mkt_sum = 0.;
		       unsigned aux_mkt_id = 0;
		       unsigned aux_init_jt = 0;
		       for (unsigned jt = 0; jt < S.size(); ++jt) {
		         if (aux_mkt_id == mkt_id[jt]) {
			   mkt_sum += s_aux[th][jt];
		         } else {
			   for (unsigned jt2 = aux_init_jt; jt2 < jt; ++jt2) {
			     s_aux[th][jt2] /= mkt_sum;
			   }
			   mkt_sum = 0.;
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
			 s_ijt(i, jt) = s_aux[th][jt];
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

void BLP::contraction()
{
  // increase tol params
  unsigned iters_nbr1 = 100;  // iterations before first increase
  unsigned iters_nbr2 = 50;  // number of iters for further increases
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
    if (ctol_inc) {
      if (iter_count == iters_nbr1 || (iter_count > iters_nbr1 &&\
				       (iter_count - iters_nbr1) %\
				       iters_nbr2 == 0)) {
	contract_tol *= tol_factor;
      }
    }
  }
}

void BLP::calc_phi_inv()
{
  phi_inv = ublas::prod(ublas::trans(Z), Z);
  Eigen::MatrixXd phi_inv_(phi_inv.size1(), phi_inv.size2());
  for (unsigned i = 0; i < phi_inv.size1(); ++i) {
    for (unsigned j = 0; j < phi_inv.size2(); ++j) {
      phi_inv_(i, j) = phi_inv(i, j);
    }
  }
  phi_inv_ = phi_inv_.inverse();
  for (unsigned i = 0; i < phi_inv.size1(); ++i) {
    for (unsigned j = 0; j < phi_inv.size2(); ++j) {
      phi_inv(i, j) = phi_inv_(i, j);
    }
  }
}

void BLP::calc_theta1()
{
  this->calc_phi_inv();
  ublas::matrix<double> aux_mat1;
  ublas::matrix<double> aux_mat1_inv;
  ublas::matrix<double> aux_mat2;
  // aux_mat1 = (X1'Z*phi_inv*Z'X1)^(-1)
  aux_mat1 = ublas::prod(ublas::trans(X1), Z);
  aux_mat1 = ublas::prod(aux_mat1, phi_inv);
  aux_mat1 = ublas::prod(aux_mat1, ublas::trans(Z));
  aux_mat1 = ublas::prod(aux_mat1, X1);
  Eigen::MatrixXd aux_mat1_(aux_mat1.size1(), aux_mat1.size2());
  for (unsigned i = 0; i < aux_mat1.size1(); ++i) {
    for (unsigned j = 0; j < aux_mat1.size2(); ++j) {
      aux_mat1_(i, j) = aux_mat1(i, j);
    }
  }
  aux_mat1_ = aux_mat1_.inverse();
  for (unsigned i = 0; i < aux_mat1.size1(); ++i) {
    for (unsigned j = 0; j < aux_mat1.size2(); ++j) {
      aux_mat1(i, j) = aux_mat1_(i, j);
    }
  }
  // aux_mat2 = aux_mat1_inv*X1'Z*phi_inv*Z'
  aux_mat2 = ublas::prod(aux_mat1, ublas::trans(X1));
  aux_mat2 = ublas::prod(aux_mat2, Z);
  aux_mat2 = ublas::prod(aux_mat2, phi_inv);
  aux_mat2 = ublas::prod(aux_mat2, ublas::trans(Z));
  // theta1 = aux_mat2*delta
  theta1 = ublas::prod(aux_mat2, delta);
}

void BLP::calc_Ddelta()
{
  std::vector<std::thread> threads;
  unsigned th, j, k, block_size;
  block_size = ns / num_threads;
  auto Ddelta_L = [&] (unsigned th, unsigned begin, unsigned end) {
		    for (unsigned i =  begin; i < end; ++i) {
		      /* fill Ddelta1 (share derivatives w/ respect to delta) &
		         compute Ddelta2 (share derivatives w/ respect to\
                         theta2) summation terms */
		      for (unsigned jt = 0; jt < S.size(); ++jt) {
			for (unsigned jt2 = 0; jt2 < S.size(); ++jt2) {
			  if (jt == jt2) {
			    Ddelta1_aux[th](jt, jt2) = s_ijt(i, jt) * (1 -\
			      s_ijt(i, jt));
			  } else {
			    Ddelta1_aux[th](jt, jt2) = -1 * s_ijt(i, jt) *\
			      s_ijt(i, jt2);
			  }
			}
			// Ddelta2 first terms
			Ddelta2a_aux[th](jt, 0) = v[i][0][area_id[jt]] *\
			  s_ijt(i, jt) * X2(jt);
			Ddelta2a_aux[th](jt, 1) = D[i][0][area_id[jt]] *\
			  s_ijt(i, jt) * X2(jt);
			Ddelta2a_aux[th](jt, 2) = D[i][1][area_id[jt]] *\
			  s_ijt(i, jt) * X2(jt);
			Ddelta2a_aux[th](jt, 3) = D[i][2][area_id[jt]] *\
			  s_ijt(i, jt) * X2(jt);
			// second terms (under j summation, done below)
			Ddelta2b_aux[th](jt, 0) = s_ijt(i, jt) * X2(jt);
		      }
		      // fill Ddelta2 
		      unsigned aux_mkt_id = 0;
		      unsigned aux_init_jt = 0;
		      unsigned aux_end_jt;
		      bool mkt_end = false;
		      for (unsigned jt = 0; jt < S.size(); ++jt) {
			if (mkt_id[jt] != aux_mkt_id || jt == 0) {
			  aux_mkt_id = mkt_id[jt];
			  aux_init_jt = jt;
			  mkt_end = false;
			  while (!mkt_end) {
			    ++jt;
			    if (mkt_id[jt] != aux_mkt_id || jt == S.size() - 1) {
			      aux_end_jt = jt;
			      mkt_end = true;
			    }
			  }
			}
			for (unsigned jt2 = aux_init_jt; jt2 < aux_end_jt;\
			     ++jt2) {
			  if (jt2 != jt) {
			    Ddelta2a_aux[th](jt, 0) -= v[i][0][area_id[jt]] *\
			      s_ijt(i, jt) * Ddelta2b_aux[th](jt2, 0);
			    Ddelta2a_aux[th](jt, 1) -= D[i][0][area_id[jt]] *\
			      s_ijt(i, jt) * Ddelta2b_aux[th](jt2, 0);
			    Ddelta2a_aux[th](jt, 2) -= D[i][1][area_id[jt]] *\
			      s_ijt(i, jt) * Ddelta2b_aux[th](jt2, 0);
			    Ddelta2a_aux[th](jt, 3) -= D[i][2][area_id[jt]] *\
			      s_ijt(i, jt) * Ddelta2b_aux[th](jt2, 0);
			  }
			}
		      }
		      for (unsigned jt = 0; jt < S.size(); ++jt) {
			if (i == begin) {
			  for (unsigned jt2 = 0; jt2 < S.size(); ++jt2) {
			    Ddelta1[th](jt, jt2) = Ddelta1_aux[th](jt, jt2);
			  }
			  for (unsigned jt2 = 0; jt2 < theta2.size(); ++jt2) {
			    Ddelta2[th](jt, jt2) = Ddelta2a_aux[th](jt, jt2);
			  }
			} else {
			  for (unsigned jt2 = 0; jt2 < S.size(); ++jt2) {
			    Ddelta1[th](jt, jt2) += Ddelta1_aux[th](jt, jt2);
			  }
			  for (unsigned jt2 = 0; jt2 < theta2.size(); ++jt2) {
			    Ddelta2[th](jt, jt2) += Ddelta2a_aux[th](jt, jt2);
			  }
			}
		      }
		    }
		  };
  j = 0;
  for (th = 0; th < (num_threads - 1); ++th) {
    k = j + block_size;
    threads.push_back(std::thread(Ddelta_L, th, j, k));
    j = k;
  }
  threads.push_back(std::thread(Ddelta_L, th, j, ns));
  for (auto& thread : threads) {
    thread.join();
  }
  // take average from threads
  for (th = 1; th < num_threads; ++th) {
    Ddelta1[0] += Ddelta1[th];
    Ddelta2[0] += Ddelta2[th];
  }
  Ddelta1[0] /= ns;
  Ddelta2[0] /= ns;
  // invert first matrix using Eigen and multiply
  Ddelta = Ddelta1[0];
  Eigen::MatrixXd Ddelta_(Ddelta.size1(), Ddelta.size2());
  for (unsigned i = 0; i < Ddelta.size1(); ++i) {
    for (unsigned j = 0; j < Ddelta.size2(); ++j) {
      Ddelta_(i, j) = Ddelta(i, j);
    }
  }
  Ddelta_ = Ddelta_.inverse();
  for (unsigned i = 0; i < Ddelta.size1(); ++i) {
    for (unsigned j = 0; j < Ddelta.size2(); ++j) {
      Ddelta(i, j) = Ddelta_(i, j);
    }
  }
  Ddelta = ublas::prod(Ddelta, Ddelta2[0]);
}

void BLP::grad_calc()
{
  this->contraction();
  this->calc_theta1();
  // calc error term
  omega = delta - ublas::prod(X1, theta1);
  this->calc_Ddelta();
  grad_aux = 2 * ublas::prod(ublas::trans(Ddelta), Z);
  grad_aux = ublas::prod(grad_aux, phi_inv);
  grad_aux = ublas::prod(grad_aux, ublas::trans(Z));
  grad = ublas::prod(grad_aux, omega);
}

void BLP::objective_calc()
{
  obj_aux = ublas::prod(ublas::trans(omega), Z);
  obj_aux = ublas::prod(obj_aux, phi_inv);
  obj_aux = ublas::prod(obj_aux, ublas::trans(Z));
  obj_value = ublas::inner_prod(obj_aux, omega);
}

void BLP::gmm(double nr_tol, double step_size, const unsigned max_iter)
{
  // contraction tol increase params
  ctol_inc = true;
  double ctol_inc_size = .01;
  
  double grad_size;
  unsigned iter_nbr = 0;
  while (true) {
    this->grad_calc();
    bool halt_check = true;
    for (unsigned i = 0; i < grad.size(); ++i) {
      if (halt_check == true && std::abs(grad(i)) > nr_tol){
	halt_check = false;
      }
    }
    if (halt_check || ctol_inc == false) {
      break;
    }
    ctol_inc = false;
    for (unsigned i = 0; i < theta2.size(); ++i) {
      theta2(i) -= grad(i) * step_size;
      if (std::abs(grad(i)) > ctol_inc_size) {
	ctol_inc = true;
      }
    }
    grad_size = 0;
    for (unsigned i = 0; i < grad.size(); ++i) {
      grad_size += std::abs(grad(i));
    }
    grad_size /= grad.size();
    this->objective_calc();
    std::cout << "NR #iter: " << iter_nbr << "  Gradient size: " << grad_size\
	      << "  Objective value: " << obj_value << std::endl;
    if (iter_nbr == max_iter) {
      break;
    }
    ++iter_nbr;
  }
}

void BLP::persist(const std::string persist_file2)
{
  std::ofstream fdesc;
  fdesc.open(persist_file2);
  assert(fdesc.is_open());
  for (unsigned i = 0; i < theta1.size(); ++i) {
    fdesc << "theta1_" << i << ": " << theta1(i) << '\n';
  }
  for (unsigned i = 0; i < theta2.size(); ++i) {
    fdesc << "theta2_" << i << ": " << theta2(i) << '\n';
  }
  std::cout << "Finished params persistance in file " << persist_file2 <<\
    std::endl;
}
