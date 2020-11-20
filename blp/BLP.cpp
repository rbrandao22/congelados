#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "BLP.hpp"

namespace ublas = boost::numeric::ublas;


BLP::BLP(const std::string initguess_f, const double min_share_, const\
	 double contract_tol_, const unsigned max_iter_contract_, const double\
	 penalty_param1_, const unsigned penalty_param2_, const double\
	 init_tetra_size1_, const double init_tetra_size2_, const unsigned\
	 max_threads)
{
  min_share = min_share_;
  contract_tol = contract_tol_;
  penalty_param1 = penalty_param1_;
  penalty_param2 = penalty_param2_;
  init_tetra_size1 = init_tetra_size1_;
  init_tetra_size2 = init_tetra_size2_;
  max_iter_contract = max_iter_contract_;
  // count params
  params_nbr = 0;
  std::ifstream ifs_ig(initguess_f);
  std::string aux_line;
  while (std::getline(ifs_ig, aux_line)) {
    ++params_nbr;
  }
  // allocate P's and y
  ublas::vector<double> auxP;
  auxP.resize(params_nbr);
  for (unsigned i = 0; i <= params_nbr+3; ++i) {
    P.push_back(auxP); /* init P0, ..., Pn+1 + (NelderMead only)
                          P^bar (P[params_nbr+1]), P* (...+2), P** (...+3) */
    y.push_back(0.); // init y
  }
  // initialize P[0]
  unsigned i = 0;
  std::ifstream ifs_ig2(initguess_f);
  while (std::getline(ifs_ig2, aux_line)) {
    if (!aux_line.empty()) {
	P[0][i] = std::stod(aux_line);
    } else {
	throw std::runtime_error("check init guess file, aborting...");
    }
    ++i;
  }
  // init parallel params
  unsigned hardware_threads = std::thread::hardware_concurrency();
  num_threads = std::min(hardware_threads != 0 ? hardware_threads : 1,\
			 max_threads);
}

void BLP::allocate()
{
  grad.resize(params_nbr);
  // Initialize N, unobs util, and allocate other ublas objs
  N = s_obs_wg.size();
  ublas::vector<double> auxV;
  for (unsigned i = 0; i <= params_nbr+3; ++i) {
    xi0.push_back(auxV);
    xi0[i].resize(N);
    std::fill(xi0[i].data().begin(), xi0[i].data().end(), 0.);
    xi1.push_back(auxV);
    xi1[i].resize(N);
    std::fill(xi1[i].data().begin(), xi1[i].data().end(), 0.);
    s_aux1.push_back(auxV);
    s_aux1[i].resize(N);
    s_aux2.push_back(auxV);
    s_aux2[i].resize(N);
    D1.push_back(auxV);
    D1[i].resize(N);
    D2.push_back(auxV);
    D2[i].resize(N);
    s_calc.push_back(auxV);
    s_calc[i].resize(N);
    std::fill(s_calc[i].data().begin(), s_calc[i].data().end(), .1);
    ln_s_obs.push_back(auxV);
    ln_s_obs[i].resize(N);
  }
  // calc initial y0
  this->calc_objective(0);
}

void BLP::calc_objective(unsigned pt, bool add_penalty)
{
  // initialization
  std::fill(xi0[pt].data().begin(), xi0[pt].data().end(), 0.);
  std::vector<std::thread> threads;
  unsigned j, k, block_size;
  block_size = N / num_threads;
  // observe s_obs = s_obs_wg * pop_ave * mu (s_obs_wg is within group share)
  auto ln_s_obs_L = [&] (unsigned begin, unsigned end) {
		      for (unsigned i = begin; i < end; ++i) {
			ln_s_obs[pt][i] = std::max(std::log(s_obs_wg[i] *\
							    (pop_ave[i] / 1e6)\
							    * P[pt][22]),\
						   std::log(min_share *\
							    (pop_ave[i] / 1e6)\
							    * P[pt][22]));
		    }
		  };
  j = 0;
  for (unsigned i = 0; i < (num_threads - 1); ++i) {
    k = j + block_size;
    threads.push_back(std::thread(ln_s_obs_L, j, k));
    j = k;
  }
  threads.push_back(std::thread(ln_s_obs_L, j, N));
  for (auto& thread : threads) {
    thread.join();
  }
  // other lambda (_L) functions for use in while loop below
  auto xi1_L = [&] (unsigned begin, unsigned end) {
		 for (unsigned i = begin; i < end; ++i) {
		   xi1[pt][i] = xi0[pt][i] + P[pt][21] *\
		     (ln_s_obs[pt][i] -	std::log(s_calc[pt][i]));
		   if (std::isnan(xi1[pt][i])) {
		     xi1[pt][i] = xi0[pt][i];
		     break;
		   }
		 }
	       };
  auto s_aux_L = [&] (unsigned begin, unsigned end) {
		   for (unsigned i = begin; i < end; ++i) {
		     s_aux1[pt][i] = std::exp((X(i, 0) * P[pt][0] + X(i, 1) *\
					       P[pt][1] + X(i, 2) * P[pt][2] +\
					       X(i, 3) * P[pt][3] + X(i, 4) *\
					       P[pt][4] + X(i, 5) * P[pt][5] +\
					       X(i, 6) * P[pt][6] + X(i, 7) *\
					       P[pt][7] + X(i, 8) * P[pt][8] +\
					       X(i, 9) * P[pt][9] +\
					       xi1[pt][i]) / P[pt][21]);
		     s_aux2[pt][i] = std::exp((X(i, 0) * P[pt][10] + X(i, 1) *\
					       P[pt][11] + X(i, 2) * P[pt][12]\
					       + X(i, 3) * P[pt][13] + X(i, 4)\
					       * P[pt][14] + X(i, 5) *\
					       P[pt][15] + X(i, 6) *\
					       P[pt][16] + X(i, 7) * P[pt][17]\
					       + X(i, 8) * P[pt][18] + X(i, 9)\
					       * P[pt][19] + xi1[pt][i]) /\
					      P[pt][21]);
		   }
		 };
  auto s_calc_L = [&] (unsigned begin, unsigned end) {
		    for (unsigned i = begin; i < end; ++i) {
		      s_calc[pt][i] = P[pt][20] * ((s_aux1[pt][i] / D1[pt][i]) *\
						   (std::pow(D1[pt][i],\
							     P[pt][21]) / \
						    (1 + std::pow(D1[pt][i],\
								  P[pt][21]))))\
			+ (1 - P[pt][20]) * ((s_aux2[pt][i] / D2[pt][i]) *\
					     (std::pow(D2[pt][i], P[pt][21]) /\
					      (1 + std::pow(D2[pt][i],\
							    P[pt][21]))));
		    }
		  };
  /// Berry's Contraction 
  bool conv_check = 0;
  unsigned iter_contract = 0;
  while (!conv_check) {
    threads.clear();
    j = 0;
    for (unsigned i = 0; i < (num_threads - 1); ++i) {
      k = j + block_size;
      threads.push_back(std::thread(xi1_L, j, k));
      j = k;
    }
    threads.push_back(std::thread(xi1_L, j, N));
    for (auto& thread : threads) {
      thread.join();
    }
    // check for convergence
    for (unsigned i = 0; i <= N; ++i) {
      if (i == N || iter_contract == max_iter_contract) {
        conv_check = 1;
        break;
      }
      if (std::abs(xi1[pt][i] - xi0[pt][i]) < contract_tol) {
        continue;
      } else {
        break;
      }
    }
    /// calc shares
    // calc s_ind1 & s_ind2 (exp(Xbeta....))
    // check header file for params mapping to P
    threads.clear();
    j = 0;
    for (unsigned i = 0; i < (num_threads - 1); ++i) {
      k = j + block_size;
      threads.push_back(std::thread(s_aux_L, j, k));
      j = k;
    }
    threads.push_back(std::thread(s_aux_L, j, N));
    for (auto& thread : threads) {
      thread.join();
    }
    // calc D1 and D2
    double aux_D1 = 0;
    double aux_D2 = 0;
    unsigned initial_aux_i = 0;
    unsigned aux_mkt_id = 0;
    for (unsigned i = 0; i < N; ++i) {
      if (aux_mkt_id == mkt_id[i]) {
        aux_D1 += s_aux1[pt][i];
        aux_D2 += s_aux2[pt][i];
      } else {
        for (unsigned j = initial_aux_i; j < i; ++j) {
          D1[pt][j] = aux_D1;
          D2[pt][j] = aux_D2;
        }
        aux_D1 = s_aux1[pt][i];
        aux_D2 = s_aux2[pt][i];
        initial_aux_i = i;
        aux_mkt_id = mkt_id[i];
      }
      if (i == N - 1) {
        for (unsigned j = initial_aux_i; j <= i; ++j) {
          D1[pt][j] = aux_D1;
          D2[pt][j] = aux_D2;
        }
      }
    }
    // compute model shares
    threads.clear();
    j = 0;
    for (unsigned i = 0; i < (num_threads - 1); ++i) {
      k = j + block_size;
      threads.push_back(std::thread(s_calc_L, j, k));
      j = k;
    }
    threads.push_back(std::thread(s_calc_L, j, N));
    for (auto& thread : threads) {
      thread.join();
    }
    // Update unobs util
    xi0[pt] = xi1[pt];
    ++iter_contract;
  }
  /// Compute objective function: y = (1/N xi'Z)I(1/N Z'xi)
  ublas::identity_matrix<double> I (Z.size2());
  y[pt] = ublas::inner_prod(ublas::prod(1./N *\
					ublas::prod(ublas::trans(xi0[pt]), Z),\
					I), 1./N * ublas::prod(ublas::trans(Z),\
							       xi0[pt]));
  /// Add penalty (Nelder Mead only)
  if (add_penalty) {
    // (gamma or lambda not in [0,1], or mu < 0)
    double penalty = {0.};
    if (P[pt][20] < 0.) {
      penalty += std::pow(1 + P[pt][20], penalty_param2);
    } else if (P[pt][20] > 1.) {
      penalty += std::pow(P[pt][20], penalty_param2);
    } else if (P[pt][21] < 0.) {
      penalty += std::pow(1 + P[pt][21], penalty_param2);
    } else if (P[pt][21] > 1.) {
      penalty += std::pow(P[pt][21], penalty_param2);
    } else if (P[pt][22] < 0.) {
      penalty += std::pow(1 + P[pt][22], penalty_param2);
    }
    penalty *= penalty_param1;
    y[pt] += penalty;
  }
  
  // Ensure finiteness
  if (std::isnan(y[pt]))
    y[pt] = std::numeric_limits<double>::max();
}

void BLP::updatePs(const double inc)
{

  for (unsigned i = 1; i <= params_nbr; ++i) {
    P[i] = P[0];
    P[i][i-1] += inc;
  }
}

void BLP::grad_calc(const double inc, const double tol)
{
  halt_check = true;
  for (unsigned i = 0; i < params_nbr; ++i) {
    grad[i] = (y[i+1] - y[0]) / inc;
    if (halt_check == true && grad[i] > tol)
      halt_check = false;
  }
}

void BLP::step(const double step_size, const double max_step, const double\
	       step_factor, unsigned iter_nbr)
{
  ublas::vector<double> P_aux = P[0];
  double y_aux = y[0];
  std::vector<double> step;
  auto step_P0 = [&] (unsigned i) {
		   if (step[i] < -max_step) {
		     P[0][i] -= max_step;
		   } else if (step[i] > max_step) {
		     P[0][i] += max_step;
		   } else {
		     P[0][i] += step[i];
		   }
		   // constraints (gamma and lambda in [0,1], mu < 0)
		   if (i >= 20 && P[0][i] < 0.)
		     P[0][i] = .1;
		   if ((i == 20 || i == 21) && P[0][i] > 1.)
		     P[0][i] = .9;
		 };
		    
  for (unsigned i = 0; i < params_nbr; ++i) {
    step.push_back(- grad[i] * step_size);
    step_P0(i);
  }
  this->calc_objective(0);
  std::cout << "y value: " << y[0] << '\t' << "# of iterations: " << iter_nbr <<\
    '\r' << std::flush;
  /* step adjustment eliminated below
  if (y[0] < y_aux) {
    std::cout << "y value (increasing step size): " << y[0] << '\t' <<\
      "# of iterations: " << iter_nbr << '\r' << std::flush;
  } else {
    std::cout << "NR step failed (moving to NM), current y value: " << y_aux <<\
      '\t' << "# of iterations: " << iter_nbr << '\r' << std::flush;
  }

  if (y[0] < y_aux) {
    do_NelderMead = false;
    while (y[0] < y_aux) {
      y_aux = y[0];
      P[params_nbr] = P[0]; // save current best in aux var
      P[0] = P_aux; // revert P0 to increase step below
      P_aux = P[params_nbr];
      for (unsigned i = 0; i < params_nbr; ++i) {
	step[i] *= step_factor;
	step_P0(i);
      }
      this->calc_objective(0);
      if (y[0] < y_aux) {
	std::cout << "y value (increased step size): " << y[0] << '\t' <<\
	  "# of iterations: " << iter_nbr << '\r' << std::flush;
      } else {
	std::cout << "NR step increase failed, current y value: " << y_aux <<\
	  '\t' << "# of iterations: " << iter_nbr << '\r' << std::flush;
      }
    }
    y[0] = y_aux;
    P[0] = P_aux;
  } else if (y[0] >= y_aux) {
    y[0] = y_aux;
    P[0] = P_aux;
    do_NelderMead = true;
  } */
}

void BLP::updatePs_NM()
{
  // initialize P's, P[0] at center (tetrahedron expansion)
  for (unsigned i = 0; i < params_nbr+1; ++i) {
    for (unsigned j = 0; j < params_nbr; ++j) {
      if (i == 0 || j != i) {
	P[i][j] = P[0][j];
      } else {
	if (i < params_nbr-3) {
	  P[i][j] = P[0][j] + init_tetra_size1;
    	} else {
    	  P[i][j] = P[0][j] + init_tetra_size2;
    	}
      }
    }
  }
  P[params_nbr][0] += init_tetra_size1;
}

void BLP::nelder_mead(const double alpha, const double beta, const double gamma)
{
  // allocate & initialize
  unsigned h;
  unsigned l;
  bool add_penalty = true;

  // Basic functions
  auto get_h = [&] () {
		 std::vector<double>::iterator y_high;
		 y_high = std::max_element(y.begin(), y.begin()+params_nbr+1);
		 h = std::distance(y.begin(), y_high);
		 return h;
	       };
  auto get_l = [&] () {
		 std::vector<double>::iterator y_low;
		 y_low = std::min_element(y.begin(), y.begin()+params_nbr+1);
		 l = std::distance(y.begin(), y_low);
		 return l;
	       };
  auto P_bar = [&] () {
		 for (unsigned i =0; i < params_nbr; ++i) {
		   P[params_nbr+1][i] = 0.;
		 }
		 for (unsigned i = 0; i < params_nbr+1; ++i) {
		   if (i != h)
		     P[params_nbr+1] += P[i];
		 }
		 P[params_nbr+1] /= params_nbr;
		 return;
	       };
  auto reflection = [&] () {
		      P[params_nbr+2] = (1 + alpha) * P[params_nbr+1] - alpha *\
			P[h];
		      this->calc_objective(params_nbr+2, add_penalty);
		      return;
		    };
  auto expansion = [&] () {
		     P[params_nbr+3] = gamma * P[params_nbr+2] + (1 - gamma) *\
		       P[params_nbr+1];
		     this->calc_objective(params_nbr+3, add_penalty);
		     return;
		   };
  auto contraction = [&] () {
		       P[params_nbr+3] = beta * P[h] + (1 - beta) *\
			 P[params_nbr+1];
		       this->calc_objective(params_nbr+3, add_penalty);
		       return;
		     };
  /* where P^bar = P[params_nbr+1],
     P* = P[params_nbr+2],
     P** = P[params_nbr+3] */
  
  // procedural functions
  auto step_1 = [&] () {
		  get_h();
		  P_bar();
		  reflection();
		};

  // NM algorithm (see Computer Journal 1965 page 2)
  step_1();
  get_l();
  if (y[params_nbr+2] < y[l]) { // y* < yl ?
    expansion();
    if (y[params_nbr+3] < y[l]) { // y** < yl ?
      P[h] = P[params_nbr+3]; // Ph = P**
    } else {
      P[h] = P[params_nbr+2]; // Ph = P*
    }
  } else {
    bool test = true;
    for (unsigned i = 0; i < params_nbr+1; ++i) {
      if (y[params_nbr+2] < y[i])
	test = false;
    }
    if (test) { // y* > yi all i ?
      if (y[params_nbr+2] < y[h]) { // y* < yh ?
	P[h] = P[params_nbr+2]; // replace Ph by P* (intermediate)
	this->calc_objective(h, add_penalty);
      }
      contraction();
      if (y[params_nbr+3] > y[h]) { // y** > yh ?
	for (unsigned i = 0; i < params_nbr+1; ++i) {
	  P[i] = (P[i] + P[l]) / 2; // Pi = (Pi + Pl) /2 all i
	}
      } else {
	P[h] = P[params_nbr+3];  // Ph = P**
      }
    } else {
      P[h] = P[params_nbr+2]; // Ph = P*
    }
  }
}

void BLP::updateP0_NM(unsigned iter_nbr)
{
  // update P0 w/ tetrahedron centroid
  for (unsigned i = 0; i <= params_nbr; ++i) {
    P[params_nbr+1] += P[i];
  }
  P[params_nbr+1] /= (params_nbr+1);
  P[0] = P[params_nbr+1];
  std::cout << "Nelder Mead procedure finished at iter " << iter_nbr << '\r' <<\
    std::flush;
}

void BLP::variance()
{
  // TODO
}

void BLP::persist(const std::string persist_file2)
{
  std::ofstream fdesc;
  fdesc.open(persist_file2);
  assert(fdesc.is_open());
  for (unsigned i = 0; i < params_nbr; ++i) {
    fdesc << "P" << i << ": " << P[0][i] << '\n';
  }
  fdesc << "y: " << y[0] << '\n';
  fdesc.close();
  std::cout << "Finished params persistance in file " << persist_file2 <<\
    std::endl;
}

void BLP::persist_ig(const std::string initguess_f)
{
  std::ofstream fdesc;
  fdesc.open(initguess_f);
  assert(fdesc.is_open());
  for (unsigned i = 0; i < params_nbr; ++i) {
      fdesc << P[0][i] << '\n';
  }
  fdesc.close();
}
