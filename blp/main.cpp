#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <limits>
#include <iostream>
#include <string>
#include <thread>
#include <valarray>
#include <vector>
#include <boost/archive/text_oarchive.hpp>

#include "BLP.hpp"
#include "GenArrays.hpp"


/* current run options:
   1) argv[1] genarrays
   2) argv[1] estimation
   3) argv[1] genarrays, argv[2] estimation */


int main(int argc, char* argv[])
{
  auto chrono_start = std::chrono::steady_clock::now();

  /* PARAMETERS */

  // estimation periods - currently quarters only
  const std::vector<std::string> dates = {"201801", "201802", "201803"};
  // run identifier
  const std::string run_id = "2018Q1";

  /// Genarrays
  // price bins
  const std::valarray<double> bins = {0, 100, 200, 300, 400, 500, 600, 700, 800,\
				      900, 1e3, 1.25e3, 1.5e3, 1.75e3, 2e3, 3e3,\
				      5e5};
  // population threshold
  const unsigned pop_thres = {500000};

  // results directory
  const std::string results_dir = "results/";
  const std::string persist_file = results_dir + "arrays/" + run_id;
  const std::string persist_file2 = results_dir + "est_params/" + run_id;
  // initial guess file ((alpha, beta)_r, gamma, lambda, mu)
  const std::string initguess_f = results_dir + "init_guess";
  
  // maximum number of iterations
  const unsigned max_iter = {1000000};
  
  //// Estimation params
  // minimum 'observed shares' for numerical feasibility
  const double min_share = {1e-20};
  /// BLP contraction
  // tolerance (BJ10 suggests 1e-12)
  const double contract_tol = {1e-12};
  // maximum number of iterations
  const unsigned max_iter_contract = {1000};
  /// Newton Raphson params
  const double inc = {1e-4}; // dx for numerical gradient
  const double step_size = {1e-7};
  const double max_step = {1e-1};
  const double step_factor = {10};
  const double tol = {1e-8};
  /// Nelder Mead params
  // constrained optimization penalty
  const double penalty_param1 = {1e6};
  const unsigned penalty_param2 = {4}; // (must be even)
  // initial tetrahedron "size" for Nelder Mead procedure
  const double init_tetra_size1 = {.1};
  const double init_tetra_size2 = {.05}; // for constrained params (last 3)
  // NM coefficients
  const double alpha = {5}; // reflection, alpha > 0
  const double beta = {.5}; // contraction, beta in [0,1]
  const double gamma = {15}; // expansion, gamma > 1

  /* END OF PARAMETERS */

  if (argc > 1 && std::strcmp(argv[1], "genarrays") == 0) {
    GenArrays inst_GA(dates, bins, pop_thres);
    inst_GA.gen_instruments();
    inst_GA.gen_arrays();    
    // serialize
    {
        std::remove(persist_file.c_str());
        std::ofstream ofs(persist_file);
        assert(ofs.is_open());
        boost::archive::text_oarchive oa(ofs);
        oa << inst_GA;
    }

  } else if ((argc > 1 && std::strcmp(argv[1], "estimation") == 0) ||\
	     (argc > 2 && std::strcmp(argv[1], "genarrays") == 0 &&\
	      std::strcmp(argv[2], "estimation") == 0)) {
    // instantiate
    BLP inst_BLP(initguess_f, min_share, contract_tol, max_iter_contract,\
		 penalty_param1, penalty_param2, init_tetra_size1,\
		 init_tetra_size2);
    
    // deserialize
    {
        std::ifstream ifs(persist_file);
        assert(ifs.is_open());
        boost::archive::text_iarchive ia(ifs);
        ia >> inst_BLP;
	ifs.close();
    }
    inst_BLP.allocate();
    // GMM
    unsigned iter_nbr = 0;
    std::vector<std::thread> threads = {};
    while (true) {
      inst_BLP.updatePs(inc);
      threads.clear();
      for (unsigned pt = 1; pt <= inst_BLP.params_nbr; ++pt) {
        threads.push_back(std::thread(&BLP::calc_objective, std::ref(inst_BLP),\
				      pt, false));
      }
      for (auto& thread : threads) {
        thread.join();
      }
      inst_BLP.grad_calc(inc, tol);
      if (inst_BLP.halt_check) 
        break;
      inst_BLP.step(step_size, max_step, step_factor, iter_nbr);
      /*
      if (inst_BLP.do_NelderMead) {
	inst_BLP.updatePs_NM();
	threads.clear();
	for (unsigned pt = 1; pt <= inst_BLP.params_nbr; ++pt) {
	  threads.push_back(std::thread(&BLP::calc_objective,\
					std::ref(inst_BLP), pt, true));
	}
	for (auto& thread : threads) {
	  thread.join();
	}
	inst_BLP.nelder_mead(alpha, beta, gamma);
	inst_BLP.updateP0_NM(iter_nbr);
      }*/
      if (iter_nbr == max_iter)
	break;
      ++iter_nbr;
      std::remove(initguess_f.c_str());
      inst_BLP.persist_ig(initguess_f);
    }
    // compute variance
    inst_BLP.variance();
    // persist results
    inst_BLP.persist(persist_file2);
    std::cout << "# of iterations: " << iter_nbr << std::endl;

  } else {
    std::cout << "Invalid args!" << std::endl;
    throw std::runtime_error("aborting");
  }

  // finish chrono
  auto chrono_end = std::chrono::steady_clock::now();
  auto time_diff = chrono_end - chrono_start;
  std::string time_fpersist = results_dir + "ellapsed_time";
  std::remove(time_fpersist.c_str());
  std::ofstream fdesc_time;
  fdesc_time.open(time_fpersist);
  fdesc_time << "Last run duration: " << \
    std::chrono::duration_cast<std::chrono::minutes> (time_diff).count() << \
    " mins" << std::endl;
  fdesc_time.close();
  std::cout << "Elapsed time: " << \
    std::chrono::duration_cast<std::chrono::minutes> (time_diff).count() << \
    " mins" << std::endl;

  return 0;
}
