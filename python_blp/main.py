import time

from GenArrays import GenArrays
from BLP import BLP
from Diagnostics import Diagnostics
from SupplyAnalysis import SupplyAnalysis

def main(run_option):
    ## common params
    ns = 200
    num_periods = 18

    ## GenArrays params
    areas = [[10, 11, 12, 13, 14, 15, 16], [17, 18, 19], [19], [20], [20],\
             [21, 22, 23], [24, 26, 27]]
    num_bins_renda = 7
    num_bins_idade = 12
    save_dir = "results_2e2/arrays/"

    ## other params
    num_prods = 7
    areas_size = 7
    arrays_dir = "results_2e2/arrays/"
    params_dir = "results_2e2/"
    theta2 = [-1.267939467881756599e+01, 9.280703169049563428e+01,\
              -2.276655479729512521e+01, -3.401315716142687506e+01]
    # Berry's contraction tolerance
    contract_tol = 1e-12
    # Direct Search params
    opt_tol = 1e-4
    step_size = 1e-2
    # Nelder-Mead params
    max_iter = 1000
    # Root finding params
    root_tol = 5e-3
    root_max_iter = 100
    
    if run_option == "genArrays":
        inst = GenArrays(num_periods, areas)
        inst.elim_nans()
        inst.gen_random_vars(areas, ns, num_bins_renda, num_bins_idade)
        inst.save_arrays(save_dir)

    # Nelder-Mead
    elif run_option == "estimationNM":
        inst = BLP(ns, num_periods, num_prods, areas_size, arrays_dir,\
                    params_dir, theta2, contract_tol)
        inst.gmm_NM(max_iter)

    # Direct Search
    elif run_option == "estimationDS":
        inst = BLP(ns, num_periods, num_prods, areas_size, arrays_dir,\
                    params_dir, theta2, contract_tol)
        inst.gmm_DS(opt_tol, step_size)

    elif run_option == "diagnostics":
        inst = Diagnostics(ns, num_periods, num_prods, areas_size, arrays_dir,\
                           params_dir, theta2, contract_tol)
        inst.covariance("theta2_NM")
        #inst.desc_stats()

    elif run_option == "supplyAnalysis":
        inst = SupplyAnalysis(ns, num_periods, num_prods, areas_size, arrays_dir,\
                           params_dir, theta2, contract_tol)
        inst.supp_call("theta2_NM", root_tol, root_max_iter)
    else:
        print("invalid run option")

        
if __name__ == "__main__":
    start_time = time.perf_counter()
    main("supplyAnalysis")
    print('Ellapsed time: ', (time.perf_counter() - start_time)/60)
