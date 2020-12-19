import os
from contextlib import redirect_stdout

import numpy as np
from numpy.linalg import inv
from scipy.optimize import root

from BLP import BLP
from BLP import obj_load
from BLP import persist

class SupplyAnalysis(BLP):
    def __init__(self, ns, num_periods, num_prods, areas_size, arrays_dir,\
                 params_dir, theta2, contract_tol):
        super().__init__(ns, num_periods, num_prods, areas_size,\
                         arrays_dir, params_dir, theta2, contract_tol)

    # calc share derivatives w/ respect to prices (enter firms' FOCs) 
    def calc_share_derivs(self):
        # allocate share_derivs
        i = 0
        alpha_s_jti = self.s_jti[self.idxs, i]
        area = self.area_id[self.idxs[0]] # get mkt area
        alpha_s_jti *= self.theta1[0] + self.theta2[1] * self.D_0[i, area] +\
            self.theta2[2] * self.D_1[i, area] + self.theta2[3] *\
            self.D_2[i, area] + self.theta2[0] * self.v[i, area]
        self.share_derivs[str(self.mkt)] = np.outer(alpha_s_jti,\
                                               self.s_jti[self.idxs, i])
        np.fill_diagonal(self.share_derivs[str(self.mkt)],\
                         self.share_derivs[str(self.mkt)].diagonal() -\
                         alpha_s_jti)
        for i in range(1, self.ns):
            alpha_s_jti = self.s_jti[self.idxs, i]
            alpha_s_jti *= self.theta1[0] + self.theta2[1] *\
                self.D_0[i, area] + self.theta2[2] * self.D_1[i, area] +\
                self.theta2[3] * self.D_2[i, area] + self.theta2[0] *\
                self.v[i, area]
            self.share_derivs[str(self.mkt)] += np.outer(alpha_s_jti,\
                                                    self.s_jti[self.idxs, i])
            np.fill_diagonal(self.share_derivs[str(self.mkt)],\
                             self.share_derivs[str(self.mkt)].diagonal() -\
                             alpha_s_jti)
        self.share_derivs[str(self.mkt)] /= self.ns

    # brands ownership current structure
    def calc_Omega(self):
        self.Omega[str(self.mkt)] = np.outer(self.brf_id[self.idxs],\
                                             self.brf_id[self.idxs])
        np.fill_diagonal(self.Omega[str(self.mkt)],\
                         np.ones(self.Omega[str(self.mkt)].shape[0]))

    # marginal costs
    def calc_mc(self):
        self.calc_share_derivs()
        self.calc_Omega()
        self.mc[str(self.mkt)] = self.X2[self.idxs] - inv(self.Omega[str(self.mkt)] *\
                                                self.share_derivs[str(self.mkt)])\
                                                @ self.s_calc[self.idxs]

    # optimal pricing under new offer structure (one firm for each brand)
    def firms_FOCs(self):
        Omega = np.identity(self.Omega[str(self.mkt)].shape[0])
        self.X2[self.idxs] = self.mc[str(self.mkt)] +\
            inv(Omega * self.share_derivs[str(self.mkt)]) @\
            self.s_calc[self.idxs]
        self.X1[self.idxs, 0] = self.X2[self.idxs][:, 0]
        
    def supp_obj(self, X_p):
        self.X_p = X_p[:, np.newaxis] # (newaxis needed for hybr & lm methods)
        self.X2[self.idxs] = np.copy(self.X_p)
        self.X1[self.idxs, 0] = np.copy(self.X_p[:, 0])
        self.contraction()
        self.calc_share_derivs()
        self.firms_FOCs()
        diff = self.X2[self.idxs] - self.X_p
        max_diff = abs(max(diff.min(), diff.max(), key=abs))
        with open(self.params_dir+"hybr_diff", 'a') as diff_file:
            diff_file.write(str(self.mkt)+'\t'+str(max_diff)+"\n")
        print(max_diff)
        return (self.X2[self.idxs] - self.X_p)[:, 0]

    def supp_call(self, params_file, root_tol, root_max_iter):
        # initilize vars and clear log file
        self.theta2 = obj_load(self.params_dir, params_file)
        try:
            os.remove(self.params_dir+"hybr_diff")
        except OSError:
            pass
        X_0 = np.copy(self.X2) # initial prices
        X_f = np.copy(self.X2) # final prices
        self.contraction()
        self.calc_phi_inv()
        self.calc_theta1()
        self.Omega = {}
        self.mc = {}
        self.share_derivs = {}
        # calc marginal costs under prevailing offer structure
        for self.mkt in np.unique(self.mkt_id):
            self.idxs = np.sort(np.where(self.mkt_id == self.mkt)[0])[:-1]
            self.calc_mc()
        # calc new eq prices in each mkt
        for self.mkt in np.unique(self.mkt_id):
            self.idxs = np.sort(np.where(self.mkt_id == self.mkt)[0])[:-1]
            self.X2 = np.copy(X_0)
            self.X1[:, 0] = np.copy(X_0[:, 0])
            self.contraction()
            X_p = np.copy(self.X2[self.idxs])
            sol = root(self.supp_obj, X_p, method="hybr", tol = root_tol,\
                       options={"maxiter": root_max_iter})
            print("Results for mkt = " + str(self.mkt) + ":")
            print(sol.success)
            print(sol.x)
            print(sol.message)
            X_f[self.idxs] = np.copy(self.X2[self.idxs])
        persist(self.params_dir + "X_f_hybr", self.X_f)
