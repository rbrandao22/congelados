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
        self.contraction()
        self.calc_phi_inv()
        self.calc_theta1()
        # allocate share_derivs
        num_mkts = len(np.unique(self.mkt_id))
        self.share_derivs = {}
        for mkt in np.unique(self.mkt_id):
            idxs = np.sort(np.where(self.mkt_id == mkt)[0])[:-1] # outgood excluded
            i = 0
            alpha_s_jti = self.s_jti[idxs, i]
            area = self.area_id[idxs[0]] # get mkt area
            alpha_s_jti *= self.theta1[0] + self.theta2[1] * self.D_0[i, area] +\
                self.theta2[2] * self.D_1[i, area] + self.theta2[3] *\
                self.D_2[i, area] + self.theta2[0] * self.v[i, area]
            self.share_derivs[str(mkt)] = np.outer(alpha_s_jti,\
                                                   self.s_jti[idxs, i])
            np.fill_diagonal(self.share_derivs[str(mkt)],\
                             self.share_derivs[str(mkt)].diagonal() -\
                             alpha_s_jti)
            for i in range(1, self.ns):
                alpha_s_jti = self.s_jti[idxs, i]
                alpha_s_jti *= self.theta1[0] + self.theta2[1] *\
                    self.D_0[i, area] + self.theta2[2] * self.D_1[i, area] +\
                    self.theta2[3] * self.D_2[i, area] + self.theta2[0] *\
                    self.v[i, area]
                self.share_derivs[str(mkt)] += np.outer(alpha_s_jti,\
                                                        self.s_jti[idxs, i])
                np.fill_diagonal(self.share_derivs[str(mkt)],\
                                 self.share_derivs[str(mkt)].diagonal() -\
                                 alpha_s_jti)
            self.share_derivs[str(mkt)] /= self.ns

    # brands ownership current structure
    def calc_Omega(self):
        self.Omega = {}
        for mkt in np.unique(self.mkt_id):
            idxs = np.sort(np.where(self.mkt_id == mkt)[0])[:-1]
            self.Omega[str(mkt)] = np.outer(self.brf_id[idxs], self.brf_id[idxs])
            np.fill_diagonal(self.Omega[str(mkt)],\
                             np.ones(self.Omega[str(mkt)].shape[0]))

    # marginal costs
    def calc_mc(self):
        self.calc_share_derivs()
        self.calc_Omega()
        self.mc = {}
        for mkt in np.unique(self.mkt_id):
            idxs = np.sort(np.where(self.mkt_id == mkt)[0])[:-1]
            self.mc[str(mkt)] = self.X2[idxs] - inv(self.Omega[str(mkt)] *\
                                                    self.share_derivs[str(mkt)])\
                                                    @ self.s_calc[idxs]

    # optimal pricing under new offer structure (one firm for each brand)
    def firms_FOCs(self):
        for mkt in np.unique(self.mkt_id):
            idxs = np.sort(np.where(self.mkt_id == mkt)[0])[:-1]
            Omega = np.identity(self.Omega[str(mkt)].shape[0])
            self.X2[idxs] = self.mc[str(mkt)] + inv(Omega *\
                                                    self.share_derivs[str(mkt)])\
                                                    @ self.s_calc[idxs]
            self.X1[idxs, 0] = self.X2[idxs][:, 0]
        
    def supp_obj(self, X_p):
        self.X_p = X_p #[:, np.newaxis] (term needed for hybr & lm methods)
        self.X2 = np.copy(self.X_p)
        self.X1[:, 0] = np.copy(self.X_p[:, 0])
        self.calc_share_derivs()
        self.firms_FOCs()
        diff = self.X2 - self.X_p
        max_diff = abs(max(diff.min(), diff.max(), key=abs))
        with open(self.params_dir+"krylov_diff", 'a') as diff_file:
            diff_file.write(str(max_diff)+"\n")
        self.iter_counter += 1
        if self.iter_counter % 20 == 0:
            persist(self.params_dir + "X_p", self.X_p)
        return (self.X2 - self.X_p)[:, 0]

    def supp_call(self, params_file):
        self.theta2 = obj_load(self.params_dir, params_file)
        self.calc_mc()
        try:
            os.remove(self.params_dir+"krylov_diff")
        except OSError:
            pass
        self.iter_counter = 0
        X_p = obj_load(self.params_dir, "X_p")
        sol = root(self.supp_obj, X_p, method="krylov")
        with open(self.params_dir+"root_sol", 'w') as f:
            with redirect_stdout(f):
                print(sol.x)
        self.X_p = np.copy(self.X2)
        persist(self.params_dir + "X_p", self.X_p)
