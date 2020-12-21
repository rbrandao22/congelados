import numpy as np
from numpy.linalg import inv

from BLP import BLP
from BLP import obj_load
from BLP import persist

class Diagnostics(BLP):
    def __init__(self, ns, num_periods, num_prods, areas_size, arrays_dir,\
                 params_dir, theta2, contract_tol):
        super().__init__(ns, num_periods, num_prods, areas_size,\
                         arrays_dir, params_dir, theta2, contract_tol)

    def covariance(self, params_file):
        self.theta2 = obj_load(self.params_dir, params_file)
        # recalc other vars
        self.calc_phi_inv()
        self.contraction()
        self.calc_theta1()
        self.omega = self.delta - self.X1 @ self.theta1
        # variance calc
        n = int(self.delta.shape[0])
        omega_matrix = np.zeros([self.Z.shape[1], self.Z.shape[1]])
        Q_xz = 1./n * self.X1.T @ self.Z
        Q_zz = 1./n * self.Z.T @ self.Z
        Q_zx = 1./n * self.Z.T @ self.X1
        for jt in range(n):
            Z_jt = self.Z[jt].T[:, np.newaxis]
            omega_matrix += Z_jt @ Z_jt.T * self.omega[jt] ** 2
            
        omega_matrix /= n
        var_theta1 = inv(Q_xz @ inv(Q_zz) @ Q_zx) @\
            (Q_xz @ inv(Q_zz) @ omega_matrix @ inv(Q_zz) @ Q_zx) @\
            inv(Q_xz @ inv(Q_zz) @ Q_zx)
        std_dev_theta1 = np.diagonal(var_theta1) ** .5
        persist(self.params_dir + "std_dev_theta1", std_dev_theta1)

    def desc_stats(self):
        # init averages and std dev from data
        avg_prices = {} # avgs for each area
        avg_prod_prices = np.zeros(self.num_prods) # overall prod avgs
        avg_shares = {}
        std_prices = {}
        avg_std_prices = np.zeros(self.num_prods)
        
        # averages from counterfactual new eq prices and marginal costs
        X_f = obj_load(self.params_dir, "X_f_hybr")
        cfeq_idxs = obj_load(self.params_dir, "cfeq_idxs_hybr")
        cfeq_idxs = np.concatenate(cfeq_idxs)
        aux_diff = X_f - self.X2
        avg_cfeq_prices = {}
        diff_cfeq_prices = np.zeros(self.num_prods)
        mc = obj_load(self.params_dir, "mc")
        avg_mcs = np.zeros([self.num_prods, self.areas_size])
        avg_mcs_counter = np.zeros([self.num_prods, self.areas_size])
        
        for prod in np.sort(np.unique(self.prod_id))[:-1]:
            avg_prices[str(prod)] = {}
            avg_shares[str(prod)] = {}
            std_prices[str(prod)] = {}
            avg_cfeq_prices[str(prod)] = {}

            # average prices and shares
            for area in np.unique(self.area_id):
                avg_prices[str(prod)][str(area)] =\
                    np.mean(self.X2[np.intersect1d(np.where(self.prod_id ==\
                                                            prod),\
                                                   np.where(self.area_id ==\
                                                            area))])
                avg_shares[str(prod)][str(area)] =\
                    np.mean(self.S[np.intersect1d(np.where(self.prod_id ==\
                                                            prod),\
                                                   np.where(self.area_id ==\
                                                            area))])
                avg_cfeq_prices[str(prod)][str(area)] =\
                    np.mean(X_f[np.intersect1d(np.where(self.prod_id ==\
                                                        prod),\
                                               np.where(self.area_id ==\
                                                        area))])

            # std deviation in prices
            for period in np.sort(np.unique(self.period_id))[:-1]:
                std_prices[str(prod)][str(period)] =\
                    np.std(self.X2[np.intersect1d(np.where(self.prod_id ==\
                                                           prod),\
                                                  np.where(self.period_id ==\
                                                           period))])
            for period, std_price in std_prices[str(prod)].items():
                if not(np.isnan(std_price)):
                    avg_std_prices[prod] += std_price
                else:
                    avg_std_prices[prod] += 0.
            avg_std_prices[prod] /= len(std_prices[str(prod)])

            # diff between couterfactual eq prices and actual prices
            diff_cfeq_prices[prod] =\
                np.mean(aux_diff[np.intersect1d(np.where(self.prod_id==prod),\
                                                cfeq_idxs)])
            for area, price in avg_prices[str(prod)].items():
                avg_prod_prices[prod] += price
            avg_prod_prices[prod] /= len(avg_prices[str(prod)])

        # average marginal costs
        for mkt in np.unique(self.mkt_id):
            mkt_rows = np.where(self.mkt_id==mkt)[0][:-1]
            for i in range(mc[str(mkt)].shape[0]):
                prod = self.prod_id[mkt_rows[i]]
                area = self.area_id[mkt_rows[i]]
                avg_mcs[prod, area] += mc[str(mkt)][i]
                avg_mcs_counter[prod, area] += 1
        avg_mcs /= avg_mcs_counter

        persist(self.params_dir + "avg_prices", avg_prices, False)
        persist(self.params_dir + "avg_shares", avg_shares, False)
        persist(self.params_dir + "avg_cfeq_prices", avg_cfeq_prices,\
                False)
        persist(self.params_dir + "std_prices", std_prices, False)
        persist(self.params_dir + "avg_std_prices", avg_std_prices)
        persist(self.params_dir + "diff_cfeq_prices", diff_cfeq_prices)
        persist(self.params_dir + "avg_prod_prices", avg_prod_prices)
        persist(self.params_dir + "avg_mcs", avg_mcs)
