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
        # init averages and std from data
        avg_prices = {}
        avg_shares = {}
        std_prices = {}
        avg_std_prices = np.zeros(self.num_prods)
        # average from counterfactual new eq prices
        X_p = obj_load(self.params_dir, "X_p")
        avg_cfeq_prices = {}
        for prod in np.sort(np.unique(self.prod_id))[:-1]:
            avg_prices[str(prod)] = {}
            avg_shares[str(prod)] = {}
            std_prices[str(prod)] = {}
            avg_cfeq_prices[str(prod)] = {}
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
                    np.mean(X_p[np.intersect1d(np.where(self.prod_id ==\
                                                        prod),\
                                               np.where(self.area_id ==\
                                                        area))])
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
        persist(self.params_dir + "avg_prices", avg_prices, False)
        persist(self.params_dir + "avg_shares", avg_shares, False)
        persist(self.params_dir + "avg_cfeq_prices", avg_cfeq_prices,\
                False)
        persist(self.params_dir + "std_prices", std_prices, False)
        persist(self.params_dir + "avg_std_prices", avg_std_prices)
