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
        # average prices
        avg_prices = {}
        for prod in np.sort(np.unique(self.prod_id))[:-1]:
            for area in np.unique(self.area_id):
                avg_prices[str(prod)+"_"+str(area)] = np.average(self.X2)
