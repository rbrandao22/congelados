import os
import pickle

import numpy as np
from numpy.linalg import inv
from scipy.optimize import minimize


def obj_load(directory, obj_name):
    obj_file = directory + obj_name + ".pkl"
    try:
        os.path.exists(obj_file)
        with open(obj_file, 'rb') as pkl_file:
            x = pickle.load(pkl_file)
    except:
        raise Exception("File " + obj_file + " not found or load failed")

    return x

def persist(filename, obj, numpy=True):
    with open(filename + ".pkl", 'wb') as pkl_file:
        pickle.dump(obj, pkl_file)
    if numpy:
        with open(filename, 'w') as np_file:
            np.savetxt(np_file, obj)
    print("File " + filename + " saved")
    

class BLP:
    '''
    Estimation of demand params and supply side analysis using empirical model
    '''
    def __init__(self, ns, num_periods, num_prods, areas_size, arrays_dir,\
                 params_dir, theta2, contract_tol):
        # load objects from GenArrays
        self.S = obj_load(arrays_dir, "S")
        self.delta = obj_load(arrays_dir, "delta")
        self.X1 = obj_load(arrays_dir, "X1")
        self.X2 = obj_load(arrays_dir, "X2")
        self.Z = obj_load(arrays_dir, "Z")
        self.mkt_id = obj_load(arrays_dir, "mkt_id")
        self.period_id = obj_load(arrays_dir, "period_id")
        self.prod_id = obj_load(arrays_dir, "prod_id")
        self.brf_id = obj_load(arrays_dir, "brf_id")
        self.area_id = obj_load(arrays_dir, "area_id")
        self.outgood_id = obj_load(arrays_dir, "outgood_id")
        self.v = obj_load(arrays_dir, "v")
        self.D_0 = obj_load(arrays_dir, "D_0")
        self.D_1 = obj_load(arrays_dir, "D_1")
        self.D_2 = obj_load(arrays_dir, "D_2")

        # initialize other vars
        self.ns = ns
        self.num_prods = num_prods
        self.areas_size = areas_size
        self.mkts_size = areas_size * num_periods
        self.theta2 = np.array(theta2)
        self.contract_tol = contract_tol
        self.params_dir = params_dir

    def calc_shares(self):
        self.s_jti = np.empty([self.S.shape[0], self.v.shape[0]])
        for i in range(self.ns):
            for area in range(self.areas_size):
                idxs = np.where(self.area_id == area)[0]
                # calc shares numerator
                self.s_jti[idxs, i] = np.exp(self.delta[idxs] + self.X2[idxs] *\
                                        (self.theta2[0] * self.v[i, area] +\
                                         self.theta2[1] * self.D_0[i, area] +\
                                         self.theta2[2] * self.D_1[i, area] +\
                                         self.theta2[3] * self.D_2[i, area])\
                                        )[:, 0]
            for mkt in range(self.mkts_size):
                idxs = np.where(self.mkt_id == mkt)[0]
                mkt_sum = np.sum(self.s_jti[idxs, i])
                self.s_jti[idxs, i] /= mkt_sum

        self.s_calc = self.s_jti.sum(axis=1) / self.ns
        self.s_calc = self.s_calc[:, np.newaxis]

    def contraction(self):
        conv_check = False
        idxs = np.where(self.outgood_id == 0)[0]
        exp_delta0 = np.exp(self.delta[idxs])
        while not conv_check:
            self.calc_shares()
            exp_delta1 = exp_delta0 * self.S[idxs] / self.s_calc[idxs]
            exp_delta0 = exp_delta1
            new_delta = np.log(exp_delta1)
            diff = new_delta - self.delta[idxs]
            self.delta[idxs] = new_delta
            if abs(max(diff.min(), diff.max(), key=abs)) < self.contract_tol:
                conv_check = True

    def calc_phi_inv(self):
        self.phi_inv = inv(self.Z.T @ self.Z)

    def calc_theta1(self):
        # (X1'Z*phi_inv*Z'X1)^(-1)*X1'Z*phi_inv*Z'*delta
        self.theta1 = inv(self.X1.T @ self.Z @ self.phi_inv @ self.Z.T @\
                          self.X1) @ self.X1.T @ self.Z @ self.phi_inv @\
                          self.Z.T @ self.delta

    # calc derivatives of mean values w/ respect to params (for grad eval)
    def calc_Ddelta(self):
        i = 0
        # initialize first matrix
        Ddelta1 = -1. * np.outer(self.s_jti[:, i], self.s_jti[:, i])
        np.fill_diagonal(Ddelta1, Ddelta1.diagonal() + self.s_jti[:, i])
        # init second
        Ddelta2 = np.empty([self.S.shape[0], 4])
        for area in range(self.areas_size):
            idxs = np.where(self.area_id == area)[0]
            Ddelta2[idxs, 0] = self.v[i, area] * self.s_jti[idxs, i] *\
                self.X2[idxs][:, 0]
            Ddelta2[idxs, 1] = self.D_0[i, area] * self.s_jti[idxs, i] *\
                self.X2[idxs][:, 0]
            Ddelta2[idxs, 2] = self.D_1[i, area] * self.s_jti[idxs, i] *\
                self.X2[idxs][:, 0]
            Ddelta2[idxs, 3] = self.D_2[i, area] * self.s_jti[idxs, i] *\
                self.X2[idxs][:, 0]
            for idx in idxs:
                idxs_mkt = np.where(self.mkt_id == self.mkt_id[idx])[0]
                idxs_mkt = idxs_mkt[:, np.newaxis]
                idxs_mkt = np.delete(idxs_mkt, np.where(idxs_mkt == idx))
                aux_sum = np.sum(self.X2[idxs_mkt][:, 0] *\
                                 self.s_jti[idxs_mkt, i])
                Ddelta2[idx, 0] -= self.v[i, area] * self.s_jti[idx, i] *\
                    aux_sum
                Ddelta2[idx, 1] -= self.D_0[i, area] * self.s_jti[idx, i] *\
                    aux_sum
                Ddelta2[idx, 2] -= self.D_1[i, area] * self.s_jti[idx, i] *\
                    aux_sum
                Ddelta2[idx, 3] -= self.D_2[i, area] * self.s_jti[idx, i] *\
                    aux_sum
                              
        for i in range(1, self.ns):
            # compute first
            Ddelta1 += -1 * np.outer(self.s_jti[:, i], self.s_jti[:, i])
            np.fill_diagonal(Ddelta1, Ddelta1.diagonal() + self.s_jti[:, i])
            # compute second
            for area in range(self.areas_size):
                idxs = np.where(self.area_id == area)[0]
                Ddelta2[idxs, 0] += self.v[i, area] * self.s_jti[idxs, i] *\
                    self.X2[idxs][:, 0]
                Ddelta2[idxs, 1] += self.D_0[i, area] * self.s_jti[idxs, i] *\
                    self.X2[idxs][:, 0]
                Ddelta2[idxs, 2] += self.D_1[i, area] * self.s_jti[idxs, i] *\
                    self.X2[idxs][:, 0]
                Ddelta2[idxs, 3] += self.D_2[i, area] * self.s_jti[idxs, i] *\
                    self.X2[idxs][:, 0]
                for idx in idxs:
                    idxs_mkt = np.where(self.mkt_id == self.mkt_id[idx])[0]
                    idxs_mkt = idxs_mkt[:, np.newaxis]
                    idxs_mkt = np.delete(idxs_mkt, np.where(idxs_mkt == idx))
                    aux_sum = np.sum(self.X2[idxs_mkt][:, 0] *\
                                     self.s_jti[idxs_mkt, i])
                    Ddelta2[idx, 0] -= self.v[i, area] * self.s_jti[idx, i] *\
                        aux_sum
                    Ddelta2[idx, 1] -= self.D_0[i, area] * self.s_jti[idx, i] *\
                        aux_sum
                    Ddelta2[idx, 2] -= self.D_1[i, area] * self.s_jti[idx, i] *\
                        aux_sum
                    Ddelta2[idx, 3] -= self.D_2[i, area] * self.s_jti[idx, i] *\
                        aux_sum

        Ddelta1 /= ns
        Ddelta2 /= ns
        self.Ddelta = -1 * inv(Ddelta1) @ Ddelta2

    def grad_calc(self):
        self.contraction()
        self.calc_theta1()
        # compute error term
        self.omega = self.delta - self.X1 @ self.theta1
        self.calc_Ddelta()
        self.grad = 2 * self.Ddelta.T @ self.Z @ self.phi_inv @ self.Z.T @\
            self.omega

    def objective_calc(self):
        self.obj_value = self.omega.T @ self.Z @ self.phi_inv @ self.Z.T @\
            self.omega

    # Direct Search GMM implementation
    def gmm_DS(self, opt_tol, step_size):
        self.calc_phi_inv() # homoscedastic assumption takes only one calc
        while True:
            self.grad_calc()
            print("Grad: ", self.grad[:, 0])
            self.objective_calc()
            print("Objective: ", self.obj_value[0][0])
            if np.all(np.absolute(self.grad) < np.ones(np.shape(self.grad)) *\
                      opt_tol):
                print("Finished optimization")
                break
            self.theta2 -= step_size * self.grad[:, 0]
            print("theta2: ", self.theta2)

        persist(self.params_dir + "theta1", self.theta1)
        persist(self.params_dir + "theta2", self.theta2)

    def objective_calc2(self, theta2):
        self.theta2 = theta2
        self.contraction()
        self.calc_theta1()
        # compute error term
        self.omega = self.delta - self.X1 @ self.theta1
        self.obj_value = self.omega.T @ self.Z @ self.phi_inv @ self.Z.T @\
            self.omega
        return self.obj_value

    def callbackF(self, theta2):
        print("Objective: ", self.obj_value[0][0])
        print("theta2: ", self.theta2)

    # Nelder Mead GMM implementation using SciPy
    def gmm_NM(self, max_iter):
        self.calc_phi_inv()
        sol = minimize(self.objective_calc2, self.theta2, method="Nelder-Mead",\
                       callback=self.callbackF, options={"maxiter": max_iter,\
                                                         "disp": True})
        sol.x
        persist(self.params_dir + "theta1_NM", self.theta1)
        persist(self.params_dir + "theta2_NM", self.theta2)
