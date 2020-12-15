import os

import pickle
import numpy as np
from scipy import optimize


def obj_load(save_dir, obj):
    obj_file = save_dir + obj + ".pkl"
    try:
        os.path.exists(obj_file)
        with open(obj_file, 'rb') as pkl_file:
            x = pickle.load(pkl_file)
    except:
        print("File " + obj_file + " not found or load failed")

    return x


class EmpModel:
    '''
    Estimation of demand params and supply side analysis using empirical model
    '''
    def __init__(self, ns, num_periods, areas_size, save_dir, theta2):
        # load objects from GenArrays
        self.S = obj_load(save_dir, "S")
        self.delta = obj_load(save_dir, "delta")
        self.X1 = obj_load(save_dir, "X1")
        self.X2 = obj_load(save_dir, "X2")
        self.Z = obj_load(save_dir, "Z")
        self.mkt_id = obj_load(save_dir, "mkt_id")
        self.area_id = obj_load(save_dir, "area_id")
        self.v = obj_load(save_dir, "v")
        self.D_0 = obj_load(save_dir, "D_0")
        self.D_1 = obj_load(save_dir, "D_1")
        self.D_2 = obj_load(save_dir, "D_2")

        # initialize other vars
        self.ns = ns
        self.areas_size = areas_size
        self.mkts_size = areas_size * num_periods
        self.theta2 = np.array(theta2)

    def calc_shares(self):
        self.s_ijt = np.empty([self.S.shape[0], self.v.shape[0]])
        for i in range(self.ns):
            for area in range(self.areas_size):
                idxs = np.where(self.area_id == area)[0]
                # calc shares numerator
                self.s_ijt[idxs, i] = np.exp(self.delta[idxs] + self.X2[idxs] *\
                                        (self.theta2[0] * self.v[i, area] +\
                                         self.theta2[1] * self.D_0[i, area] +\
                                         self.theta2[2] * self.D_1[i, area] +\
                                         self.theta2[3] * self.D_2[i, area])\
                                        )[:, 0]
            for mkt in range(self.mkts_size):
                idxs = np.where(self.mkt_id == mkt)[0]
                mkt_sum = np.sum(self.s_ijt[idxs, i])
                self.s_ijt[idxs, i] /= mkt_sum

        self.s_calc = self.s_ijt.sum(axis=1) / ns
        self.s_calc = self.s_calc[:, np.newaxis]
        print("debug")
                                                

    # test scipy inside class
    def fun(self, x):
        self.x = x
        return [self.x[0] + .5 * (self.x[0] - self.x[1])**3 - 1.,
               .5 * (self.x[1] - self.x[0])**3 + self.x[1]]

    def call(self):
        x = np.zeros(2)
        sol = optimize.root(self.fun, x)
        print(sol.x)
        

    
if __name__ == "__main__":
    # params
    ns = 100 # changing requires running GenArrays w/ same ns first
    num_periods = 18
    areas_size = 7
    save_dir = "results/arrays/"
    theta2 = [.1, .1, .1, .1]

    # instantiate and call
    inst = EmpModel(ns, num_periods, areas_size, save_dir, theta2)
    inst.calc_shares()
#    inst.call()
