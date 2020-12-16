import math
import os
import pdb
import random

import numpy as np
import pickle
import psycopg2


class GenArrays:
    '''
    Obtains data from db, transforms for econometric use as specified in model
    and saves for use in BLP class. Also generates populational variables
    related to random coefficients
    '''
    def __init__(self, num_periods, areas):
        try:
            conn = psycopg2.connect("dbname='congelados' user='postgres' " +\
                                    "host='localhost' password='passwd'")
        except:
            print("Unable to connect do database")
        cur = conn.cursor()
        
        # fill observed shares - S, and mkt_id
        shares_query = "SELECT * FROM lanches_shares;"
        cur.execute(shares_query)
        shares_data = cur.fetchall()
        num_mkts = num_periods * len(areas)
        num_prods = len(shares_data)
        self.S = np.empty([num_prods * num_mkts])
        self.S = self.S[:, np.newaxis]
        self.mkt_id = np.empty([num_prods * num_mkts], dtype=int)
        i = 0
        mkt_counter = 0
        for t in range(num_mkts):
            for row in shares_data:
                self.S[i] = row[t+1]
                self.mkt_id[i] = mkt_counter
                i += 1
            mkt_counter += 1

        # fill delta (initial value log S_jt = log S_0t)
        self.delta = np.empty([num_prods * num_mkts])
        self.delta = self.delta[:, np.newaxis]
        for t in range(num_mkts):
            for j in range(num_prods - 1):
                self.delta[t * num_prods + j] = np.log(self.S[t * num_prods +\
                                                              j]) -\
                    np.log(self.S[(t+1) * num_prods - 1])
            self.delta[(t+1) * num_prods - 1] = 0.

        # fill X2 - precos, and outgood_id
        self.X2 = np.empty([num_prods * num_mkts])
        self.X2 = self.X2[:, np.newaxis]
        self.outgood_id = np.zeros([num_prods * num_mkts], dtype=int)
        precos_query = "SELECT * FROM lanches_precos;"
        cur.execute(precos_query)
        precos_data = cur.fetchall()
        i = 0
        for t in range(num_mkts):
            for row in precos_data:
                self.X2[i] = row[t+1]
                i += 1
            self.X2[i] = 0.
            self.outgood_id[i] = 1
            i += 1

        # fill X1 - precos, area dummies e brand dummies & area_id
        self.X1 = np.empty([num_prods * num_mkts, 1 + len(areas) + num_prods -\
                            2])
        self.area_id = np.empty([num_prods * num_mkts], dtype=int)
        i = 0
        for t in range(num_mkts):
            area = t // num_periods
            for j in range(num_prods):
                self.area_id[i] = area
                # precos
                self.X1[i, 0] = self.X2[i]
                # area dummies
                for col in range(1, 1 + len(areas)):
                    if area == 0:
                        self.X1[i, col] = 0. # coeffs relative to Area I
                    elif area != 0 and col == area:
                        self.X1[i, col] = 1.
                    else:
                        self.X1[i, col] = 0.
                # brand dummies
                for col in range(len(areas), self.X1.shape[1]):
                    if j == col - len(areas) and j != num_prods -1:
                        self.X1[i, col] = 1.
                    else:
                        self.X1[i, col] = 0.
                i += 1

        # fill Z - area dummies, brand dummies e instrumentos: precos de outras\
        # areas no mesmo periodo
        self.Z = np.empty([num_prods * num_mkts, len(areas) + num_prods - 2 +\
                            len(areas) - 1])
        i = 0
        for t in range(num_mkts):
            period = t % num_periods
            for j in range(num_prods):
                for col in range(1, self.X1.shape[1]):
                    self.Z[i, col-1] = self.X1[i, col]
                col = self.X1.shape[1] - 1
                for k in range(len(areas)):
                    aux_row = j + period * num_prods + k * (num_prods *\
                                                            num_periods)
                    if i != aux_row:
                        if math.isnan(self.X1[aux_row, 0]):
                            self.Z[i, col] = 0.
                        else:
                            self.Z[i, col] = self.X1[aux_row, 0]
                        col += 1
                i += 1

        # close database connection
        conn.commit()
        cur.close()
        conn.close()

    def elim_nans(self):
        nan_rows = np.argwhere(np.isnan(self.S))
        self.S = np.delete(self.S, nan_rows, axis=0)
        self.delta = np.delete(self.delta, nan_rows, axis=0)
        self.X1 = np.delete(self.X1, nan_rows, axis=0)
        self.X2 = np.delete(self.X2, nan_rows, axis=0)
        self.Z = np.delete(self.Z, nan_rows, axis=0)
        self.mkt_id = np.delete(self.mkt_id, nan_rows, axis=0)
        self.area_id = np.delete(self.area_id, nan_rows, axis=0)
        self.outgood_id = np.delete(self.outgood_id, nan_rows, axis=0)

        # rescale prices
        self.X1[:, 0] /= 100.
        self.X2[:, 0] /= 100.
    
    def gen_random_vars(self, areas, ns, num_bins_renda, num_bins_idade):
        try:
            conn = psycopg2.connect("dbname='congelados' user='postgres' " +\
                                    "host='localhost' password='passwd'")
        except:
            print("Unable to connect do database")
        cur = conn.cursor()
        
        # allocate v and D & grab data from db
        self.v = np.empty([ns, len(areas)])
        self.D_0 = np.empty([ns, len(areas)])
        self.D_1 = np.empty([ns, len(areas)])
        self.D_2 = np.empty([ns, len(areas)])
        renda_query = "SELECT * FROM renda;"
        cur.execute(renda_query)
        renda_data = cur.fetchall()
        idade_query = "SELECT * FROM idade;"
        cur.execute(idade_query)
        idade_data = cur.fetchall()
        dict_idades_query = "SELECT * FROM dict_idades;"
        cur.execute(dict_idades_query)
        dict_idades_data = cur.fetchall()

        #iterate areas
        for i in range(len(areas)):
            pop = np.empty(len(areas[i]))
            pop_share = np.empty(len(areas[i]))
            aux_pop_query = "SELECT populacao FROM populacao WHERE id = "
            for j in range(len(areas[i])):
                pop_query = aux_pop_query + str(areas[i][j]) + ";"
                cur.execute(pop_query)
                pop[j] = float(cur.fetchall()[0][0])
            pop_total = np.sum(pop)
            pop_share = pop / pop_total
            pop_cum_share = np.concatenate([np.zeros(1), np.cumsum(pop_share)])

            # take draws
            for draw_counter in range(ns):
                # select estado from a given area or mkt
                draw_estado = random.uniform(0., 1.)
                for j in range(pop_cum_share.shape[0]-1):
                    if draw_estado >= pop_cum_share[j] and draw_estado <\
                       pop_cum_share[j+1]:
                        num_estado = areas[i][j]
                        
                # take draw for renda
                rendas_estado = renda_data[num_estado-1]
                draw_renda = random.uniform(0., 1.)
                for j in range(num_bins_renda):
                    if j == 0 and draw_renda < rendas_estado[num_bins_renda+1]:
                        renda = rendas_estado[j+1]
                        j = num_bins_renda
                    elif j == num_bins_renda - 1 and draw_renda >=\
                         rendas_estado[num_bins_renda+j]:
                        renda = rendas_estado[j+1]
                    elif draw_renda >= rendas_estado[num_bins_renda+j] and\
                         draw_renda < rendas_estado[num_bins_renda+j+1]:
                        renda = rendas_estado[j+1]
                        j = num_bins_renda

                # take draw for idade
                idades_estado = idade_data[num_estado-1]
                draw_idade = random.uniform(0., 1.)
                for j in range(num_bins_idade):
                    if j == 0 and draw_idade < idades_estado[j+1]:
                        idade = dict_idades_data[j][1]
                        j = num_bins_idade
                    elif j == num_bins_idade - 1 and draw_idade >=\
                         idades_estado[j]:
                        idade = dict_idades_data[j][1]
                    elif draw_idade >= idades_estado[j] and draw_idade <\
                         idades_estado[j+1]:
                        idade = dict_idades_data[j][1]
                        j = num_bins_idade

                # rescale variables
                renda /= 1000.
                idade /= 100.

                # fill arrays
                self.v[draw_counter][i] = random.normalvariate(0., 1.)
                self.D_0[draw_counter][i] = np.log(renda)
                self.D_1[draw_counter][i] = np.log(renda)**2
                self.D_2[draw_counter][i] = idade

        # close database connection
        conn.commit()
        cur.close()
        conn.close()

    def save_arrays(self,save_dir):
        save_data = {"S": self.S, "delta": self.delta, "X1": self.X1, "X2":\
                     self.X2, "Z": self.Z, "mkt_id": self.mkt_id, "area_id":\
                     self.area_id, "outgood_id": self.outgood_id, "v": self.v,\
                     "D_0": self.D_0, "D_1": self.D_1, "D_2": self.D_2}
        for name, array in save_data.items():
            filename = save_dir + name + ".pkl"
            try:
                os.remove(filename)
            except OSError:
                pass
            with open(filename, 'wb') as array_file:
                pickle.dump(array, array_file)
            print(name + " array made persistent in file")

        
if __name__ == "__main__":
    num_periods = 18;
    areas = [[10, 11, 12, 13, 14, 15, 16], [17, 18, 19], [19], [20], [20],\
             [21, 22, 23], [24, 26, 27]]
    ns = 100  # number of draws from population
    num_bins_renda = 7
    num_bins_idade = 12
    save_dir = "results/arrays/"
    inst = GenArrays(num_periods, areas)
    inst.elim_nans()
    inst.gen_random_vars(areas, ns, num_bins_renda, num_bins_idade)
    inst.save_arrays(save_dir)
