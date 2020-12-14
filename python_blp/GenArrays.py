import numpy as np
import os
import psycopg2
import pdb


class GenArrays:
    '''
    Obtains data from db, transforms for econometric use as specified in model
    and saves for use in BLP class
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
        self.mkt_id = np.empty([num_prods * num_mkts])
        self.mkt_id = self.mkt_id[:, np.newaxis]
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

        # fill X2 - precos
        self.X2 = np.empty([num_prods * num_mkts])
        self.X2 = self.X2[:, np.newaxis]
        precos_query = "SELECT * FROM lanches_precos;"
        cur.execute(precos_query)
        precos_data = cur.fetchall()
        i = 0
        for t in range(num_mkts):
            for row in precos_data:
                self.X2[i] = row[t+1]
                i += 1
            self.X2[i] = 0.
            i += 1

        # fill X1 - precos, area dummies e brand dummies & area_id
        self.X1 = np.empty([num_prods * num_mkts, 1 + len(areas) + num_prods -\
                            2])
        self.area_id = np.empty([num_prods * num_mkts])
        self.area_id = self.area_id[:, np.newaxis]
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

        # fill Z
        conn.commit()
        cur.close()
        conn.close()

        
if __name__ == "__main__":
    num_periods = 18;
    areas = [[10, 11, 12, 13, 14, 15, 16], [17, 18, 19], [19], [20], [20],\
             [21, 22, 23], [24, 26, 27]]
    inst = GenArrays(num_periods, areas)
