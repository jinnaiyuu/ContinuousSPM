# Dimensionality reduction preprocessing for Continuous LAMP
# 
# 1. Convert to ranking
# 2. run some reduction algorithm
# 3. Output as a file

import sys
import numpy as np
# from sklearn import cluster
from pfa import PFA
# import csv


def load_data(in_file_name):
    data = np.loadtxt(in_file_name, delimiter=",")
    data = data.transpose()
    # data = np.around(data)
    return data


def rank(data):
    d = np.copy(data)
    ret = np.zeros_like(d)
    for i, array in enumerate(d):
        order = array.argsort()
        ranks = order.argsort()
        ranks_f = ranks.astype(float)
        ret[i] = ranks_f / len(array)
    return ret


def reduction(data, n_features):
    d = data.transpose()
    pfa = PFA(n_features=n_features)
    # red_model = cluster.FeatureAgglomeration(n_clusters=n_features)
    # print red_model
    pfa.fit(d)
    data_reduced = pfa.features_
    d_reduced = data_reduced.transpose()
    column_indices = pfa.indices_
    return (d_reduced, column_indices)


argv = sys.argv
argc = len(argv)

input_file = argv[1]
output_file = argv[2]

d = load_data(input_file)
# print d[1]
r = rank(d)
# print d[1]
print "RANK"
print r
red, column_indices = reduction(r, 2)
print "SHRINKED RANK"
print red

print "INDICES"
print column_indices
