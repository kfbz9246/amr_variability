'''
Created on Apr 18, 2017

@author: kfbz9246
'''
from math import log, sqrt
from copy import deepcopy
from itertools import product, combinations
from pandas import DataFrame
from numpy import array, complex128, mean
from numpy.polynomial.polynomial import polyfromroots, polyroots
from numpy.random import dirichlet, multinomial
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap

'''
[.13,.1,.1,.19,.1,.15,.16,.07] : [0.099482, 0.130994, 0.130887, 0.158560, 0.132423, 0.117017, 0.127126, 0.103350]
[.42,.17,.01,.1,.15,.03,.11,.01] : [0.388828, 0.190676, 0.051094, 0.065458, 0.167363, 0.020109, 0.078639, 0.024684]
'''
cell_probabilites = [0.388828, 0.190676, 0.051094, 0.065458, 0.167363, 0.020109, 0.078639, 0.024684]


dimention = 3

print 1-sum(cell_probabilites)
expected_cell_probabilites = []
for a_product in product([0,1], repeat=dimention):
    next_expected_cell_probabilites_row_dict = {}
    for variable, value in zip(range(dimention), a_product):
        next_expected_cell_probabilites_row_dict[str(variable)] = value
    expected_cell_probabilites.append(next_expected_cell_probabilites_row_dict)
expected_cell_probabilites = DataFrame(expected_cell_probabilites)
     
for row_name, probability  in zip(list(expected_cell_probabilites.index), cell_probabilites):
    expected_cell_probabilites.loc[row_name, 'p'] = probability

print expected_cell_probabilites

for variable in range(dimention):
    expected_cell_probabilites_copy = deepcopy(expected_cell_probabilites)
    print DataFrame(expected_cell_probabilites_copy.groupby(str(variable))['p'].sum())


for variable_set_size in range(2, dimention+1):
    for variable_set in combinations(range(dimention), variable_set_size):
        expected_cell_probabilites_copy = deepcopy(expected_cell_probabilites)
        variable_set = [str(variable) for variable in variable_set]
        expected_cell_probabilites_marginalized = DataFrame(expected_cell_probabilites_copy.groupby(variable_set)['p'].sum())
        odds_ratio = 1
        for row_name, row_value  in expected_cell_probabilites_marginalized.iterrows():
            if sum(row_name)%2==0:
                odds_ratio *= row_value['p']
            else:
                odds_ratio /= row_value['p']
        print str(variable_set) + ': ' + str(odds_ratio)


for row_name, row in expected_cell_probabilites.iterrows():
    if sum(row[range(dimention)])%2==0:
        expected_cell_probabilites.loc[row_name, 'even'] = True
    else:
        expected_cell_probabilites.loc[row_name, 'even'] = False
        
