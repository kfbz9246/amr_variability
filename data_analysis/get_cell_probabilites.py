'''
Created on Apr 17, 2017

@author: kfbz9246
'''
from math import log10, sqrt
from copy import deepcopy
from itertools import product, combinations
from pandas import DataFrame
from numpy import array, complex128, zeros
from numpy.polynomial.polynomial import polyfromroots, polyroots
from numpy.random import dirichlet, multinomial
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap

from statsmodels.api import GLM, families



marginal_probabilites = [.5,.4,.6]
odds_ratio_lists = [[1.5, 1.5, 1.5], [1]]#[[1.2, .7,1.1], [2]]
dimention = len(marginal_probabilites)

def get_expected_cell_probabilites(marginal_probabilites, odds_ratio_lists):
    simulated_data = []
    for a_product in product([0,1], repeat=dimention):
        next_simulated_data_row_dict = {}
        for variable, value in zip(range(dimention), a_product):
            next_simulated_data_row_dict[str(variable)] = value
        simulated_data.append(next_simulated_data_row_dict)
    simulated_data = DataFrame(simulated_data)
    simulated_data.insert(0,'p', 1)
    
    for row_name, row in simulated_data.iterrows():
        for variable in range(dimention):
            variable_value = row[str(variable)]
            if variable_value == 0:
                simulated_data.loc[row_name, 'p'] *= marginal_probabilites[variable]
            else:
                simulated_data.loc[row_name, 'p'] *= 1-marginal_probabilites[variable]
    
    
    for variable_set_size, odds_ratio_list in zip(range(2, dimention), odds_ratio_lists):
        simulated_data_copy = deepcopy(simulated_data)
        for variable_set, odds_ratio in  zip(combinations(range(dimention), variable_set_size), odds_ratio_list):
            variable_set = [str(variable) for variable in variable_set]
            simulated_data_marginalized = DataFrame(simulated_data_copy.groupby(variable_set)['p'].sum())
            for simulated_data_margin in simulated_data_marginalized.index:
                if sum(simulated_data_margin)%2==0:
                    simulated_data_marginalized.loc[simulated_data_margin, 'even'] = True
                else:
                    simulated_data_marginalized.loc[simulated_data_margin, 'even'] = False
            
            odds_ratio_numerator_roots = []
            odds_ratio_denominator_roots = []
            for simulated_data_margin, simulated_data_margin_data in simulated_data_marginalized.iterrows():
                if simulated_data_margin_data['even']:
                    odds_ratio_numerator_roots.append(-1*simulated_data_margin_data['p'])
                else:
                    odds_ratio_denominator_roots.append(simulated_data_margin_data['p'])       
            odds_ratio_numerator_polynomial = polyfromroots(odds_ratio_numerator_roots)
            odds_ratio_denominator_polynomial = polyfromroots(odds_ratio_denominator_roots)*odds_ratio
            fudge_factors = list(polyroots(odds_ratio_numerator_polynomial- odds_ratio_denominator_polynomial))
            fudge_factor = None
            fudge_factor_tmp = 2
            for next_fudge_factor in fudge_factors:
                fudge_factor_real_component = None
                if isinstance(next_fudge_factor, complex128):
                    if next_fudge_factor.imag ==0:
                        fudge_factor_real_component = next_fudge_factor.real
                else:
                    fudge_factor_real_component = next_fudge_factor
                if fudge_factor_real_component != None:
                    if abs(fudge_factor_real_component) < abs(fudge_factor_tmp):
                        fudge_factor_tmp = fudge_factor_real_component
            
            if abs(fudge_factor_tmp) < 1:
                fudge_factor = fudge_factor_tmp
            
            fudge_factor_multipliers = {}
            for simulated_data_margin, simulated_data_margin_data in simulated_data_marginalized.iterrows():
                fudge_factor_proportion = fudge_factor/simulated_data_margin_data['p']
                if simulated_data_margin_data['even']:
                    fudge_factor_multipliers[simulated_data_margin] = 1 + fudge_factor_proportion
                else:
                    fudge_factor_multipliers[simulated_data_margin] = 1 - fudge_factor_proportion
                
            
            for row_name, row in simulated_data.iterrows():
                simulated_data.loc[row_name, 'p'] *= fudge_factor_multipliers[tuple(row.loc[variable_set])]
    
    variable_set = [str(variable) for variable in range(dimention)]
    odds_ratio = odds_ratio_lists[-1][0]
    for row_name, row in simulated_data.iterrows():
        if sum(row[variable_set])%2==0:
            simulated_data.loc[row_name, 'even'] = True
        else:
            simulated_data.loc[row_name, 'even'] = False
    
    odds_ratio_numerator_roots = []
    odds_ratio_denominator_roots = []
    for row_name, row in simulated_data.iterrows():
        if row['even']:
            odds_ratio_numerator_roots.append(-1*row['p'])
        else:
            odds_ratio_denominator_roots.append(row['p'])       
    
    odds_ratio_numerator_polynomial = polyfromroots(odds_ratio_numerator_roots)
    odds_ratio_denominator_polynomial = polyfromroots(odds_ratio_denominator_roots)*odds_ratio
    fudge_factors = list(polyroots(odds_ratio_numerator_polynomial- odds_ratio_denominator_polynomial))
    fudge_factor = None
    fudge_factor_tmp = 2
    for next_fudge_factor in fudge_factors:
        fudge_factor_real_component = None
        if isinstance(next_fudge_factor, complex128):
            if next_fudge_factor.imag ==0:
                fudge_factor_real_component = next_fudge_factor.real
        else:
            fudge_factor_real_component = next_fudge_factor
        if fudge_factor_real_component != None:
            if abs(fudge_factor_real_component) < abs(fudge_factor_tmp):
                fudge_factor_tmp = fudge_factor_real_component
    
    if abs(fudge_factor_tmp) < 1:
        fudge_factor = fudge_factor_tmp
    
    for row_name, row in simulated_data.iterrows():
        fudge_factor_proportion = fudge_factor/row['p']
        if row['even']:
            simulated_data.loc[row_name, 'p'] *= 1 + fudge_factor_proportion
        else:
            simulated_data.loc[row_name, 'p'] *= 1 - fudge_factor_proportion
        
    return simulated_data

expected_cell_probabilities = get_expected_cell_probabilites(marginal_probabilites, odds_ratio_lists)#.to_csv('/home/workspace/Test/data.csv')
get_expected_cell_probabilites(marginal_probabilites, odds_ratio_lists).to_csv('/home/workspace/Test/data.csv')
print expected_cell_probabilities
exit()

contingency_table = zeros([2]*dimention,dtype=float)

for index, row in expected_cell_probabilities.iterrows():
    contingency_table[tuple([i for i in row[range(1, dimention+1)]])] = row['p']
print contingency_table