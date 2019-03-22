'''
Created on Apr 15, 2017

@author: kfbz9246
'''
from math import log, sqrt
from copy import deepcopy
from itertools import product, combinations
from pandas import DataFrame
from numpy import array, complex128, mean, sum, inner
from numpy.polynomial.polynomial import polyfromroots, polyroots
from numpy.random import dirichlet, multinomial
from scipy.stats import norm, uniform
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap
from matplotlib.gridspec import GridSpec

from patsy import dmatrix
from numpy import matrix, array
from math import exp

null_coeffs = [-1.836, -.4055, .04886,-.8598,.4055,.4055,.4055, 0]
dimention = 3
n = 500
number_simulations = 500


variable_list = []
for variable_number in range(dimention):
    variable_name = 'X' + str(variable_number)
    variable_list.append(variable_name)

expected_cell_probabilites = []
for a_product in product([0,1], repeat=dimention):
    next_expected_cell_probabilites_row_dict = {}
    for variable_name, value in zip(variable_list, a_product):
        next_expected_cell_probabilites_row_dict[variable_name] = value
    expected_cell_probabilites.append(next_expected_cell_probabilites_row_dict)
expected_cell_probabilites = DataFrame(expected_cell_probabilites)


for row_name, row in expected_cell_probabilites.iterrows():
    if sum(row[range(dimention)])%2==0:
        expected_cell_probabilites.loc[row_name, 'even'] = True
    else:
        expected_cell_probabilites.loc[row_name, 'even'] = False


design_matrix = matrix(dmatrix('('+'+'.join(variable_list)+')**'+str(dimention), expected_cell_probabilites))
design_matrix_column_names = []
for design_matrix_column_id in range(design_matrix.shape[1]):
    column_name = 'D'+str(design_matrix_column_id)
    design_matrix_column_names.append(column_name)
    expected_cell_probabilites[column_name] = design_matrix[:, design_matrix_column_id]

gridspec = GridSpec(3,3)
temp_fig = figure(figsize=(9.25, 7), dpi=300, frameon=True)

for interaction_parameter_index in range(3):
    interaction_parameter = interaction_parameter_index *2   
    coeffs = deepcopy(null_coeffs)
    coeffs[-1] = interaction_parameter
    cell_data = expected_cell_probabilites['even'].to_frame()
    for row_id, row_data in expected_cell_probabilites.iterrows():
        cell_data.loc[row_id, 'p'] =  exp(inner(coeffs, row_data.loc[design_matrix_column_names])) 
    sum_p = sum(cell_data['p'])
    cell_data['p'] = cell_data['p']/sum_p
    print expected_cell_probabilites
    print interaction_parameter

    
    for over_dispersion_factor_index in range(3):
        over_dispersion_factor = 1 + 3*over_dispersion_factor_index
        
        cell_counts = (array(cell_data['p'])*n)+.5
        mean = 1.0
        for (cell_id, cell_is_even), cell_count in zip(cell_data['even'].iteritems(), cell_counts):
                if cell_is_even:
                    mean *= cell_count
                else:
                    mean /= cell_count
        mean = log(mean)
        
        standard_deviation = 0
        for cell_count in cell_counts:
            standard_deviation += 1.0/cell_count
        standard_deviation = sqrt(standard_deviation*over_dispersion_factor)
        
        
        theoretical_wald_parameters = []
        number_simulations_float = float(number_simulations)
        for index in range(1, number_simulations+1):
            theoretical_wald_parameters.append(norm.ppf(index/number_simulations_float, mean, standard_deviation))
        

        
        simulated_wald_parameters = []
        if over_dispersion_factor == 1:
            for i in range(number_simulations):
                sample = multinomial(n, cell_data['p'])+.5
                
               
                odds_ratio = 1.0
                for (cell_id, cell_is_even), cell_count in zip(cell_data['even'].iteritems(), sample):
                    if cell_is_even:
                        odds_ratio *= cell_count
                    else:
                        odds_ratio /= cell_count
                odds_ratio = log(odds_ratio)
       
                simulated_wald_parameters.append(odds_ratio)
        else:
            alpha = (n-over_dispersion_factor)/float(over_dispersion_factor-1)
            direchlet_parameters = array(cell_data['p'])*alpha
        
            for i in range(number_simulations):
                p = dirichlet(direchlet_parameters)
                sample = multinomial(n, p)+.5
                
                odds_ratio = 1.0
                for (cell_id, cell_is_even), cell_count in zip(cell_data['even'].iteritems(), sample):
                    if cell_is_even:
                        odds_ratio *= cell_count
                    else:
                        odds_ratio /= cell_count
                odds_ratio = log(odds_ratio)
            
                simulated_wald_parameters.append(odds_ratio)
        simulated_wald_parameters.sort()
        
        
        


        fig_ax = temp_fig.add_subplot(gridspec[over_dispersion_factor_index, interaction_parameter_index])
        
        
        fig_ax.scatter(theoretical_wald_parameters, simulated_wald_parameters, s=3, c='k')
        x_limits = fig_ax.get_xlim()
        fig_ax.plot(x_limits, x_limits, 'k', linewidth=1)
        fig_ax.set_xlim(x_limits)
        
        fig_ax.set_title('Model Parameter: ' + str(interaction_parameter) + '   Dispersion: '+ str(over_dispersion_factor), fontsize=8)
        if over_dispersion_factor_index == 2:
            fig_ax.set_xlabel('Theoretical', fontsize=8)
        else:
            fig_ax.set_xticklabels([])
            
        if interaction_parameter_index == 0:      
            fig_ax.set_ylabel('Empirical', fontsize=8)
        else:
            fig_ax.set_yticklabels([])

            
        
temp_fig.suptitle('Quantile-Quantile Plot')
temp_fig.savefig('qq.png', bbox_inches='tight')    