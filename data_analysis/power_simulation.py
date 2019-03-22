'''
Created on Feb 22, 2017

@author: kfbz9246
'''
from math import log10, sqrt, log
from copy import deepcopy
from itertools import product, combinations
from pandas import DataFrame
from numpy import array, complex128
from numpy.polynomial.polynomial import polyfromroots, polyroots
from numpy.random import dirichlet, multinomial
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap
from matplotlib.gridspec import GridSpec

from statsmodels.api import GLM, families
'''
[.13,.1,.1,.19,.1,.15,.16,.07] : [0.099482, 0.130994, 0.130887, 0.158560, 0.132423, 0.117017, 0.127126, 0.103350]
[.42,.17,.01,.1,.15,.03,.11,.01] : [0.388828, 0.190676, 0.051094, 0.065458, 0.167363, 0.020109, 0.078639, 0.024684]
'''
cell_probabilities = [([.13,.1,.1,.19,.1,.15,.16,.07], [0.099482, 0.130994, 0.130887, 0.158560, 0.132423, 0.117017, 0.127126, 0.103350]),
([.42,.17,.01,.1,.15,.03,.11,.01],  [0.388828, 0.190676, 0.051094, 0.065458, 0.167363, 0.020109, 0.078639, 0.024684])]


ns = [100, 500, 1000, 2500]
over_dispersion_factors = [1.01, 2 , 3, 4, 5, 6, 7]
null_cell_probabilites = [0.388828, 0.190676, 0.051094, 0.065458, 0.167363, 0.020109, 0.078639, 0.024684]
alternative_cell_probabilites = [.42,.17,.01,.1,.15,.03,.11,.01]
number_simulations = 50000
dimention = 3


printable_interaction = 3.995498728
def get_cell_probabilites_dataframe(cell_probabilites):
    expected_cell_probabilites = []
    for a_product in product([0,1], repeat=dimention):
        next_expected_cell_probabilites_row_dict = {}
        for variable, value in zip(range(dimention), a_product):
            next_expected_cell_probabilites_row_dict[str(variable)] = value
        expected_cell_probabilites.append(next_expected_cell_probabilites_row_dict)
    expected_cell_probabilites = DataFrame(expected_cell_probabilites)
         
    
    
    
    for row_name, probability  in zip(list(expected_cell_probabilites.index), cell_probabilites):
        expected_cell_probabilites.loc[row_name, 'p'] = probability
    
    
    for row_name, row in expected_cell_probabilites.iterrows():
        if sum(row[range(dimention)])%2==0:
            expected_cell_probabilites.loc[row_name, 'even'] = True
        else:
            expected_cell_probabilites.loc[row_name, 'even'] = False
    return expected_cell_probabilites

temp_fig = figure(figsize=(9.25, 7), dpi=300)#, frameon=True, tight_layout=True)
grid_spec = GridSpec(2,1)
legend_lines = []
legend_labels = []
for alternative_index, (alternative_cell_probabilites, null_cell_probabilites) in enumerate(cell_probabilities):
    null_expected_cell_probabilites = get_cell_probabilites_dataframe(null_cell_probabilites)
    null_cell_index_is_even = null_expected_cell_probabilites['even']
    
    alternative_expected_cell_probabilites = get_cell_probabilites_dataframe(alternative_cell_probabilites)
    alternative_cell_index_is_even = alternative_expected_cell_probabilites['even']
    
    powers_lists = []
    for n in ns:
        powers = []
        for over_dispersion_factor in over_dispersion_factors:
            alpha = (n-over_dispersion_factor)/float(over_dispersion_factor-1)
            null_direchlet_parameters = array(null_expected_cell_probabilites['p'])*alpha
            alternative_direchlet_parameters = array(alternative_expected_cell_probabilites['p'])*alpha
        
        
            null_log_odds_ratios = []
            for i in range(number_simulations):
                p = dirichlet(null_direchlet_parameters)
                sample = multinomial(n, p)+.5
                odds_ratio = 1.0
                for is_even, sample_value in zip(null_cell_index_is_even, sample):
                    if is_even:
                        odds_ratio *= sample_value
                    else:
                        odds_ratio /= sample_value
            
                null_log_odds_ratios.append(log10(odds_ratio))
            null_log_odds_ratios.sort()
            lower_critical_value = null_log_odds_ratios[int(.025*number_simulations)]
            upper_critical_value = null_log_odds_ratios[int(.975*number_simulations)]
            
            
            number_rejected_null_hypothses = 0.0
            for i in range(number_simulations):
                p = dirichlet(alternative_direchlet_parameters)
                sample = multinomial(n, p)+.5
                odds_ratio = 1.0
                for is_even, sample_value in zip(null_cell_index_is_even, sample):
                    if is_even:
                        odds_ratio *= sample_value
                    else:
                        odds_ratio /= sample_value
                
                log_odds_ratio = log10(odds_ratio)
                if log_odds_ratio < lower_critical_value or log_odds_ratio > upper_critical_value:
                    number_rejected_null_hypothses += 1
            powers.append(number_rejected_null_hypothses/number_simulations)
        powers_lists.append(powers)
    
    
    temp_fig = figure(figsize=(9.25, 7), dpi=300, frameon=True)
    fig_ax = temp_fig.add_subplot(111)
    for power_list, sample_size in zip(powers_lists, ns):
        fig_ax.plot(over_dispersion_factors, power_list, label=str(sample_size))
    box = fig_ax.get_position()
    fig_ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    fig_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Sample Size')
    fig_ax.set_title('3-Way Interaction Parameter=%.2f' % printable_interaction)
    fig_ax.set_xlabel('Overdispersion Factor')
    fig_ax.set_ylabel('Power')

savefig('power.png', bbox_inches='tight')

