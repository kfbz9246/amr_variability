'''
Created on Mar 15, 2017

@author: kfbz9246
'''
from math import log, sqrt
from copy import deepcopy
from itertools import product, combinations
from pandas import DataFrame
from numpy import array, complex128
from numpy.polynomial.polynomial import polyfromroots, polyroots
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, savefig, figlegend
from matplotlib.gridspec import GridSpec

cell_probabilities = [([.13,.1,.1,.19,.1,.15,.16,.07], [0.099482, 0.130994, 0.130887, 0.158560, 0.132423, 0.117017, 0.127126, 0.103350]),
([.42,.17,.01,.1,.15,.03,.11,.01],  [0.388828, 0.190676, 0.051094, 0.065458, 0.167363, 0.020109, 0.078639, 0.024684])]




ns = [100, 500, 1000, 2500]
over_dispersion_factors = [1.01, 2 , 3, 4, 5, 6, 7]
dimention = 3


def get_log_odds_ratio(cell_probabilites):
    expected_cell_probabilites = []
    for a_product in product([0,1], repeat=int(log(len(cell_probabilites), 2))):
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
    
    odds_ratio = 1.0
    for row_name, row_data in expected_cell_probabilites.iterrows():
        if row_data['even']:
            odds_ratio *= row_data['p']
        else:
            odds_ratio /= row_data['p']
    return log(odds_ratio)



temp_fig = figure(figsize=(9.25, 7), dpi=300)#, frameon=True, tight_layout=True)
grid_spec = GridSpec(2,1)
legend_lines = []
legend_labels = []
for alternative_index, (alternative_cell_probabilites, null_cell_probabilites) in enumerate(cell_probabilities):
    null_cell_probabilites = array(null_cell_probabilites)
    alternative_cell_probabilites = array(alternative_cell_probabilites)
    interaction = get_log_odds_ratio(alternative_cell_probabilites)


    
    powers_lists = []
    for n in ns:
        null_cell_counts = null_cell_probabilites*n
        alternative_cell_counts = alternative_cell_probabilites*n
        
        null_standard_deviation = 0
        for cell_count in null_cell_counts:
            null_standard_deviation += 1.0/cell_count
        null_standard_deviation = sqrt(null_standard_deviation)
        
        alternative_standard_deviation = 0
        for cell_count in alternative_cell_counts:
            alternative_standard_deviation += 1.0/cell_count
        alternative_standard_deviation = sqrt(alternative_standard_deviation)
        print alternative_standard_deviation
         
        powers = []
        for over_dispersion_factor in over_dispersion_factors:
            overdispersed_null_standard_deviation = sqrt(over_dispersion_factor) * null_standard_deviation
            overdispersed_alternative_standard_deviation = sqrt(over_dispersion_factor) * alternative_standard_deviation
            
            null_lower_critical_value = norm.ppf(.025, 0, overdispersed_null_standard_deviation)
            null_upper_critical_value = norm.ppf(.975, 0, overdispersed_null_standard_deviation)
    
            power = norm.cdf(null_lower_critical_value, interaction, overdispersed_alternative_standard_deviation)
            power += 1 - norm.cdf(null_upper_critical_value, interaction, overdispersed_alternative_standard_deviation)
            powers.append(power)
        powers_lists.append(powers)
    
    
    fig_ax = temp_fig.add_subplot(grid_spec[alternative_index,0])
    for power_list, sample_size in zip(powers_lists, ns):
        next_legend_line = fig_ax.plot(over_dispersion_factors, power_list)[0]
        if alternative_index == 0:
            legend_lines.append(next_legend_line)
            legend_labels.append(str(sample_size))

    
    #fig_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Sample_size')
    fig_ax.set_title('3-Way Interaction Parameter=%.2f' % interaction)
    fig_ax.set_xlabel('Overdispersion Factor')
    fig_ax.set_ylabel('Power')
    fig_ax.text(0.05, 1.05, chr(97+alternative_index)+'.)', fontsize=16)#, va='top')
temp_fig.legend(legend_lines, legend_labels, 'center right', title='Sample Size')
temp_fig.suptitle('Wald Test Power', fontsize=16, fontweight='bold')
temp_fig.subplots_adjust(top=.9, bottom=0, right=.85, hspace=.3)
savefig('wald_power.png', bbox_inches='tight')