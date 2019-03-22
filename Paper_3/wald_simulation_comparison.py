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
from numpy.random import dirichlet, multinomial, normal
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, savefig


n = 1000
over_dispersion_factor = 4
marginal_probabilites = [.5,.6,.4]
odds_ratio_lists = [[1.2, .7,1.1], [2]]
num_simulations = 500000
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

mean_cell_probabilites = get_expected_cell_probabilites(marginal_probabilites, odds_ratio_lists)
wald_cell_counts = array(mean_cell_probabilites['p'])*n

wald_standard_deviation = 0
for cell_count in wald_cell_counts:
    wald_standard_deviation += 1.0/cell_count
wald_standard_deviation = sqrt(wald_standard_deviation)
#overdispersed_wald_standard_deviation = sqrt(over_dispersion_factor) * wald_standard_deviation
overdispersed_wald_standard_deviation = over_dispersion_factor * wald_standard_deviation
wald_observations = normal(log(odds_ratio_lists[-1][0]), overdispersed_wald_standard_deviation, num_simulations)


simulated_log_odds_ratios = []
alpha = (n-over_dispersion_factor)/float((over_dispersion_factor*dimention)-dimention)
direchlet_parameters = array(mean_cell_probabilites['p'])*alpha
cell_index_is_even = mean_cell_probabilites['even']
for i in range(num_simulations):
    p = dirichlet(direchlet_parameters)
    sample = multinomial(n, p)+.5
    odds_ratio = 1.0
    for is_even, sample_value in zip(cell_index_is_even, sample):
        if is_even:
            odds_ratio *= sample_value
        else:
            odds_ratio /= sample_value

    simulated_log_odds_ratios.append(log(odds_ratio))

hist_min = min([min(wald_observations), min(simulated_log_odds_ratios)])
hist_max = max([max(wald_observations), max(simulated_log_odds_ratios)])
hist_range = (hist_min, hist_max)
temp_fig = figure(figsize=(9.25, 7), dpi=300, frameon=True)
fig_ax = temp_fig.add_subplot(111)
fig_ax.hist(wald_observations, 500, hist_range, color='g', alpha=.5)
fig_ax.hist(simulated_log_odds_ratios, 500, hist_range, color='b', alpha=.5)
#box = fig_ax.get_position()
#fig_ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#fig_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Sample_size')
#fig_ax.set_title('3-Way Interaction Parameter=%.3f' % alternative_odds_ratio_lists[-1][0])
fig_ax.set_xlabel('Overdispersion Factor')
fig_ax.set_ylabel('Power')

savefig('power_hists.png', bbox_inches='tight')