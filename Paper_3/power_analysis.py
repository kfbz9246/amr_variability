'''
Created on Apr 22, 2017

@author: kfbz9246
'''

from math import log, sqrt
from copy import deepcopy
from itertools import product, combinations
from pandas import DataFrame, Panel
from numpy import array, complex128, mean, sum, inner
from numpy.polynomial.polynomial import polyfromroots, polyroots
from numpy.random import dirichlet, multinomial
from scipy.stats import norm, uniform
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap, tight_layout
from matplotlib.gridspec import GridSpec
from patsy import dmatrix
from numpy import matrix, array
from math import exp

null_coeffs = [-1.836, -.4055, .04886,-.8598,.4055,.4055,.4055, 0]
over_dispersion_factors = [1,2,3,4,5,6,7,8,9,10]#[3,10]#
interaction_parameters = [.5,4]#[0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]#
sample_sizes = [250, 500, 2500, 7500]
dimention = 3

number_simulations = 50000




variable_list = []
for variable_number in range(dimention):
    variable_name = 'X' + str(variable_number)
    variable_list.append(variable_name)

base_cell_data = []                           
for a_product in product([0,1], repeat=dimention):
    next_base_cell_data_row_dict = {}
    for variable_name, value in zip(variable_list, a_product):
        next_base_cell_data_row_dict[variable_name] = value
    base_cell_data.append(next_base_cell_data_row_dict)
base_cell_data = DataFrame(base_cell_data)

for row_name, row in base_cell_data.iterrows():
    if sum(row[range(dimention)])%2==0:
        base_cell_data.loc[row_name, 'even'] = True
    else:
        base_cell_data.loc[row_name, 'even'] = False
        
design_matrix = matrix(dmatrix('('+'+'.join(variable_list)+')**'+str(dimention), base_cell_data))
design_matrix_column_names = []
for design_matrix_column_id in range(design_matrix.shape[1]):
    column_name = 'D'+str(design_matrix_column_id)
    design_matrix_column_names.append(column_name)
    base_cell_data[column_name] = design_matrix[:, design_matrix_column_id]

simulated_powers = {}
wald_powers = {}
for sample_size in sample_sizes:
    simulated_powers_at_sample_size = {}
    wald_powers_at_sample_size = {}
    for over_dispersion_factor in over_dispersion_factors:
        simulated_powers_at_sample_size_overdispersion = {}
        wald_powers_at_sample_size_overdispersion = {}
        coeffs = deepcopy(null_coeffs)
        cell_data = base_cell_data['even'].to_frame()
        for row_id, row_data in base_cell_data.iterrows():
            cell_data.loc[row_id, 'p'] =  exp(inner(coeffs, row_data.loc[design_matrix_column_names])) 
        sum_p = sum(cell_data['p'])
        cell_data['p'] = cell_data['p']/sum_p


        cell_counts = (array(cell_data['p'])*sample_size)+.5
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
        
        lower_wald_critical_value = norm.ppf(.025, mean, standard_deviation)
        upper_wald_critical_value = norm.ppf(.975, mean, standard_deviation)
        
    
        simulated_null_distrobution = []
        if over_dispersion_factor == 1:
            for i in range(number_simulations):
                sample = multinomial(sample_size, cell_data['p'])+.5
                
                odds_ratio = 1.0
                for sample_value, (cell_id, cell_is_even) in zip(sample, cell_data['even'].iteritems()):
                    if cell_is_even:
                        odds_ratio *= sample_value
                    else:
                        odds_ratio /= sample_value
                simulated_null_distrobution.append(log(odds_ratio))
    
        else:
            alpha = (sample_size-over_dispersion_factor)/float(over_dispersion_factor-1)
            direchlet_parameters = array(cell_data['p'])*alpha
        
            for i in range(number_simulations):
                p = dirichlet(direchlet_parameters)
                sample = multinomial(sample_size, p)+.5
                odds_ratio = 1.0
                for sample_value, (cell_id, cell_is_even) in zip(sample, cell_data['even'].iteritems()):
                    if cell_is_even:
                        odds_ratio *= sample_value
                    else:
                        odds_ratio /= sample_value
            
                simulated_null_distrobution.append(log(odds_ratio))
        simulated_null_distrobution.sort()
        lower_simulated_critical_value = simulated_null_distrobution[int(number_simulations*.025)]
        upper_simulated_critical_value = simulated_null_distrobution[int(number_simulations*.975)]
        
        
        for interaction_parameter in interaction_parameters:
            coeffs = deepcopy(null_coeffs)
            coeffs[-1] = interaction_parameter
            cell_data = base_cell_data['even'].to_frame()
            for row_id, row_data in base_cell_data.iterrows():
                cell_data.loc[row_id, 'p'] =  exp(inner(coeffs, row_data.loc[design_matrix_column_names])) 
            sum_p = sum(cell_data['p'])
            cell_data['p'] = cell_data['p']/sum_p
            
            cell_counts = (array(cell_data['p'])*sample_size)+.5
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

            wald_powers_at_sample_size_overdispersion[interaction_parameter] = norm.cdf(lower_wald_critical_value, mean, standard_deviation)+ 1 -norm.cdf(upper_wald_critical_value, mean, standard_deviation)
            

            samples_rejected_count = 0.0
            if over_dispersion_factor == 1:
                for i in range(number_simulations):
                    sample = multinomial(sample_size, cell_data['p'])+.5
                    
                    odds_ratio = 1.0
                    for (cell_id, cell_is_even), sample_value in zip(cell_data['even'].iteritems(), sample):
                        if cell_is_even:
                            odds_ratio *= sample_value
                        else:
                            odds_ratio /= sample_value
                    odds_ratio = log(odds_ratio) 
                    if odds_ratio < lower_simulated_critical_value  or  upper_simulated_critical_value < odds_ratio:
                        samples_rejected_count += 1 
        
            else:
                alpha = (sample_size-over_dispersion_factor)/float(over_dispersion_factor-1)
                direchlet_parameters = array(cell_data['p'])*alpha
            
                for i in range(number_simulations):
                    p = dirichlet(direchlet_parameters)
                    sample = multinomial(sample_size, p)+.5
                    odds_ratio = 1.0
                    for (row_id, row_is_even), sample_value in zip(cell_data['even'].iteritems(), sample):
                        if row_is_even:
                            odds_ratio *= sample_value
                        else:
                            odds_ratio /= sample_value
                
                    odds_ratio = log(odds_ratio) 
                    if odds_ratio < lower_simulated_critical_value  or  upper_simulated_critical_value < odds_ratio:
                        samples_rejected_count += 1 
            simulated_powers_at_sample_size_overdispersion[interaction_parameter]=samples_rejected_count/number_simulations
    
        wald_powers_at_sample_size[over_dispersion_factor] = wald_powers_at_sample_size_overdispersion
        simulated_powers_at_sample_size[over_dispersion_factor] = simulated_powers_at_sample_size_overdispersion
    wald_powers[sample_size] = DataFrame(wald_powers_at_sample_size)
    simulated_powers[sample_size] = DataFrame(simulated_powers_at_sample_size)
wald_powers = Panel(wald_powers)
simulated_powers = Panel(simulated_powers)


num_subplots = None
powers_generator = None
title_base = None
x_label = ''

if len(over_dispersion_factors) < len(interaction_parameters):
    num_subplots = len(over_dispersion_factors)
    powers_generator = ((over_dispersion_factor, wald_powers.loc[:, :, over_dispersion_factor].transpose(), simulated_powers.loc[:, :, over_dispersion_factor].transpose()) for over_dispersion_factor in over_dispersion_factors)
    title_base = 'Dispersion = '
    x_label = 'Model Parameter'
else:
    num_subplots = len(interaction_parameters)
    powers_generator = ((interaction_parameter, wald_powers.loc[:, interaction_parameter, :].transpose(), simulated_powers.loc[:, interaction_parameter, :].transpose()) for interaction_parameter in interaction_parameters)
    title_base = 'Model Parameter = '
    x_label = 'Dispersion'

grid_spec = GridSpec(num_subplots, 1)
marker_list = ['o', 'v', 's', 'P','s', 'h']
temp_fig = figure(figsize=(9.25, 7), dpi=300, frameon=True)#, tight_layout=True)
lines_for_legend = None
labels_for_legend = None
for plot_number, power_data in enumerate(powers_generator):
    fig_ax = temp_fig.add_subplot(grid_spec[plot_number, 0])
    outer_param_val = power_data[0]
    wald_data = power_data[1]
    simulated_data = power_data[2]
    lines_for_legend = []
    labels_for_legend = []
    
    for marker_type, (sample_size, power_list) in zip(marker_list, wald_data.iterrows()):
        plot_lines = fig_ax.plot(power_list.index.values, power_list, c='k', marker=marker_type)
        lines_for_legend.append(plot_lines[0])
        labels_for_legend.append(str(sample_size)+' Theoretical')
    legend_location = 1
    for marker_type, (sample_size, power_list) in zip(marker_list, simulated_data.iterrows()):
        plot_lines = fig_ax.plot(power_list.index.values, power_list, c='k', marker=marker_type, linestyle='--')
        lines_for_legend.insert(legend_location, plot_lines[0])
        labels_for_legend.insert(legend_location, str(sample_size)+' Simulated')
        legend_location += 2
    fig_ax.set_title(title_base+str(outer_param_val))
    fig_ax.set_ylabel('Power')
    box = fig_ax.get_position()
    fig_ax.set_position([box.x0, box.y0, box.width * 0.85, box.height*.9])
    #fig_ax.tick_params(axis='both', which='major', labelsize=12)
    fig_ax.set_xlabel(x_label)
    fig_ax.set_ylim(0,1.05)
temp_fig.suptitle('Power Analysis', y=.925, fontsize=16)
temp_fig.legend(lines_for_legend, labels_for_legend, loc='center right')
temp_fig.savefig('power_plots.png', bbox_inches='tight')
    




#box = fig_ax.get_position()
#fig_ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#fig_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)#, title='Sample Size')
    
