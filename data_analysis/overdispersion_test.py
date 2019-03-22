'''
Created on Feb 15, 2017

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series, merge, DataFrame, pivot_table, read_csv
from numpy import std, array, mean, var, zeros, sum, linspace
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, savefig
from matplotlib.gridspec import GridSpec
from math import sqrt, log10, isnan
from cStringIO import StringIO
from scipy.stats import chi2
from itertools import combinations

data = read_csv('test.csv')
combination_string_list = data.columns.values

def run_analysis(data_table, source):
    print source
    data_table.to_csv(source+'.csv')
    column_totals = data_table.sum(0)
    row_totals = data_table.sum(1)
    column_marginals = column_totals/sum(column_totals)
    for combination_id in combination_string_list:
        data_table[combination_id] /= row_totals

    data_table = data_table.transpose()
    heat_map_fig = figure(figsize=(9.25, 7), dpi=300, frameon=True)
    heat_map_ax = heat_map_fig.add_subplot(111)
    heat_map_ax.imshow(data_table, cmap='Blues', interpolation='nearest', aspect='auto', vmin=0, vmax=1)
    heat_map_ax.set_yticks(range(data_table.shape[0]))
    heat_map_ax.set_yticklabels(list(data_table.index))
    heat_map_ax.set_xticks(range(data_table.shape[1]))
    heat_map_ax.set_xticklabels(data_table.columns.values, rotation=45)
    heat_map_ax.set_title(source)
    heat_map_ax.set_ylabel('Resistance Pattern')
    heat_map_ax.set_xlabel('Years')
    savefig(source+'_heatmap.png', bbox_inches='tight')
    data_table = data_table.transpose()

    multinomial_over_dispersion_test_statistic = 0
    for row_name, row_data in data_table.iterrows():
        for combination_id in combination_string_list:
            multinomial_over_dispersion_test_statistic += (((row_data[combination_id]-column_marginals[combination_id])**2)/column_marginals[combination_id])*row_totals[row_name]
    print 'multinomial test: ' +  str(1-chi2.cdf(multinomial_over_dispersion_test_statistic, (data_table.shape[0]-1)*(data_table.shape[1]-1)))

    odds_ratio_expected_variance = 0
    for column_name, column_total in column_totals.iteritems():
        odds_ratio_expected_variance += 1.0/column_total
    
    overall_odds_ratio = 1
    for column_name, marginal_column_proportion in column_marginals.iteritems():
        if column_name.count('R')%2 == 0:
            overall_odds_ratio *= marginal_column_proportion
        else:
            overall_odds_ratio /= marginal_column_proportion
    overall_odds_ratio = log10(overall_odds_ratio)
            
    odds_ratio_figure = figure(figsize=(9.25, 7), dpi=300, frameon=True)
    odds_ratio_scatter_plot_axis = odds_ratio_figure.add_subplot(111)
    ratio_test_statistic = 0
    for year, year_data in data_table.iterrows():
        year_odds_ratio = 1
        for column_name, column_value in year_data.iteritems():
            if column_name.count('R')%2 == 0:
                year_odds_ratio *= column_value
            else:
                year_odds_ratio /= column_value
        year_odds_ratio = log10(year_odds_ratio)
                
        ratio_test_statistic += (((year_odds_ratio-overall_odds_ratio)**2)*row_totals[year])/odds_ratio_expected_variance
        odds_ratio_scatter_plot_axis.errorbar([year], [year_odds_ratio], yerr=sqrt(odds_ratio_expected_variance), fmt='.-', color='k', linewidth=2)
        odds_ratio_scatter_plot_axis.scatter([year], [year_odds_ratio], color='b')
    min_year = min(data_table.index)
    max_year = max(data_table.index)
    odds_ratio_scatter_plot_axis.plot([min_year-1, max_year +1], [overall_odds_ratio, overall_odds_ratio], color='k', linestyle='--', linewidth=2)
    odds_ratio_scatter_plot_axis.text(min_year -.9,overall_odds_ratio+.1, 'Average Odds Ratio')
    odds_ratio_scatter_plot_axis.set_xlim([min_year-1, max_year +1])
    odds_ratio_scatter_plot_axis.set_xlabel('Year')
    odds_ratio_scatter_plot_axis.set_ylabel('Difference in Log Odds Ratio')
    odds_ratio_scatter_plot_axis.set_title(source)
    savefig(source+'_odds_ratio_plot.png', bbox_inches='tight')
    print 'odds-ratio test: ' + str(1-chi2.cdf(ratio_test_statistic, (data_table.shape[0]-1)))

run_analysis(data, 'total')








































