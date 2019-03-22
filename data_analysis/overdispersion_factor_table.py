'''
Created on Apr 10, 2018

@author: kfbz9246
'''
from antimicrobail_name_lookup import AntimicrobialNameLookup
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series, merge, DataFrame, pivot_table
from numpy import std, array, mean, var, zeros, sum, linspace
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, savefig
from matplotlib.gridspec import GridSpec
from math import sqrt, log, isnan
from cStringIO import StringIO
from scipy.stats import chi2
from itertools import combinations

base_antibiotic_set = ['amp', 'gen', 'fis', 'tet']#['gen', 'fis', 'tet', 'amp']
microbe_genus = 'EC'
microbe_species = 'coli'
host = 'Chickens'
start_year = 2004
end_year = 2012
exclude_years = (2007,)

query_years = []
for year in range(start_year, end_year+1):
    if year not in exclude_years:
        query_years.append(year)



engine = create_engine('sqlite:///narms_database.sqlite')
session = sessionmaker(bind=engine)()

resistant_query = session.query(Samples.year.label('Year'), Samples.source.label('Source'), Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('1').label('Resistance')).\
                         join(Isolates, Microbes, Tests).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic>=Breakpoints.mic_r, Samples.year.in_(query_years), Microbes.genus==microbe_genus, Samples.host_species==host)
suscptable_query = session.query(Samples.year.label('Year'), Samples.source.label('Source'), Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('0').label('Resistance')).\
                         join(Isolates, Microbes, Tests).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic<Breakpoints.mic_r, Samples.year.in_(query_years), Microbes.genus==microbe_genus, Samples.host_species==host)

if len(base_antibiotic_set) > 0:
    resistant_query = resistant_query.filter(Tests.drug_id.in_(base_antibiotic_set))
    suscptable_query = suscptable_query.filter(Tests.drug_id.in_(base_antibiotic_set))

data_query = resistant_query.union(suscptable_query)
data = read_sql(data_query.statement, engine)

data = pivot_table(data, index=['Source', 'Year', 'IsolateID'], columns='DrugID')
data.columns = [str(col[0]) if str(col[1]) == '' else str(col[1])  for col in data.columns.values]
data.dropna(how='any',inplace=True)



overall_odds_ratio_table = {}
full_table_overdispersion_value = {}
full_table_overdispersion_p_value = {}
overdispersion_table = {}
p_value_table = {}

def run_analysis(data_table, antibiotic_set_string, source, odds_ratio_scatter_plot_axis):
    
    data_table = data_table.add(.5)
    column_totals = data_table.sum(0)
    row_totals = data_table.sum(1)
    total_observations = sum(column_totals)
    column_marginals = column_totals/total_observations
    num_years = data_table.shape[0]
    num_cells = data_table.shape[1]
    
    odds_ratio_expected_variance_delta_part = 0
    for column_name, marginal in column_marginals.iteritems():
        odds_ratio_expected_variance_delta_part += 1.0/marginal
    
    overall_odds_ratio = 1.0
    for column_name, column_total in column_totals.iteritems():
        if column_name.count('R')%2 == 0:
            overall_odds_ratio *= column_total
        else:
            overall_odds_ratio /= column_total
    overall_odds_ratio = log(overall_odds_ratio)
    


         
    ratio_test_statistic = 0
    full_table_overdispersion_test_statistic = 0
    for (year, year_counts), (yearagain, year_total) in zip(data_table.iterrows(), row_totals.iteritems()):
        estimated_year_cell_frequencies = year_counts/year_total
        cell_frequency_differences = estimated_year_cell_frequencies - column_marginals
        full_table_overdispersion_test_statistic +=  sum((year_total * (cell_frequency_differences**2))/column_marginals)
    
        year_odds_ratio = 1
        for column_name, column_value in year_counts.iteritems():
            if column_name.count('R')%2 == 0:
                year_odds_ratio *= column_value
            else:
                year_odds_ratio /= column_value
        year_odds_ratio = log(year_odds_ratio)
        ratio_test_statistic += (((year_odds_ratio-overall_odds_ratio)**2)*year_total)/odds_ratio_expected_variance_delta_part
        
        year_odds_ratio_expected_variance = 0
        for cell_name, cell_data in year_counts.iteritems():
            year_odds_ratio_expected_variance += 1.0/cell_data
                
        odds_ratio_scatter_plot_axis.errorbar([year], [year_odds_ratio], yerr=(2*sqrt(year_odds_ratio_expected_variance)), fmt='o', capsize=7, color='k', linewidth=2)
    

    
    overall_odds_ratio_table[antibiotic_set_string][source] = overall_odds_ratio
    full_table_overdispersion_value[antibiotic_set_string][source] = full_table_overdispersion_test_statistic/((num_years-1)*(num_cells-1))
    full_table_overdispersion_p_value[antibiotic_set_string][source] =1-chi2.cdf(full_table_overdispersion_test_statistic, (num_years-1)*(num_cells-1))
    over_dispersion_p_value = 1-chi2.cdf(ratio_test_statistic, (num_years-1))
    p_value_table[antibiotic_set_string][source] = over_dispersion_p_value
    overdispersion_factor = ratio_test_statistic/(num_years-1)
    overdispersion_table[antibiotic_set_string][source] = overdispersion_factor
    

    
    

    min_year = min(data_table.index)
    max_year = max(data_table.index)
    odds_ratio_scatter_plot_axis.plot([min_year-1, max_year +1], [overall_odds_ratio, overall_odds_ratio], color='k', linestyle='--', linewidth=2)
    y_lim = odds_ratio_scatter_plot_axis.get_ylim()
    top_of_plot_coord = y_lim[1]
    odds_ratio_scatter_plot_axis.set_ylim(y_lim[0], top_of_plot_coord+(y_lim[1] - y_lim[0])*.1)
    odds_ratio_scatter_plot_axis.text(min_year -.9,overall_odds_ratio+(.025*(y_lim[1]-y_lim[0])), 'Average')
    if over_dispersion_p_value >= .00001:
        odds_ratio_scatter_plot_axis.text(max_year-3, top_of_plot_coord, 'Overdispersion Test p-Value %.3g'%(over_dispersion_p_value))
    else:
        odds_ratio_scatter_plot_axis.text(max_year-3, top_of_plot_coord, 'Overdispersion Test p-Value <.00001')
    odds_ratio_scatter_plot_axis.set_xlim([min_year-1, max_year +1])
    odds_ratio_scatter_plot_axis.set_ylabel('Highest Order\nModel Parameter')
    odds_ratio_scatter_plot_axis.set_title(source)


for antibiotic_set_size in range(1, len(base_antibiotic_set)+1):
    for antibiotic_set in combinations(base_antibiotic_set, antibiotic_set_size):
        antibiotic_set_string = ''
        antimicrobial_name_lookuper = AntimicrobialNameLookup()
        for antibiotic in antibiotic_set:
            antibiotic_set_string += antimicrobial_name_lookuper.lookup(antibiotic) + ' - '
        antibiotic_set_string = antibiotic_set_string[0:-3]
        
        overall_odds_ratio_table[antibiotic_set_string] = {}
        full_table_overdispersion_value[antibiotic_set_string] = {}
        full_table_overdispersion_p_value[antibiotic_set_string] = {}
        overdispersion_table[antibiotic_set_string]= {}
        p_value_table[antibiotic_set_string] = {}
        
        antibiotic_set_data = data.loc[:, antibiotic_set]

        antibiotic_set_data.insert(0, 'Combination', None)
        for row_name, row_data in antibiotic_set_data.iterrows():
            combination_string = StringIO()
            for antibiotic in antibiotic_set:
                if row_data[antibiotic] == 0:
                    combination_string.write('S')
                else:
                    combination_string.write('R')
            antibiotic_set_data.loc[row_name, 'Combination'] = combination_string.getvalue()
            combination_string.close()
        
        antibiotic_set_data.insert(0,'Count', 1)
        antibiotic_set_data.reset_index(inplace=True)
        
        
        antibiotic_set_data = pivot_table(antibiotic_set_data, 'Count', ['Source','Year'], 'Combination', sum, 0)


        combination_string_list  = []
        for number_resistant in range(len(antibiotic_set)+1):
            for resistant_combination in combinations(range(len(antibiotic_set)), number_resistant):
                combination_string = StringIO()
                for combination_string_index in range(len(antibiotic_set)):
                    if combination_string_index in resistant_combination:
                        combination_string.write('R')
                    else:
                        combination_string.write('S')
                combination_string_list.append(combination_string.getvalue())
                combination_string.close()
        
        for combination_id in combination_string_list:
            if not combination_id in antibiotic_set_data.columns:
                antibiotic_set_data.insert(0, combination_id, 0)
        antibiotic_set_data = antibiotic_set_data[combination_string_list]
        antibiotic_set_data.reset_index(inplace=True)

        odds_ratio_grispec = GridSpec(3,1)
        odds_ratio_figure = figure(figsize=(9.25, 7), dpi=300, frameon=True)
        odds_ratio_figure.suptitle(antibiotic_set_string)
        odds_ratio_subplot_counter = 0
        prev_axis = None
        data_grouped_by_source = antibiotic_set_data.groupby('Source')#, squeeze=True)
        for source in ['Slaughter', 'Retail']:
            source_data = data_grouped_by_source.get_group(source)
            source_data = source_data.drop('Source', 1)
            source_data.set_index('Year', inplace=True)
            if prev_axis == None:
                odds_ratio_scatter_plot_axis = odds_ratio_figure.add_subplot(odds_ratio_grispec[odds_ratio_subplot_counter, 0])
            else:
                odds_ratio_scatter_plot_axis = odds_ratio_figure.add_subplot(odds_ratio_grispec[odds_ratio_subplot_counter, 0], sharex=prev_axis)
            run_analysis(source_data, antibiotic_set_string, source, odds_ratio_scatter_plot_axis)
            [label.set_visible(False) for label in odds_ratio_scatter_plot_axis.get_xticklabels()]
            prev_axis = odds_ratio_scatter_plot_axis
            odds_ratio_subplot_counter += 1
        antibiotic_set_data = antibiotic_set_data.groupby('Year').sum()
        odds_ratio_scatter_plot_axis = odds_ratio_figure.add_subplot(odds_ratio_grispec[odds_ratio_subplot_counter, 0], sharex=prev_axis)
        odds_ratio_scatter_plot_axis.set_xlabel('Year')

        run_analysis(antibiotic_set_data, antibiotic_set_string, 'Total', odds_ratio_scatter_plot_axis)
        odds_ratio_figure.savefig(antibiotic_set_string+'.png', bbox_inches='tight')
DataFrame(overdispersion_table).transpose().to_csv('overdispersions.csv')
DataFrame(p_value_table).transpose().to_csv('pvlaues.csv')
DataFrame(full_table_overdispersion_value).transpose().to_csv('overall_dispersion.csv')
DataFrame(full_table_overdispersion_p_value).transpose().to_csv('full_table_p_values.csv')
DataFrame(overall_odds_ratio_table).transpose().to_csv('odds_ratio_table.csv')

data = data.reset_index().groupby(['Year', 'Source']).size().unstack()
data['Total'] = data['Retail'] + data['Slaughter']
data.to_csv('counts.csv')
