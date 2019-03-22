'''
Created on Jul 15, 2018

@author: kfbz9246
'''
from antimicrobail_name_lookup import AntimicrobialNameLookup
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series, merge, DataFrame, pivot_table, read_csv
from numpy import std, array, mean, var, zeros, sum, linspace
from numpy.random import multinomial
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, savefig
from matplotlib.gridspec import GridSpec
from math import sqrt, log, isnan
from cStringIO import StringIO
from scipy.stats import norm
from itertools import combinations

overdispersion_table = read_csv('overdispersions.csv', index_col=0)
full_table_overdispersion_table = read_csv('overall_dispersion.csv', index_col=0)


num_simulations = 1000
base_antibiotic_set = ['gen', 'fis', 'tet', 'amp']#['gen', 'fis', 'tet', 'amp']
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

data = pivot_table(data, index=['Year', 'IsolateID', 'Source'], columns='DrugID')


data.columns = [str(col[0]) if str(col[1]) == '' else str(col[1])  for col in data.columns.values]
data.dropna(how='any',inplace=True)



number_observations = data.reset_index()[['Source', 'Year']]
number_observations.insert(0, 'Count', 1)
number_observations = number_observations.groupby(['Source', 'Year']).count()

data.index = data.index.droplevel(2)


p_value_table = {}

def run_analysis(data_table, antibiotic_set_string):
    print antibiotic_set_string
    row_totals = data_table.sum(1)
    data_table = data_table.divide(row_totals,0)

    overdispersion_factor_ratios = []
    full_table_overdispersion_value_ratios = []
    for simulation_index in range(num_simulations):
        overdispersion_factors = []
        full_table_overdispersion_values = []
        for source in ['Slaughter', 'Retail']:
            bootstrap_data = {}
            for year, data_row in data_table.iterrows():
                bootstrap_data[year] =  dict(zip(data_row.index, multinomial(number_observations.loc[(source, year), 'Count'], data_row)))
            bootstrap_data =  DataFrame(bootstrap_data).transpose()[data_table.columns.values]
            
            

            
            bootstrap_data = bootstrap_data[list(data_table.columns.values)]

            
            bootstrap_data = bootstrap_data.add(.5)
            column_totals = bootstrap_data.sum(0)
            row_totals = bootstrap_data.sum(1)
            total_observations = sum(column_totals)
            column_marginals = column_totals/total_observations
            num_years = bootstrap_data.shape[0]
            num_cells = bootstrap_data.shape[1]
            
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
            for (year, year_counts), (yearagain, year_total) in zip(bootstrap_data.iterrows(), row_totals.iteritems()):
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
                    
    
            overdispersion_factors.append(ratio_test_statistic/(num_years-1))
            full_table_overdispersion_values.append(full_table_overdispersion_test_statistic/((num_years-1)*(num_cells-1)))
    

        overdispersion_factor_ratios.append(log(overdispersion_factors[0]/overdispersion_factors[1]))
        full_table_overdispersion_value_ratios.append(log(full_table_overdispersion_values[0]/full_table_overdispersion_values[1]))
    
    fig = figure(figsize=(9.25, 7), dpi=300, frameon=True)
    a_grispec = GridSpec(2,1)
    fig_ax1 = fig.add_subplot(a_grispec[0,0])
    fig_ax1.hist(overdispersion_factor_ratios, 50)
    fig_ax2 = fig.add_subplot(a_grispec[1,0])
    fig_ax2.hist(full_table_overdispersion_value_ratios, 50)
    fig.savefig(antibiotic_set_string+'.png', bbox_inches='tight')
    
    #overdispersion_factor_ratios_stdev = std(overdispersion_factor_ratios)
    observed_overdispersion_factor_ratio = abs(log(overdispersion_table.loc[antibiotic_set_string, 'Slaughter']/overdispersion_table.loc[antibiotic_set_string, 'Retail']))
    observed_overdispersion_factor_ratio *= -1
    overdispersion_factor_ratios.sort()
    overdispersion_factor_ratios_index = 0
    while overdispersion_factor_ratios[overdispersion_factor_ratios_index] < observed_overdispersion_factor_ratio:
        overdispersion_factor_ratios_index += 1
    overdispersion_factor_ratios_p_value = overdispersion_factor_ratios_index/float(num_simulations)
    observed_overdispersion_factor_ratio *= -1
    while overdispersion_factor_ratios[overdispersion_factor_ratios_index] <= observed_overdispersion_factor_ratio:
        overdispersion_factor_ratios_index += 1
    overdispersion_factor_ratios_p_value += (num_simulations - overdispersion_factor_ratios_index)/float(num_simulations)
    p_value_table[antibiotic_set_string]['model_term_overdispersion'] = overdispersion_factor_ratios_p_value
    #p_value_table[antibiotic_set_string]['model_term_overdispersion'] = 2*(1-norm.cdf(observed_overdispersion_factor_ratio, scale = overdispersion_factor_ratios_stdev))

    #full_table_overdispersion_value_ratios_stdev = std(full_table_overdispersion_value_ratios)
    observed_full_table_overdispersion_value_ratio = abs(log(full_table_overdispersion_table.loc[antibiotic_set_string, 'Slaughter']/full_table_overdispersion_table.loc[antibiotic_set_string, 'Retail']))
    observed_full_table_overdispersion_value_ratio *= -1
    full_table_overdispersion_value_ratios.sort()
    full_table_overdispersion_value_ratios_index = 0
    while full_table_overdispersion_value_ratios[full_table_overdispersion_value_ratios_index] < observed_full_table_overdispersion_value_ratio:
        full_table_overdispersion_value_ratios_index += 1
    full_table_overdispersion_value_ratios_p_value = full_table_overdispersion_value_ratios_index/float(num_simulations)
    observed_full_table_overdispersion_value_ratio *= -1
    while full_table_overdispersion_value_ratios[full_table_overdispersion_value_ratios_index] <= observed_full_table_overdispersion_value_ratio:
        full_table_overdispersion_value_ratios_index += 1
    full_table_overdispersion_value_ratios_p_value += (num_simulations - full_table_overdispersion_value_ratios_index)/float(num_simulations)
    p_value_table[antibiotic_set_string]['multinomail_overdispersion'] = full_table_overdispersion_value_ratios_p_value
    #p_value_table[antibiotic_set_string]['multinomail_overdispersion'] = 2*(1-norm.cdf(observed_full_table_overdispersion_value_ratio, scale = full_table_overdispersion_value_ratios_stdev))

    



for antibiotic_set_size in range(1, len(base_antibiotic_set)+1):
    for antibiotic_set in combinations(base_antibiotic_set, antibiotic_set_size):
        antibiotic_set_string = ''
        antimicrobial_name_lookuper = AntimicrobialNameLookup()
        for antibiotic in antibiotic_set:
            antibiotic_set_string += antimicrobial_name_lookuper.lookup(antibiotic) + ' - '
        antibiotic_set_string = antibiotic_set_string[0:-3]

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
        

        antibiotic_set_data = pivot_table(antibiotic_set_data, 'Count', 'Year', 'Combination', sum, 0)


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

        
        

        run_analysis(antibiotic_set_data, antibiotic_set_string)
DataFrame(p_value_table).transpose().to_csv('overdispersion_difference_pvalues.csv')

