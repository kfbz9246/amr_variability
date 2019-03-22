'''
Created on Jan 17, 2017

@author: kfbz9246
'''
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
from math import sqrt, log10, isnan
from cStringIO import StringIO
from scipy.stats import chi2
from itertools import combinations

antibiotic_set = ['str','tet','fis']#['amc', 'fox', 'axo', 'tio']#['amc', 'axo', 'fox']#['amp', 'tio', 'amc', 'fox', 'axo', 'chl']
microbe_genus = 'EC'
microbe_species = 'coli'
host = 'Chickens'
start_year = 2005
end_year = 2013


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

engine = create_engine('sqlite:////home/workspace/amr_db/narms_database.sqlite')
session = sessionmaker(bind=engine)()

resistant_query = session.query(Samples.year.label('Year'), Samples.source.label('Source'), Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('1').label('Resistance')).\
                         join(Isolates, Microbes, Tests).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic>=Breakpoints.mic_r, Samples.year>=start_year, Samples.year<=end_year, Microbes.genus==microbe_genus, Samples.host_species==host)
suscptable_query = session.query(Samples.year.label('Year'), Samples.source.label('Source'), Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('0').label('Resistance')).\
                         join(Isolates, Microbes, Tests).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic<Breakpoints.mic_r, Samples.year>=start_year, Samples.year<=end_year, Microbes.genus==microbe_genus, Samples.host_species==host)

if len(antibiotic_set) > 0:
    resistant_query = resistant_query.filter(Tests.drug_id.in_(antibiotic_set))
    suscptable_query = suscptable_query.filter(Tests.drug_id.in_(antibiotic_set))

data_query = resistant_query.union(suscptable_query)
data = read_sql(data_query.statement, engine)

data = pivot_table(data, index=['Source', 'Year', 'IsolateID'], columns='DrugID')


data.columns = [str(col[0]) if str(col[1]) == '' else str(col[1])  for col in data.columns.values]

if len(antibiotic_set) == 0:
    antibiotic_set = tuple(data.columns.values)
    print antibiotic_set

data.dropna(1,thresh=int(len(data)*.2), inplace=True)
data.dropna(how='any',inplace=True)



data.insert(0, 'Combination', None)
for row_name, row_data in data.iterrows():
    combination_string = StringIO()
    for antibiotic in antibiotic_set:
        if row_data[antibiotic] == 0:
            combination_string.write('S')
        else:
            combination_string.write('R')
    data.loc[row_name, 'Combination'] = combination_string.getvalue()
    combination_string.close()

data.insert(0,'Count', 1)
data.reset_index(inplace=True)


data = pivot_table(data, 'Count', ['Source','Year'], 'Combination', sum, 0)

for combination_id in combination_string_list:
    if not combination_id in data.columns:
        data.insert(0, combination_id, 0)
data = data[combination_string_list]
data = data.add(.5)
data.reset_index(inplace=True)

grid_spec = GridSpec(20,3)
color_map = 'Blues'
temp_fig = figure(figsize=(9.25, 7), dpi=300, frameon=True)
temp_fig.suptitle('Cell Proportions')

heat_map_number = 0

for source, source_data in data.groupby('Source'):
    column_totals = source_data.sum(0)
    row_totals = source_data.sum(1)
    column_marginals = column_totals/sum(column_totals)


    for combination_id in combination_string_list:
        source_data[combination_id] /= row_totals

    test_statistic = 0
    for row_name, row_data in source_data.iterrows():
        for combination_id in combination_string_list:
            test_statistic += (((row_data[combination_id]-column_marginals[combination_id])**2)/column_marginals[combination_id])*row_totals[row_name]
    
    ratio_expected_variance = 0
    for column_name, column_value in column_totals.iteritems():
        ratio_expected_variance += 1.0/column_value
    expected = (column_marginals['SSS']*column_marginals['SRR']*column_marginals['RSR']*column_marginals['RRS'])/(column_marginals['SSR']*column_marginals['SRS']*column_marginals['RSS']*column_marginals['RRR'])
    
    ratio_test_statistic = 0
    for row_name, row_data in source_data.iterrows():
        observed= (row_data['SSS']*row_data['SRR']*row_data['RSR']*row_data['RRS'])/(row_data['SSR']*row_data['SRS']*row_data['RSS']*row_data['RRR'])
        ratio_test_statistic += ((observed-expected)**2)*row_totals[row_name]/ratio_test_statistic
    
    
    print source 
    print 1-chi2.cdf(test_statistic, (source_data.shape[0]-2)*(source_data.shape[1]-2))
    print 1-chi2.cdf(ratio_test_statistic, 8)


    source_data=source_data.transpose()





    heat_map_ax = temp_fig.add_subplot(grid_spec[0:16,heat_map_number])
    heat_map_number += 1
    heat_map_ax.imshow(source_data, cmap=color_map, interpolation='nearest', aspect='auto', vmin=0, vmax=1)
    
    heat_map_ax.set_yticks(range(source_data.shape[0]))
    heat_map_ax.set_yticklabels(list(source_data.index))
    heat_map_ax.set_xticks(range(source_data.shape[1]))
    heat_map_ax.set_xticklabels(source_data.columns.values, rotation=45)
    heat_map_ax.set_title(source)
    heat_map_ax.set_ylabel('Years')
    heat_map_ax.set_xlabel('Cell Number')
    
    
    
    
    
    
source_data = pivot_table(data, 'Count', ['Year'], 'Combination', sum, 0)

for combination_id in combination_string_list:
    if not combination_id in source_data.columns:
        source_data.insert(0, combination_id, 0)
source_data = source_data[combination_string_list]
source_data = source_data.add(.5)

column_marginals = source_data.sum(0)
row_totals = source_data.sum(1)
column_marginals/sum(column_marginals)

for combination_id in combination_string_list:
    source_data[combination_id] /= row_totals

test_statistic = 0
for row_name, row_data in source_data.iterrows():
    for combination_id in combination_string_list:
        test_statistic += (((row_data[combination_id]-column_marginals[combination_id])**2)/column_marginals[combination_id])*row_totals[row_name]

source = 'Total'
print source
print 1-chi2.cdf(test_statistic, (source_data.shape[0]-2)*(source_data.shape[1]-2))


source_data=source_data.transpose()





heat_map_ax = temp_fig.add_subplot(grid_spec[0:16,2])
heat_map_ax.imshow(source_data, cmap=color_map, interpolation='nearest', aspect='auto', vmin=0, vmax=1)

heat_map_ax.set_yticks(range(source_data.shape[0]))
heat_map_ax.set_yticklabels(list(source_data.index))
heat_map_ax.set_xticks(range(source_data.shape[1]))
heat_map_ax.set_xticklabels(source_data.columns.values, rotation=45)
heat_map_ax.set_title(source)
heat_map_ax.set_ylabel('Years')
heat_map_ax.set_xlabel('Cell Number')
    
    
    
    
    
    

color_bar_ax = temp_fig.add_subplot(grid_spec[-1,0:3])
color_bar_ax.imshow(array([linspace(0,1,256)]), cmap=color_map, aspect='auto')
color_bar_ax.set_title('Observed Proportion', {'fontsize':12})
color_bar_ax.set_yticks([])
color_bar_ax.set_xticks(linspace(0,256,5))
color_bar_ax.set_xticklabels(linspace(0,1,5))

savefig('heapmap.png', bbox_inches='tight')





















































