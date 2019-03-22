'''
Created on Dec 28, 2016

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_, distinct
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series, merge, DataFrame
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap

drug = 'str'
host ='Chickens'
genus = 'EC'

engine = create_engine('sqlite:////home/code/amr_db/narms_database.sqlite')
session = sessionmaker(bind=engine)()

data_query = session.query(Isolates.id.label('IsolateID'), Samples.year.label('Year'), Samples.state.label('State'), Tests.corrected_mic.label('MIC')).\
                    join(Samples, Microbes, Tests).\
                    filter(Samples.host_species==host, Microbes.genus ==genus, Tests.drug_id==drug, Samples.source=='Retail')
data = read_sql(data_query.statement, engine)

states_list = list(data['State'].unique())
states_list.sort()
states_positions_dict = {}
state_position = len(states_list)/-2.0+.5
for state in states_list:
    states_positions_dict[state]=state_position*.2
    state_position += 1
temp_fig = figure(figsize=(20, 7), dpi=96, frameon=True)
fig_ax = temp_fig.add_subplot(111)
for (state, year), grouped_data in data.groupby(['State', 'Year']):
    fig_ax.violinplot(list(grouped_data['MIC']), [year+states_positions_dict[state]], widths=.2, showmeans=True)
savefig('violine.png', bbox_inches='tight')


