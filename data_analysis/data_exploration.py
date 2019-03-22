'''
Created on Dec 16, 2016

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_, distinct
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series, merge, DataFrame
from numpy import std, array, mean, var, zeros
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap
from math import sqrt, log10, isnan
from numpy import sum
from networkx import Graph, draw_networkx_nodes, draw_networkx_labels, draw_networkx_edges, spring_layout, draw_networkx_edge_labels, draw_networkx, write_graphml
from itertools import combinations, izip, product
from numpy.random import dirichlet
import statsmodels.formula.api as smf
from statsmodels.api import families

data_query_num = -1


engine = create_engine('sqlite:////home/code/amr_db/narms_database.sqlite')
session = sessionmaker(bind=engine)()


if data_query_num == 0:
    data_query = session.query(Isolates.id, Samples.state.label('State'), Samples.year.label('Year'), Samples.host_species.label('Host')).\
                             join(Samples, Microbes).\
                             filter(Samples.source=='Retail', Microbes.genus=='S')#Samples.host_species=='Cattle'


if data_query_num == 1:
    data_query = session.query(Samples.state.label('State'), Samples.year.label('Year'), Samples.host_species.label('Host'), Microbes.genus.label('Genus'), func.count('*').label('Count')).\
                             join(Isolates, Microbes).\
                             filter(Samples.source=='Retail').\
                             group_by(Samples.state.label('State'), Samples.year.label('Year'), Samples.host_species.label('Host'), Microbes.genus.label('Genus'))
    data = read_sql(data_query.statement, engine)
    for key, data_subset in  data.groupby(['Host', 'Genus']):
        print key
        print data_subset.pivot('State','Year', 'Count')

if data_query_num == 2:
    data_query = session.query(Samples.source.label('Source'), Samples.host_species.label('Host'), Microbes.genus.label('Genus'), func.count('*').label('Count')).\
                             join(Isolates, Microbes).\
                             group_by(Samples.source, Samples.host_species, Microbes.genus)
    data = read_sql(data_query.statement, engine)
    data.sort_values('Count', inplace=True)
    data.to_csv('counts.csv')

if data_query_num == 3:
    host='Swine'
    bacteria='EC'
    data_query = session.query(Samples.state.label('State'), Samples.year.label('Year'), func.count('*').label('Count')).\
                             join(Isolates, Microbes).\
                             filter(Samples.source=='Retail', Samples.host_species==host, Microbes.genus==bacteria ).\
                             group_by(Samples.state.label('State'), Samples.year.label('Year'))
    data = read_sql(data_query.statement, engine)
    data=data.pivot('State','Year', 'Count')
    data.to_csv('by_state.csv')

if data_query_num == 4:
    data_query = session.query(distinct(Microbes.serotype))
    data = read_sql(data_query.statement, engine)
    print data
    
if data_query_num == 5:
    data_query = session.query(distinct(Microbes.species)).\
                        filter(Microbes.genus=='E')
    data = read_sql(data_query.statement, engine)
    print data

if data_query_num == 6:
    host='Chickens'
    bacteria='C'
    species='jejuni'#'faecalis'
    data_query = session.query(Samples.state.label('State'), Samples.year.label('Year'), func.count('*').label('Count')).\
                             join(Isolates, Microbes).\
                             filter(Samples.source=='Retail', Samples.host_species==host, Microbes.genus==bacteria, Microbes.species==species).\
                             group_by(Samples.state.label('State'), Samples.year.label('Year'))
    data = read_sql(data_query.statement, engine)
    data=data.pivot('State','Year', 'Count')
    print data
if data_query_num == 7:
    
    data_query = session.query(Samples.state.label('State'), func.count('*').label('Count')).\
                        group_by(Samples.state.label('State'))
    data = read_sql(data_query.statement, engine)
    print data

if data_query_num == 8:
    genus = 'C'
    data_query = session.query(distinct(Tests.drug_id)).\
                        join(Isolates, Microbes).\
                        filter(Microbes.genus==genus)
    data = read_sql(data_query.statement, engine)
    data.to_csv(genus+'.csv')
    print data

if data_query_num == 9:
    data_query = session.query(distinct(Samples.source))
    data = read_sql(data_query.statement, engine)
    print data

if True:
    bacteria='C'
    data_query = session.query(Samples.year.label('Year'), Samples.host_species.label('Host')).distinct().\
                             join(Isolates, Microbes).\
                             filter(Samples.source=='Slaughter', Microbes.genus==bacteria).\
                             order_by(Samples.year)
    data = read_sql(data_query.statement, engine)
    print data