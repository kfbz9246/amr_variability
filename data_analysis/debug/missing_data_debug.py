'''
Created on Apr 5, 2017

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

antibiotic_set = ['gen', 'fis', 'tet', 'amp']#['gen', 'fis', 'tet', 'amp']
microbe_genus = 'EC'
microbe_species = 'coli'
host = 'Chickens'
start_year = 2005
end_year = 2013

antibiotic_set_string = ''
antimicrobial_name_lookuper = AntimicrobialNameLookup()
for antibiotic in antibiotic_set:
    antibiotic_set_string += antimicrobial_name_lookuper.lookup(antibiotic) + ' - '
antibiotic_set_string = antibiotic_set_string[0:-3]
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
data.reset_index(inplace=True)
print data.groupby(['Year', 'Source']).count()
