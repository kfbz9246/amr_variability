'''
Created on Jul 17, 2018

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

my_drug = 'fis'#['gen', 'fis', 'tet', 'amp']#['gen', 'fis', 'tet', 'amp']
microbe_genus = 'EC'
microbe_species = 'coli'
host = 'Chickens'
start_year = 2004
end_year = 2012
"""
antibiotic_set_string = ''
antimicrobial_name_lookuper = AntimicrobialNameLookup()
for antibiotic in antibiotic_set:
    antibiotic_set_string += antimicrobial_name_lookuper.lookup(antibiotic) + ' - '
antibiotic_set_string = antibiotic_set_string[0:-3]
"""
engine = create_engine('sqlite:////home/workspace/amr_db/narms_database.sqlite')
session = sessionmaker(bind=engine)()

query = session.query(Samples.year.label('Year'), Samples.source.label('Source'), Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID')).\
                         join(Isolates, Microbes, Tests).\
                         filter(Samples.year==2007, Microbes.genus==microbe_genus, Samples.host_species==host, Tests.drug_id==my_drug, Samples.source=='Slaughter')




data = read_sql(query.statement, engine)
print data