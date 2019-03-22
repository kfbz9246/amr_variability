'''
Created on Feb 6, 2018

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased

import sqlalchemy
print sqlalchemy.__version__
exit()

engine = create_engine('sqlite:////home/workspace/amr_db/narms_database.sqlite')
session = sessionmaker(bind=engine)()

resistant_query = session.query(Isolates.microbe_id.label('isolates_microbe_id'), Microbes.id.label('microbe_id'))

from pandas import read_sql
print read_sql(resistant_query.statement, engine)
