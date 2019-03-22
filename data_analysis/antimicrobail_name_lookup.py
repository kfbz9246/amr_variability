'''
Created on Mar 30, 2017

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased

engine = create_engine('sqlite:///narms_database.sqlite')
session = sessionmaker(bind=engine)()

class AntimicrobialNameLookup(object):
    def __init__(self):
        self.lookup_table = {}
        engine = create_engine('sqlite:///narms_database.sqlite')
        session = sessionmaker(bind=engine)()
        for drug_id, drug_name in session.query(Drugs.id, Drugs.name).all():
            self.lookup_table[drug_id] = drug_name
        self.lookup_table['gen'] = 'Gentamicin'
    def lookup(self, drug_id):
        return self.lookup_table[drug_id]