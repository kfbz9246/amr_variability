'''
Created on Sep 24, 2015

@author: kfbz9246
'''
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float, Boolean, ForeignKey
from sqlalchemy.orm import relationship, backref
Base = declarative_base()


class Samples(Base):
    __tablename__ = 'Samples'
    id = Column('ID', String, primary_key=True)
    year = Column('Year', Integer)
    state = Column('State', String)
    meat_type = Column('MeatType', String)
    organic = Column('Organic', Boolean)
    agency = Column('Agency', String)
    source = Column('Source', String)
    host_species = Column('HostSpecies', String)

class Microbes(Base):    
    __tablename__ = 'Microbes'
    id = Column('ID', Integer, primary_key=True)
    genus = Column('Genus', String)
    species = Column('Species', String)
    serotype = Column('Serotype', String)
    antigenic_formula = Column('AntigenicFormula', String)
    
class Isolates(Base):
    __tablename__ = 'Isolates'
    id = Column('ID', String, primary_key=True)
    sample_id = Column('SampleID', String, ForeignKey('Samples.ID'))
    microbe_id = Column('MicrobeID', Integer, ForeignKey('Microbes.ID'))

class Drugs(Base):
    __tablename__ = 'Drugs'
    id = Column('ID', String, primary_key=True)
    name = Column('Name', String)
    type = Column('Type', String)
    type_name = Column('TypeName', String)
    family = Column('Family', String)
    family_name = Column('FamilyName', String)
    class_code = Column('Class', String)
    class_name = Column('ClassName', String)

class Breakpoints(Base):
    __tablename__ = 'Breakpoints'
    drug_id = Column('DrugID', String, ForeignKey('Drugs.ID'), primary_key=True)
    genus = Column('Genus', String, ForeignKey('Microbes.Genus'), primary_key=True)
    mic_i = Column('MICI', Float)
    mic_r = Column('MICR', Float)
    drug = relationship("Drugs", backref='breakpoint')
    microbe = relationship('Microbes', backref='breakpoint')
    
class Tests(Base):
    __tablename__ = 'Tests'
    isolate_id = Column('IsolateID', String, ForeignKey('Isolates.ID'), primary_key=True)
    drug_id = Column('DrugID', String, ForeignKey('Drugs.ID'), primary_key=True)
    mic = Column('MIC', Float)
    corrected_mic = Column('CorrectedMIC', Float)
    sign = Column('sign', String)