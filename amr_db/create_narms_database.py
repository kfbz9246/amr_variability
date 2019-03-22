'''
Created on Jun 8, 2015

@author: kfbz9246
'''

from narms_parsing_tools import check_series_equality, replace_null
from pandas import read_excel, Series, isnull
from pandas.tslib import Timestamp
from re import compile, match
from sqlite3 import connect, register_adapter
from datetime import date



narms_data_file_path = 'amr_db/NARMS Data for Public Release 021215 FDA USDA 4-10-15.xlsx'

narms_database_ddl_file = open('setup_database.sql')
narms_ddl_commands = narms_database_ddl_file.read().split(';')
narms_database_ddl_file.close()

narms_database_connection = connect('narms_database.sqlite')
narms_database_cursor = narms_database_connection.cursor()
for narms_ddl_command in narms_ddl_commands:
    narms_database_cursor.execute(narms_ddl_command + ';')


narms_data = read_excel(narms_data_file_path, 'NARMS_Data_for_Public_Release')

print 'data read in'
narms_data_column_names = narms_data.columns.values
drug_re = compile('[A-Za-z]{3}$')
#sir_re = compile('([A-Za-z]{3}) SIR$')
sign_re = compile('([A-Za-z]{3}) Sign$')



samples = {}
next_microbe_id = 0
microbe_ids = {}
microbes = {}
isolates = {}
tests_table = []


counter = 0
for index, row in narms_data.iterrows():
    counter += 1
    if counter % 250 == 0:
        print counter
        
    if not isnull(row.loc['GROWTH']):
        if row.loc['GROWTH'].lower() == 'yes':
            mic_s = {}
            #sir_s = {}
            sign_s = {}
            sample_id = row.loc['SAMPLE_ID']
            sample_data = row.loc[['STATE', 'Month','Year', 'MEAT_TYPE', 'SELLBY_DATE', 'ORGANIC_IND', 'AGENCY', 'SOURCE', 'HOST_SPECIES']]
            if samples.has_key(sample_id):
                samples_sample_data = samples[sample_id]
                if not check_series_equality(samples_sample_data, sample_data) == True:
                    print 'uh oh'
                    print sample_id
                    print samples_sample_data
                    print sample_data
                    exit()
            else:
                samples[sample_id] = sample_data
                
            microbe = str(row.loc['GENUS_NAME']) + '#' + str(row.loc['SPECIES']) + '#' + str(row.loc['SEROTYPE']) + '#' + str(row.loc['ANTIGENIC_FORMULA'])
            microbe_data = row.loc[['GENUS', 'GENUS_NAME', 'SPECIES', 'SEROTYPE', 'ANTIGENIC_FORMULA']]
            if not microbe_ids.has_key(microbe):
                microbe_ids[microbe] = next_microbe_id
                next_microbe_id += 1
            microbe_id = microbe_ids[microbe]
            if microbes.has_key(microbe_id):
                microbes_microbe_data = microbes[microbe_id]
                if not check_series_equality(microbes_microbe_data, microbe_data) == True:
                    print 'uh oh 2'
                    print microbes_microbe_data
                    print microbe_data
                    exit()
            else:
                microbes[microbe_id] = microbe_data
            
            isolate_id = str(row.loc['SAMPLE_ID']) + '-' + (row.loc['GENUS'])
            isolate_data = row.loc[['PLATE', 'SAMPLE_ID']].append(Series({'MICROBE_ID': microbe_id}))
            if isolates.has_key(isolate_id):
                isolates_isolate_data = isolates[isolate_id]
                if not check_series_equality(isolates_isolate_data, isolate_data) == True:
                    print 'uh oh 3'
                    exit()
            else:
                isolates[isolate_id] = isolate_data
            
            for column_name in narms_data_column_names:
                if not isnull(row.loc[column_name]):
                    drug_re_match = match(drug_re, column_name)
                    if drug_re_match != None:   
                        mic_s[drug_re_match.group().lower()] = row.loc[column_name]
                    
                    sir_re_match = match(sir_re, column_name)
                    if sir_re_match != None:
                        sir_s[sir_re_match.group(1).lower()] = row.loc[column_name]
                        
                    sign_re_match = match(sign_re, column_name)
                    if sign_re_match != None:
                        sign_s[sign_re_match.group(1).lower()] = row.loc[column_name]
    

            for drug in mic_s.keys():
                if sir_s.has_key(drug):
                    if sign_s.has_key(drug):
                        tests_table.append((isolate_id, drug, mic_s.pop(drug), sir_s.pop(drug), sign_s.pop(drug)))
                    else:
                        tests_table.append((isolate_id, drug, mic_s.pop(drug), sir_s.pop(drug), None))
                else:
                    if sign_s.has_key(drug):
                        tests_table.append((isolate_id, drug, mic_s.pop(drug), None, sign_s.pop(drug)))
                    else:
                        tests_table.append((isolate_id, drug, mic_s.pop(drug), None, None))
            for drug in sir_s.keys():
                if sign_s.has_key(drug):
                    tests_table.append((isolate_id, drug, None, sir_s.pop(drug), sign_s.pop(drug)))
                else:
                    tests_table.append((isolate_id, drug, None, sir_s.pop(drug), None))
            
            for drug in sign_s.keys():
                tests_table.append((isolate_id, drug, None, None, sign_s.pop(drug)))
                    
            

del narms_data
sample_table = []
for sample_id in samples.keys():
    sample_data = samples.pop(sample_id)
    sample_table.append(replace_null([sample_id, sample_data.loc['SELLBY_DATE'], sample_data.loc['Month'], sample_data.loc['Year'], sample_data.loc['STATE'], sample_data.loc['MEAT_TYPE'], sample_data.loc['ORGANIC_IND'], sample_data.loc['AGENCY'], sample_data.loc['SOURCE'], sample_data.loc['HOST_SPECIES']]))
del samples

narms_database_cursor.executemany('INSERT INTO Samples (ID, SellByDate, Month, Year, State, MeatType, Organic, Agency, Source, HostSpecies) VALUES (?,?,?,?,?,?,?,?,?,?)', sample_table)
del sample_table

microbe_table = []
for microbe_id in microbes.keys():
    microbe_data = microbes.pop(microbe_id)
    microbe_table.append(replace_null([microbe_id, microbe_data.loc['GENUS'], microbe_data.loc['GENUS_NAME'], microbe_data.loc['SPECIES'], microbe_data.loc['SEROTYPE'], microbe_data.loc['ANTIGENIC_FORMULA']]))
del microbes
narms_database_cursor.executemany('INSERT INTO Microbes (ID, Genus, GenusName, Species, Serotype, AntigenicFormula) VALUES (?,?,?,?,?,?)', microbe_table)
del microbe_table

isolate_table = []
for isolate_id in isolates.keys():
    isolate_data = isolates.pop(isolate_id)
    isolate_table.append(replace_null([isolate_id, isolate_data.loc['PLATE'], isolate_data.loc['SAMPLE_ID'], isolate_data.loc['MICROBE_ID']]))
del(isolates)
narms_database_cursor.executemany('INSERT INTO Isolates (ID, Plate, SampleID, MicrobeID) VALUES (?,?,?,?)', isolate_table)
del isolate_table

narms_database_cursor.executemany('INSERT INTO Tests (IsolateID, Drug, MIC, SIR, Sign) VALUES (?,?,?,?,?)', tests_table)

narms_database_connection.commit()
narms_database_connection.close()
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        