'''
Created on Jun 8, 2015

@author: kfbz9246
'''

from narms_parsing_tools import check_series_equality, replace_null
from pandas import read_excel, Series, isnull, read_csv
from pandas.tslib import Timestamp
from re import compile, match
from sqlite3 import connect, register_adapter
from datetime import date
from math import log, isnan

mic_correction_map_file_path = '/home/code/amr_db/data/mic correction table.xlsx'
drug_data_file = '/home/code/amr_db/data/AbxTable.csv'
sir_file = '/home/code/amr_db/data/Breakpoints_Rec21Aug.csv'
narms_data_file_path_1 = '/home/code/amr_db/data/NARMSData_USDA-ARS_1997-2005-2.xlsx'
data_sheet_name_1 = 'Years 1997-2005'
narms_data_file_path_2 = '/home/code/amr_db/data/NARMSData_USDA-ARS_2006-2013.xlsx'
data_sheet_name_2 = 'Years 2006-2013'
narms_data_file_path_3 = '/home/code/amr_db/data/NARMSData_FDA-2.xlsx'
data_sheet_name_3 = 'NARMS_Data_for_Public_Release__'
narms_data_file_path_4 = '/home/code/amr_db/data/NARMSData_CDC v2.xlsx'
data_sheet_name_4 =  'NARMS_Data_for_Public_Release__'

narms_database_ddl_file = open('setup_database.sql')
narms_ddl_commands = narms_database_ddl_file.read().split(';')
narms_database_ddl_file.close()

narms_database_connection = connect('narms_database.sqlite')
narms_database_cursor = narms_database_connection.cursor()
for narms_ddl_command in narms_ddl_commands:
    narms_database_cursor.execute(narms_ddl_command + ';')
    

drug_data = read_csv(drug_data_file)
drug_table = []
for idx, row in drug_data.iterrows():
    drug_table.append([row.loc['Drug'].lower(), row.loc['Abx_name'], row.loc['Type'], row.loc['Type_Name'], row.loc['Family'], row.loc['Family_Name'], row.loc['Class'], row.loc['Class_Name']])
narms_database_cursor.executemany('INSERT INTO Drugs (ID, Name, Type, TypeName, Family, FamilyName, Class, ClassName) VALUES (?,?,?,?,?,?,?,?)', drug_table)




drug_re = compile('[A-Za-z]{3}$')
sign_re = compile('([A-Za-z]{3}) Sign$')



samples = {}
next_microbe_id = 0
microbe_ids = {}
microbes = {}
isolates = {}
tests_table = []




for data_file_path, data_sheet_name in [(narms_data_file_path_1, data_sheet_name_1), (narms_data_file_path_2, data_sheet_name_2), (narms_data_file_path_3, data_sheet_name_3), (narms_data_file_path_4, data_sheet_name_4)]:
    narms_data = read_excel(data_file_path, data_sheet_name)
    narms_data_column_names = narms_data.columns.values
    
    for index, row in narms_data.iterrows():         
                
        mic_s = {}
        sign_s = {}
        sample_id = row.loc['SAMPLE_ID']
        sample_data = row.loc[['STATE', 'Year', 'MEAT_TYPE', 'ORGANIC_IND', 'AGENCY', 'SOURCE', 'HOST_SPECIES']]
        if samples.has_key(sample_id):
            samples_sample_data = samples[sample_id]
            if not check_series_equality(samples_sample_data, sample_data) == True:
                sample_id=sample_id + '_2' 
                
        samples[sample_id] = sample_data
            
        microbe = str(row.loc['GENUS']) + '#' + str(row.loc['SPECIES']) + '#' + str(row.loc['SEROTYPE']) + '#' + str(row.loc['ANTIGENIC_FORMULA'])
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
        
        isolate_id = str(sample_id) + '-' + (row.loc['GENUS'])
        isolate_data = Series({'SAMPLE_ID': sample_id,'MICROBE_ID': microbe_id})
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
                    
                sign_re_match = match(sign_re, column_name)
                if sign_re_match != None:
                    sign_s[sign_re_match.group(1).lower()] = row.loc[column_name]
    

        for drug in mic_s.keys():
            if sign_s.has_key(drug):
                tests_table.append([isolate_id, drug, mic_s.pop(drug), sign_s.pop(drug)])
            else:
                tests_table.append([isolate_id, drug, mic_s.pop(drug), None])
                    
            

del narms_data
sample_table = []
for sample_id in samples.keys():
    sample_data = samples.pop(sample_id)
    sample_table.append(replace_null([sample_id, sample_data.loc['Year'], sample_data.loc['STATE'], sample_data.loc['MEAT_TYPE'], sample_data.loc['ORGANIC_IND'], sample_data.loc['AGENCY'], sample_data.loc['SOURCE'], sample_data.loc['HOST_SPECIES']]))
del samples

narms_database_cursor.executemany('INSERT INTO Samples (ID, Year, State, MeatType, Organic, Agency, Source, HostSpecies) VALUES (?,?,?,?,?,?,?,?)', sample_table)
del sample_table

microbe_table = []
for microbe_id in microbes.keys():
    microbe_data = microbes.pop(microbe_id)
    microbe_table.append(replace_null([microbe_id, microbe_data.loc['GENUS'], microbe_data.loc['SPECIES'], microbe_data.loc['SEROTYPE'], microbe_data.loc['ANTIGENIC_FORMULA']]))
del microbes
narms_database_cursor.executemany('INSERT INTO Microbes (ID, Genus, Species, Serotype, AntigenicFormula) VALUES (?,?,?,?,?)', microbe_table)
del microbe_table

isolate_table = []
for isolate_id in isolates.keys():
    isolate_data = isolates.pop(isolate_id)
    isolate_table.append(replace_null([isolate_id, isolate_data.loc['SAMPLE_ID'], isolate_data.loc['MICROBE_ID']]))
del(isolates)
narms_database_cursor.executemany('INSERT INTO Isolates (ID, SampleID, MicrobeID) VALUES (?,?,?)', isolate_table)
del isolate_table


mic_correction_map_data = read_excel(mic_correction_map_file_path)
mic_correction_map = {}
for row_name, row in mic_correction_map_data.iterrows():
    mic_correction_map[row.loc['Recorded Value']] = row.loc['Corrected value']

for row in tests_table:
    corrected_mic = log(mic_correction_map.get(row[2], row[2]),2)
    if row[3].find('>') > -1:
        corrected_mic += 1
    row.append(corrected_mic)

narms_database_cursor.executemany('INSERT INTO Tests (IsolateID, DrugID, MIC, Sign, CorrectedMIC) VALUES (?,?,?,?,?)', tests_table)

breakpoint_data = read_csv(sir_file)
breakpoint_table = []
for idx, row in breakpoint_data.iterrows():
    mic_i = row.loc['MIC_I']
    if mic_i ==0:
        mic_i = -100
    elif not isnan(mic_i):
        mic_i = log(row.loc['MIC_I'],2)
    breakpoint_table.append([row.loc['Drug'].lower(), row.loc['Genus'], mic_i, log(row.loc['MIC_R'],2)])
narms_database_cursor.executemany('INSERT INTO Breakpoints (DrugID, Genus, MICI, MICR) VALUES (?,?,?,?)', breakpoint_table)


narms_database_connection.commit()
narms_database_connection.close()
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        