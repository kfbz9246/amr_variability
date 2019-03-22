'''
Created on Jun 9, 2015

@author: kfbz9246
'''
from abstract_database_table_converter import AbstractDataBaseTableConverter
from pandas import isnull
from pandas.tslib import Timestamp
from datetime import date

def check_series_equality(series_1, series_2):
    return_value= True
    for item_1, item_2 in zip(series_1, series_2):
        if not (isnull(item_1) and isnull(item_2)):
            if not item_1 == item_2:
                return_value = False
    return return_value

def replace_null(a_list):
    for index in range(len(a_list)):
        a_index_item = a_list[index]
        if isnull(a_index_item):
            a_list[index] = None
            
        elif type(a_index_item) == Timestamp:
            a_list[index] = date(a_index_item.year, a_index_item.month, a_index_item.day)
    return tuple(a_list)

class SampleTableConverter(AbstractDataBaseTableConverter):
    @staticmethod
    def row_generator_func(sample_info):
        for sample_id in sample_info.keys():
            sample_data = sample_info[sample_id]
            yield str(sample_id) + '\t' + str(sample_data.loc['RECEIVED_DATE']) + '\t' + str(sample_data.loc['SELLBY_DATE']) + '\t' + str(sample_data.loc['Month']) + '\t' + str(sample_data.loc['Year']) + '\t' + str(sample_data.loc['STATE']) + '\t' + str(sample_data.loc['MEAT_TYPE']) + '\t' + str(sample_data.loc['ORGANIC_IND']) + '\t' + str(sample_data.loc['AGENCY']) + '\t' + str(sample_data.loc['SOURCE']) + '\t' + str(sample_data.loc['HOST_SPECIES']) + '\n'


class MicrobeTableConverter(AbstractDataBaseTableConverter):
    @staticmethod
    def row_generator_func(microbe_info):
        for microbe_id in microbe_info.keys():
            microbe_data = microbe_info[microbe_id]
            yield str(microbe_id) + '\t' + str(microbe_data.loc['GENUS']) + '\t' + str(microbe_data.loc['GENUS_NAME']) + '\t' + str(microbe_data.loc['SPECIES']) + '\t' + str(microbe_data.loc['SEROTYPE'])+ '\t' + str(microbe_data.loc['ANTIGENIC_FORMULA']) + '\n'
            

class IsolateTableConverter(AbstractDataBaseTableConverter):
    @staticmethod
    def row_generator_func(isolate_info):
        for isolate_id in isolate_info.keys():
            isolate_data = isolate_info[isolate_id]
            yield str(isolate_id) + '\t' + str(isolate_data.loc['PLATE']) + '\t' + str(isolate_data.loc['SAMPLE_ID']) + '\t' + str(isolate_data.loc['MICROBE_ID']) + '\n'

class TestTableConverter(AbstractDataBaseTableConverter):
    @staticmethod
    def row_generator_func(test_table):
        for test in test_table:
            yield str(test[0])  + '\t' + str(test[1])  + '\t' + str(test[2]) + '\t' + str(test[3]) + '\n'
