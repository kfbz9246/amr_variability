'''
Created on Dec 28, 2013

@author: kfbz9246
'''
class AbstractDataBaseTableConverter(object):
    def __init__(self,*positional_args, **keyword_args):
        self.row_generator = self.row_generator_func(*positional_args, **keyword_args)
        self.row_generator_is_not_empty = True
        self.out_buffer = bytearray()
        self.data_source_not_empty = True
        
    def readline(self, size = -1):
        return_val = ''
        if size < 0:
            if len(self.out_buffer) > 0:
                return_val = str(self.out_buffer)
            elif self.row_generator_is_not_empty:
                try:
                    return_val = self.row_generator.next()
                except StopIteration:
                    self.row_generator_is_not_empty = False
        else:
            return_val = self.read(size)
        return return_val
    
    def read(self, size = -1):
        return_val = ''
        if size < 0:
            return_val = str(self.out_buffer)
            for line in self.row_generator:
                return_val += line
            self.row_generator_is_not_empty = False
        else:
            while len(self.out_buffer) < size and self.row_generator_is_not_empty:
                try:
                    self.out_buffer.extend(bytearray(self.row_generator.next()))
                except StopIteration:
                    self.row_generator_is_not_empty = False
            return_val_array = bytearray()
            while len(return_val_array) < size and self.data_source_not_empty:
                try:
                    return_val_array.append(self.out_buffer.pop(0))
                except IndexError:
                    self.data_source_not_empty = False
            return_val = str(return_val_array)
        return return_val