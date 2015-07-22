'''
array.py
This code is made for array manipulations

MODIFICATION HISTORY:
    dsj, 21, JUL, 2015: VERSION 1.00

'''

import numpy as np

def dict_sum(dict,list):

    for i, key in enumerate( list ):

        if i == 0:
            total = dict[key]
        else:
            total = total + dict[key]
            
    return np.copy( total )



