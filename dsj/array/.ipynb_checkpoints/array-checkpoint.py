'''
array.py
This code is made for array manipulations

MODIFICATION HISTORY:
    dsj, 21, JUL, 2015: VERSION 1.00
    dsj, 12, JUL, 2019: VERSION 2.00
      - add print_dict

'''

import numpy as np

def dict_sum(dict,list):

    for i, key in enumerate( list ):

        if i == 0:
            total = dict[key]
        else:
            total = total + dict[key]
            
    return np.copy( total )


'''
dict = {}
for aa in ['A1','B1','C1']:
    dict[aa] = {}
    for bb in ['A2','B2','C2']:
        dict[aa][bb] = {}
        for cc in ['A3','B3','C3']:
            dict[aa][bb][cc] = {}
            for dd in ['A4','B4','C4']:
                dict[aa][bb][cc][dd] = {}
                for ee in ['A5','B5','C5']:
                    dict[aa][bb][cc][dd][ee] = aa + bb + cc + dd + ee
'''

def print_dict(dict,pos,names):
    
    keystring1 = 'list(dict'
    keystring2 = '.keys())'
    keystring3 = '[0]'
    keystring = {}
    keystring[0] = keystring1 + keystring2
    for i in np.arange(pos)+1:
        keystring[i] = keystring1
        for ii in np.arange(i):
            keystring[i] += '[' + keystring[ii] + keystring3 + ']'
        keystring[i] += keystring2

    forstring = ''
    printstring = 'print('
    printstring_ind = 'dict'
    availstring = ''
    for i in np.arange(pos):
        forstring += '    '*i + 'for key' + str(i) + ' in ' + keystring[i] + ':\n'
        printstring += 'key' + str(i) + ','
    for i in np.arange(pos):
        printstring_ind += '[key' + str(i) + ']'
        availstring += 'key' + str(i) + ','
    printstring += 'name,' + printstring_ind + '[name])'
    forstring += '    '*(pos+1) + 'for name in names:\n'
    forstring += '    '*(pos+2) + 'try:\n'
    forstring += '    '*(pos+3) + printstring + '\n'
    forstring += '    '*(pos+2) + 'except:\n'
    forstring += '    '*(pos+3) + 'print("not available:",' + availstring + 'name' + ')'
    
    exec(forstring)


