'''
chk.py
check variable

MODIFICATION HISTORY:
    dsj, 16, MAY, 2015: VERSION 1.00
    dsj, 28, MAY, 2015: VERSION 1.10
                      - use try-except and add additional information
'''

import numpy as np

def chk(X):

    try:
        print 'type =', type(X)
    except:
        print 'no type'
        
    try:
        print 'len =', len(X)
    except:
        print 'no len'
                
    try:
        print 'shape =', np.shape(X)
    except:
        print 'no shape'
        
    try:
        print 'max =', np.max(X)
    except:
        print 'no max'
        
    try:
        print 'min =', np.min(X)
    except:
        print 'no min'
        
    try:
        print 'mean =', np.mean(X)
    except:
        print 'no mean'

    try:
        print 'Number of total size =', np.size(X)
    except:
        print 'no Number of total size'
        
    try:
        print 'Number of real =', np.sum(np.isreal(X))
    except:
        print 'no Number of real'
        
    try:
        print 'Number of Nan =', np.sum(np.isnan(X))
    except:
        print 'no Number of Nan'

