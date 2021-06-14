'''
chk.py
check variable

MODIFICATION HISTORY:
    dsj, 16, MAY, 2015: VERSION 1.00
    dsj, 28, MAY, 2015: VERSION 1.10
                      - use try-except and add additional information
    dsj, 10, JUL, 2015: VERSION 1.20
                      - Improve capability for array with nan values
    dsj, 15, MAY, 2019: VERSION 2.0
                      - Updated for python 3
    dsj, 20, AUG, 2019: VERSION 2.1
                      - add finite
'''
import numpy as np

def chk(X):

    try:
        print('type =', type(X))
    except:
        print('no type')
        
    try:
        print('len =', len(X))
    except:
        print('no len')
                
    try:
        print('shape =', np.shape(X))
    except:
        print('no shape')
        
    try:
        print('max =', np.nanmax(X))
    except:
        print('no max')
        
    try:
        print('min =', np.nanmin(X))
    except:
        print('no min')
        
    try:
        print('mean =', np.nanmean(X))
    except:
        print('no mean')

    try:
        print('Number of total size =', np.size(X))
    except:
        print('no Number of total size')
        
    try:
        print('Number of real =', np.sum(np.isreal(X)))
    except:
        print('no Number of real')

    try:
        print('Number of finite =', np.sum(np.isfinite(X)))
    except:
        print('no Number of finite')
        
    try:
        print('Number of Nan =', np.sum(np.isnan(X)))
    except:
        print('no Number of Nan')

# For python 2
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
        print 'max =', np.nanmax(X)
    except:
        print 'no max'
        
    try:
        print 'min =', np.nanmin(X)
    except:
        print 'no min'
        
    try:
        print 'mean =', np.nanmean(X)
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

'''
