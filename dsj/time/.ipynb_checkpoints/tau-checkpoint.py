'''
tau.py
handling tau values for GEOS-Chem

MODIFICATION HISTORY:
    dsj, 18, MAY, 2015: VERSION 1.00
'''

import swisseph as swe
import numpy as np


class Tau(object):
    '''
    NAME:
           Tau

    PURPOSE:
           Convert tau value to NYMD
           Convert NYMD value to tau
    '''

    def __init__(self,time,switch,firstyear=1985,firstmonth=01,firstday=01):
        self.inputtime = time
        self.firstyear = firstyear
        self.firstmonth = firstmonth
        self.firstday = firstday

        if switch == 'T' or switch == 1: # NYMD to Tau
            if isinstance(time,list):
                self.ymd = time[0]
                if len(time) == 2:
                    self.hms = time[1]
                else:
                    self.hms = 0
            else:
                self.ymd = time
                self.hms = 0
                
            self.nymd2tau()            


        if switch == 'N' or switch == 2: # Tau to NYMD
            self.tau = time

            self.tau2nymd()

            

    def nymd2tau(self):
        
        self.year = self.ymd / 10000
        self.month = ( self.ymd - self.year * 10000 ) / 100
        self.day = self.ymd - self.year * 10000 - self.month * 100 
        self.hour = self.hms / 10000
        self.minute = ( self.hms - self.hour * 10000 ) / 100
        self.second = self.hms - self.hour * 10000 - self.minute * 100

        tau0 = swe.julday( self.firstyear, self.firstmonth, self.firstday, \
                           0.0 ) * 24
        tau1 = swe.julday( self.year, self.month, self.day,
                           self.hour + self.minute / 60. + \
                             self.second / 3600. ) * 24
        self.tau = tau1 - tau0
        print self.tau
        

    def tau2nymd(self):
        tau0 = swe.julday( self.firstyear, self.firstmonth,
                                           self.firstday, 0.0 ) 
        tau1 = self.tau / 24. + tau0
        nymd = swe.revjul(tau1)
        self.year = nymd[0]
        self.month = nymd[1]
        self.day = nymd[2]
        nymd3 = np.around( nymd[3], 6 )
        self.hour = np.floor( nymd3 )
        self.minute = np.floor( ( nymd3 - self.hour ) * 3600 / 60 )
        self.second = ( nymd3 - self.hour - self.minute/60. ) * 3600.

        print self.year, self.month, self.day, \
              self.hour, self.minute, self.second
