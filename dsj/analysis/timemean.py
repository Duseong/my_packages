'''
timemean.py
This code is made for time mean (e.g. Daily mean, Diurnal variation)

MODIFICATION HISTORY:
    dsj, 03, AUG, 2015: VERSION 1.00
    dsj, 04, AUG, 2015: VERSION 1.01
                      - bug fix

'''

import numpy as np
import calendar
from dsj.time.tau import Tau
from dsj.array.chk import chk


class Dailymean(object):
    '''
    NAME:
           Dailymean

    PURPOSE:
           Calculate Daily mean

    INPUTS:
           Conc -> Concentration array as 1-d
           Taus -> Tau information for Conc
    '''

    def __init__(self,Conc,Taus):


        first = Tau( Taus[0], 'N' )
        end = Tau( Taus[-1], 'N' )

        firstTau = Tau( [ first.year * 10000 + first.month * 100 +
                        first.day ], 'T' )
        endTau = Tau( [ end.year * 10000 + end.month * 100 +
                        end.day ], 'T' )

        self.firstTau = firstTau.tau
        self.endTau = endTau.tau
        self.first = first.tau
        self.end = end.tau

        Ndays = ( self.endTau - self.firstTau ) / 24 + 1

        self.Dmean = np.zeros( (Ndays,) + np.shape(Conc)[1:] )
        Count = np.zeros( (Ndays,) + np.shape(Conc)[1:] )
        
        D = 0
        for i, t in enumerate(Taus):
            
            if i == 0:
                self.Dmean[D,:] += Conc[i,:]
                Count[D,:] += 1
                Day1 = Tau( t, 'N' ).day
            else:
                Day2 = Tau( t, 'N' ).day

                if Day1 == Day2:
                    self.Dmean[D,:] += Conc[i,:]
                    Count[D,:] += 1
                else:
                    D += 1
                    self.Dmean[D,:] += Conc[i,:]
                    Count[D,:] += 1

                Day1 = np.copy(Day2)
                

        self.Dmean = self.Dmean / Count

