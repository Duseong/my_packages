'''
budget.py
This code is made for scientific analysis especially for budget analysis

MODIFICATION HISTORY:
    dsj, 24, JUL, 2015: VERSION 1.00

'''

import numpy as np
import calendar
from dsj.time.tau import Tau



class Budget(object):
    '''
    NAME:
           Budget

    PURPOSE:
           Budget analysis

    INPUTS:
           Burdens -> Burdens calculated by Bpch class
           Wetdcv -> Convenctive wet deposition [kg/s]
           Wetdls -> Large-scale wet deposition [kg/s]
           Dryflx -> Dry deposition flux [molec/cm2/s]
           Info -> Informations given by Bpch class
           mw -> molecular weight [kg/mole]
           All above inputs should be dictionary
           Outputfile -> Outputfile for saving budget output
    '''
    
    def __init__(self,Burdens,Wetdcv,Wetdls,Dryflx,Info,mw,outputfile=None):

        Wetdep = {}
        Drydep = {}
        Avo = 6.022e23
        self.wetdep = {}
        self.drydep = {}
        self.burden = Burdens
        self.lifetime = {}


        # Construct days for a month & 3D Area array
        days = np.zeros( len(Info['tau0']), dtype='i' )
        days_4d = np.zeros( np.shape(Wetdcv[Wetdcv.keys()[0]] ), dtype='d' )
        days_3d = np.zeros( np.shape(Dryflx[Dryflx.keys()[0]]), dtype='d' )
        Area_3d = np.zeros( np.shape(Dryflx[Dryflx.keys()[0]]) )
        
        for i, t in enumerate(Info['tau0']):
            NYMD = Tau( t, 'N' )
            days[i] = calendar.monthrange(NYMD.year,NYMD.month)[1]
            days_4d[i,:,:,:] = days[i]
            days_3d[i,:,:] = days[i]
            # Area unit is [m2]
            Area_3d[i,:,:] = Info['area']

        # Calculate wet & dry deposition
        for key in Burdens.keys():
            Wetdep[key] = ( Wetdcv[key] + Wetdls[key] ) \
                        * 86400. * days_4d

            Drydep[key] = Dryflx[key+'df'] * 86400. * days_3d \
                        * Area_3d * 10000. / Avo * mw[key]

            # 1E9 -> kg to Gg
            self.wetdep[key] = np.sum( Wetdep[key] ) / 1E6
            self.drydep[key] = np.sum( Drydep[key] ) / 1E6
            
            # lifetime unit = [days-1]
            self.lifetime[key] = self.burden[key] * 365 / \
                               ( self.wetdep[key] + self.drydep[key] )

        self.wetdep['All'] = sum(self.wetdep.values())
        self.drydep['All'] = sum(self.drydep.values())
        self.lifetime['All'] = sum(self.burden.values()) * 365 / \
                               ( self.wetdep['All'] + self.drydep['All'] ) 


        if outputfile:
            f = open( outputfile, 'w' )
            f.write( 'species burden wetdep drydep lifetime\n' )
            
            for key in self.burden.keys():
                f.write( key + ' ' + str(self.burden[key]) + ' ' +
                         str(self.wetdep[key]) + ' ' +
                         str(self.drydep[key]) + ' ' +
                         str(self.lifetime[key]) + '\n' )

            f.close()
    



