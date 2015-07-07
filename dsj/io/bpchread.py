'''
bpchread.py
this code is designed for GEOS-Chem bpch file

MODIFICATION HISTORY:
    dsj, 29, APR, 2015: VERSION 1.00
    dsj, 19, JUN, 2015: VERSION 1.10
                        - update for considering ts output
'''

from bpch import bpch
import numpy as np
import pdb
from dsj.time.tau import Tau



class Bpch2(object):
    
    def __init__(self, filelist, diags=None, tracers=None,
                 tracerinfo=None, diaginfo=None):

        self.files = filelist
        self.diags = diags
        self.tracers = tracers
        self.diaginfo = diaginfo
        self.tracerinfo = tracerinfo
        self.data = {}
        self.unit = {}
        self.mw   = {}
        self.carbon = {}
        self.info = {}
        

    def read_bpch(self,elucidate=False,saveinfo=True, ugm=False,
                  ugmdiag='BXHGHT-$', ugmtracer='N(AIR)'):

        for f in np.arange(len(self.files)):
            if elucidate: print 'read file:', self.files[f]
            bcfile = bpch(self.files[f],
                          tracerinfo=self.tracerinfo,
                          diaginfo=self.diaginfo)

            for i in range(len(self.diags)):
                group = bcfile.groups[self.diags[i]]
                var = group.variables[self.tracers[i]]
                
                if f == 0:
                    self.data[self.tracers[i]] = var
                else:
                    self.data[self.tracers[i]] = np.concatenate(
                        ( self.data[self.tracers[i]], var), axis=0 )

                self.unit[self.tracers[i]] = var.units
                self.mw[self.tracers[i]] = var.kgpermole
                self.carbon[self.tracers[i]] = var.carbon

            if saveinfo:
                self.info['lons'] = group.variables['longitude']
                self.info['lats'] = group.variables['latitude']
                self.info['pres'] = group.variables['layer']
                self.info['area'] = group.variables['AREA']
                
                if f == 0:
                    self.info['tau0'] = np.array( group.variables['tau0'] )
                    self.info['tau1'] = np.array( group.variables['tau1'] )
                    self.info['taum'] = np.array( group.variables['time'] )
                else:
                    self.info['tau0'] = np.append( self.info['tau0'],
                                                   group.variables['tau0'] )
                    self.info['tau1'] = np.append( self.info['tau1'],
                                                   group.variables['tau1'] )
                    self.info['taum'] = np.append( self.info['taum'],
                                                   group.variables['time'] )

                

            if ugm:
                groupN = bcfile.groups[ugmdiag]
                varN = groupN.variables[ugmtracer]

                if f == 0:
                    Nair =  np.zeros( (len(self.files),) + \
                                      np.shape(np.squeeze(var))  )
                Nair[f,:] = np.squeeze(varN)

        if ugm:
            self.convert_mol_kg(Nair)

            
    def convert_mol_kg(self,Nair):

        for i in np.arange(len(self.data)):
            self.data[self.tracers[i]] = \
                self.data[self.tracers[i]] * Nair \
                / 1e9 / 6.022e23 * self.mw[self.tracers[i]] * 1e9
            self.unit[self.tracers[i]] = 'ug/m3'


    #def timeseries(self,elucidate=False,saveinfo=True,ugm=False,
    #               ugmdiag='BXHGHT-$', ugmtracer='N(AIR)'):
    #
    #    for f in np.arange(len(self.files)):
            
    
