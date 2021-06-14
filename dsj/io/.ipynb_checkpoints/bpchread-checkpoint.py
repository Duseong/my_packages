
'''
bpchread.py
this code is designed for GEOS-Chem bpch file

MODIFICATION HISTORY:
    dsj, 29, APR, 2015: VERSION 1.00
    dsj, 19, JUN, 2015: VERSION 1.10
                        - update for considering ts output
    dsj, 10, JUL, 2015: VERSION 1.20
                        - add seasonal mean function
    dsj, 15, JUL, 2015: VERSION 1.30
                        - add sum_tracers function
    dsj, 21, JUL, 2015: VERSION 1.40
                        - add calc_burden function
    dsj, 04, AUG, 2015: VERSION 1.41
                        - use concatenate for NAIR & BXHGHT
    dsj, 28, MAR, 2016: VERSION 1.50
                        - add diagdiff options for same tracer names 
                          of different categories
'''

from bpch import bpch
import numpy as np
import pdb
from dsj.time.tau import Tau
from dsj.io.string import list_string_add

class Bpch2(object):
    
    def __init__(self, filelist, diags=None, tracers=None,
                 tracerinfo=None, diaginfo=None, diagdiff=False,
                 sp_tracers=None ):

        self.files = filelist
        self.diags = diags
        self.read_tracers = tracers
        if diagdiff==True:
            if sp_tracers == None:
                self.tracers = list_string_add( diags, ['_']*len(diags), tracers )
            else:
                self.tracers = sp_tracers
        else:
            self.tracers = tracers
        self.diaginfo = diaginfo
        self.tracerinfo = tracerinfo
        self.data = {}
        self.unit = {}
        self.mw   = {}
        self.carbon = {}
        self.info = {}
        

    def read_bpch(self,elucidate=False,saveinfo=True,
                  ugm=False, ugmdiag='BXHGHT-$', ugmtracer='N(AIR)',
                  burden=False, burdendiag='BXHGHT-$',
                  burdentracer='BXHEIGHT', maxL=47):

        for f in np.arange(len(self.files)):
            if elucidate: print 'read file:', self.files[f]
            bcfile = bpch(self.files[f],
                          tracerinfo=self.tracerinfo,
                          diaginfo=self.diaginfo)

            for i in range(len(self.diags)):
                group = bcfile.groups[self.diags[i]]
                var = group.variables[self.read_tracers[i]]
                
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
                self.info['lons_bnds'] = group.variables['longitude_bounds']
                self.info['lats_bnds'] = group.variables['latitude_bounds']
                self.info['pres'] = group.variables['layer']
                self.info['area'] = group.variables['AREA']
                
                if f == 0:
                    self.info['tau0'] = np.array( group.variables['tau0'] )
                    self.info['tau1'] = np.array( group.variables['tau1'] )
                    self.info['taum'] = np.array( group.variables['time'] )
                    self.info['taum_bnds'] = \
                                 np.array( group.variables['time_bounds'] )
                else:
                    self.info['tau0'] = np.append( self.info['tau0'],
                                                   group.variables['tau0'] )
                    self.info['tau1'] = np.append( self.info['tau1'],
                                                   group.variables['tau1'] )
                    self.info['taum'] = np.append( self.info['taum'],
                                                   group.variables['time'] )
                    self.info['taum_bnds'] = np.append( self.info['taum_bnds'],
                                                       group.variables['time_bounds'] )

            if ugm:
                groupN = bcfile.groups[ugmdiag]
                varN = groupN.variables[ugmtracer]

                if f == 0:
                    Nair =  varN
                else:
                    Nair = np.concatenate( ( Nair, varN ), axis=0 )


            if burden:
                groupB = bcfile.groups[burdendiag]
                varB = groupB.variables[burdentracer]

                if f == 0:
                    BXHGHT =  varB
                else:
                    BXHGHT = np.concatenate( (BXHGHT, varB), axis=0 )

        if ugm:
            self.convert_mol_kg(Nair)

        # burden calculation should be after ugm
        if burden and not ugm:
            raise ValueError( 'burden calculation should turn on ugm' )
        if burden:
            self.calc_burden(BXHGHT,maxL=maxL)

            
    def convert_mol_kg(self,Nair):

        for i in np.arange(len(self.data)):
            self.data[self.tracers[i]] = \
                self.data[self.tracers[i]] * Nair \
                / 1e9 / 6.022e23 * self.mw[self.tracers[i]] * 1e9
            self.unit[self.tracers[i]] = 'ug/m3'


    def seasonal_mean(self):

        self.MAM = {}
        self.JJA = {}
        self.SON = {}
        self.DJF = {}
        MAMcount = {}
        JJAcount = {}
        SONcount = {}
        DJFcount = {}
        
        
        for i in range( len( self.data[self.tracers[0]][:,0,0,0] ) ):
            YMD = Tau( self.info['tau0'][i], 2 )

            
            for key in self.data.keys():

                if YMD.month == 3 or YMD.month == 4 or YMD.month == 5:

                    if key in self.MAM:
                        self.MAM[key] += self.data[key][i,:,:,:]
                        MAMcount[key] += 1
                        
                    else:
                        self.MAM[key] = self.data[key][i,:,:,:]
                        MAMcount[key] = 1

                elif YMD.month == 6 or YMD.month == 7 or YMD.month == 8:

                    if key in self.JJA:
                        self.JJA[key] += self.data[key][i,:,:,:]
                        JJAcount[key] += 1
                        
                    else:
                        self.JJA[key] = self.data[key][i,:,:,:]
                        JJAcount[key] = 1

                elif YMD.month == 9 or YMD.month == 10 or YMD.month == 11:

                    if key in self.SON:
                        self.SON[key] += self.data[key][i,:,:,:]
                        SONcount[key] += 1
                        
                    else:
                        self.SON[key] = self.data[key][i,:,:,:]
                        SONcount[key] = 1

                elif YMD.month == 12 or YMD.month == 1 or YMD.month == 2:

                    if key in self.DJF:
                        self.DJF[key] += self.data[key][i,:,:,:]
                        DJFcount[key] += 1
                        
                    else:
                        self.DJF[key] = self.data[key][i,:,:,:]
                        DJFcount[key] = 1

        for key in self.data.keys():

            self.MAM[key] = self.MAM[key] / MAMcount[key]
            self.JJA[key] = self.JJA[key] / JJAcount[key]
            self.SON[key] = self.SON[key] / SONcount[key]
            self.DJF[key] = self.DJF[key] / DJFcount[key]

    def sum_tracers(self):

        for i, key in enumerate( self.data.keys() ):

            if i == 0:
                self.datatotal = self.data[key]
            else:
                self.datatotal += self.data[key]

    def calc_burden(self,BXHGHT,maxL=47):

        self.burden = {}
        self.info['volume'] = np.zeros( np.shape( BXHGHT ) )


        for f in np.arange(len(self.files)):
            for l in np.arange( maxL ):
                self.info['volume'][f,l,:,:] = \
                      BXHGHT[f,l,:,:] * self.info['area']

        for key in self.data.keys():
            # 1E15 -> ug to Gg
            self.burden[key] = \
                    np.sum( self.data[key] * self.info['volume'] ) \
                    / 1e15 / len( self.info['tau0'] )



        


def read_tracerinfo( tinfofile ):

    TRCinfo = {}
    TRCinfo['Name'] = []
    TRCinfo['MW'] = []
    TRCinfo['CNumber'] = []
    TRCinfo['Number'] = []
    TRCinfo['Unit'] = []


    fid = open( tinfofile, 'r' )
    Outofloop = False
    while 1 is 1:
    
        strtmp = fid.readline()
    
        if strtmp[0] == '#':
            if Outofloop:
                break
            continue
        else:
            strtmps = strtmp.split()
            Outofloop = True

            TRCinfo['Name'].append( strtmps[0] )
            TRCinfo['MW'].append( strtmps[3] )
            TRCinfo['CNumber'].append( np.int( strtmps[4] ) )
            TRCinfo['Number'].append( np.int( strtmps[5] ) ) 
            TRCinfo['Unit'].append( strtmps[6] )

    return TRCinfo
