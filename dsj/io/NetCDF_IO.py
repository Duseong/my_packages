'''
nc_io.py
this code is designed for writing NetCDF file quickly 
from xarray or custom options

MODIFICATION HISTORY:
    Duseong Jo, 23, DEC, 2020: VERSION 1.00
    Duseong Jo, 24, DEC, 2020: VERSION 1.10
      - Deal with a time dimension array
    Duseong Jo, 26, DEC, 2020: VERSION 2.00
      - Can use attributes from xarray while providing values to be saved
    Duseong Jo, 25, MAR, 2021: VERSION 3.00
    - Add nccopy function
'''

import numpy as np
from netCDF4 import Dataset
import xarray as xr
import datetime
import cftime
import os, subprocess

class Write_NC(object):
    
    '''
    NAME:
           Write_NC

    PURPOSE:
           Writing a NetCDF file

    INPUTS:
           vnames: a list which has variable names to be saved
           dnames: a list which has dimension names to be saved
           vdims: a dictionary which has dimension info of variables
           vatts: a dictionary which has attribute info of variables
           vatts_add: in case you want to add (or change) attribute info of variables 
                      in addtion to xarray info
           datts: a dictionary which has attribute info of dimensions
           datts_add: in case you want to add (or change) attribute info of dimensions 
                      in addtion to xarray info
           gatts: a dictionary which has global attribute info
           gatts_add: in case you want to add (or change) global attribute info
                      in addition to xarray info
           values: a dictionary which has values for NetCDF (can be numpy array or xarray)
           values_add: only in case "values" is from xarray: used as "values" instead, and
                "values" is used only for saving attributes or passing xarray values without manipulation           
           tunits: an unit for the time dimension (e.g., days since 2005-01-01 00:00:00)
           tcalendar: a calendar type for the time dimension (e.g., standard, noleap, etc.)
           dimindices: in case you want to extract a portion of xarray variable (e.g. [0,-1,:,:])
                       it can be a dictionary for different dimensions for different variables
           filename: NetCDF filename to be saved
           noadd_date: in case you don't want to add the creation date at the end of the filename
           fileformat: NetCDF fileformat - can be one of 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 
                                                         'NETCDF3_64BIT_OFFSET' or 'NETCDF3_64BIT_DATA'
           datatype: should be a string/dictionary (to manually assign) or None (to use xarray values)
                     examples are i4 (INT), i8 (INT64), f4 (FLOAT), f8 (DOUBLE), etc.
           compress: compress? - True or False
           complevel: compression level from 1 to 9, 9 uses the highest compresssion
           time_calc_only: in case you don't want to write a NC file, just need to calculate time array
           verbose: verbose output? - True or False
    '''
    
    def __init__(self, vnames=[], dnames=[], vdims=None, vatts={}, datts={}, gatts={}, values=None, 
                       vatts_add={}, datts_add={}, gatts_add={}, values_add={},
                       tunits='days since 2000-01-01 00:00:00', tcalendar='standard',
                       dimindices='[:]', filename='', noadd_date=False, fileformat='NETCDF4',
                       compress=False, complevel=9, datatype=None, time_calc_only=False, verbose=False ):
        
        self.vnames = vnames
        self.dnames = dnames
        self.vdims = vdims
        self.vatts = vatts
        self.vatts_add = vatts_add
        self.vatts_final = {}
        self.datts = datts
        self.datts_add = datts_add
        self.datts_final = {}
        self.gatts = gatts
        self.gatts_add = gatts_add
        self.gatts_final = {}
        self.values = values
        self.values_add = values_add
        self.tunits = tunits
        self.tcalendar = tcalendar
        self.dimindices = dimindices
        self.noadd_date = noadd_date
        self.fileformat = fileformat
        self.compress = compress
        self.complevel = complevel
        self.datatype = datatype
        self.verbose = verbose
        
        # ===== Construct output filename =====
        if filename == '':
            raise ValueError( '"filename" must be provided!' )
         
        CTime = datetime.datetime.now()
        CDate = str(CTime.year).zfill(4) + str(CTime.month).zfill(2) + \
                str(CTime.day).zfill(2)
        self.cdate = CDate
        if filename[-3:] == '.nc':
            if noadd_date:
                self.filename = filename
            else:
                self.filename = filename.replace( '.nc', '_c' + CDate + '.nc' )
        else:
            if noadd_date:
                self.filename = filename + '.nc'
            else:
                self.filename = filename + '_c' + CDate + '.nc'
        # ===== END Construct output filename =====
                
        # ===== Check if values are from xarray =====
        if type(self.values) == type(None):
            raise ValueError('The keyword "values" must be provided')
        elif type(self.values) in [ xr.core.dataset.Dataset, 
                                    xr.core.dataarray.DataArray ]:
            self.xarray_flag = True
        else:
            self.xarray_flag = False
        # ===== END Check if values are from xarray =====
        
        # ===== Check datatype =====
        if (self.datatype == None) & ~(self.xarray_flag):
            raise ValueError( '"datatype" must be provided unless xarray is provided')
        elif (self.datatype == None) & (self.xarray_flag):
            self.datatype = {}
        elif (type(self.datatype) == str) & ~(self.xarray_flag):
            tmpdatatype = self.datatype
            self.datatype = {}
            for key in list(self.values.keys()):
                self.datatype[key] = tmpdatatype
        # ===== END Check datatype =====
        
        # ===== check errors =====
        if self.vnames == []:
            raise ValueError( '"vname" must be provided' )
        for dname in self.dnames:
            if dname not in list( values.keys() ):
                raise ValueError( 'dimension variable ' + dname + \
                                  ' is not available in ' + 'values' )
        if (vatts_add != {}) & ~(self.xarray_flag):
            raise ValueError( 'vatts_add must be used with xarray. ' + \
                               'if you provided values manually, ' + \
                               'add attributes in vatts instead' )
        if (datts_add != {}) & ~(self.xarray_flag):
            raise ValueError( 'datts_add must be used with xarray. ' + \
                               'if you provided values manually, ' + \
                               'add attributes in datts instead' )
        if (gatts_add != {}) & ~(self.xarray_flag):
            raise ValueError( 'gatts_add must be used with xarray. ' + \
                               'if you provide values manually, ' + \
                               'add attributes in gatts instead' )
        if (values_add != {}) & ~(self.xarray_flag):
            raise ValueError( 'values_add must be used with xarray. ' + \
                               'if you provide values manually, ' + \
                               'add values in values keyword instead' )
        # ===== END check errors =====
        
        # ===== Construct time array if time is included in dimension variables =====
        if self.xarray_flag:
            if 'time' in self.values.coords.keys():
                if verbose:
                    print( "CALL construct_time_array function" )
                    print( "- criteria met: xarray with time dimension" )
                self.construct_time_array()
        elif 'time' in self.dnames:
            if verbose:
                print( "CALL construct_time_array function" )
                print( '- criteria met: time in "dnames"')
            self.construct_time_array()
        # ===== END Construct time array =====
        
        if not time_calc_only:
            self.write_nc()
        

    def write_nc(self):

        # ===== Open NetCDF file for writing 
        fid = Dataset( self.filename, 'w', format=self.fileformat )
        
        # ===== Create dimensions =====
        if self.xarray_flag:
            dname_list = self.values.coords.keys()
        else:
            dname_list = self.dnames
        
        if self.verbose:
            print( '===== Dimensions =====' )
        for dimname in list( dname_list ):
            self.datts_final[dimname] = {}
            if self.verbose:
                print( '- ' + dimname )
            
            # create dimension variables
            if dimname == 'time':
                # unlimited time dimension
                fid.createDimension( dimname )
            else:
                fid.createDimension( dimname, len(self.values[dimname]) )
            
            # assign datatype from xarray
            if self.xarray_flag:
                if dimname == 'time':
                    self.datatype[dimname] = 'f8'
                else:
                    self.datatype[dimname] = np.dtype( self.values[dimname].values[0] ).name
            
            
            # write dimension variables
            dimvar = fid.createVariable( dimname, self.datatype[dimname], (dimname,) )
            if dimname == 'time':
                dimvar[:] = self.time_array[:]
                print( 'time?', self.time_array[:] )
            else:
                dimvar[:] = np.copy( self.values[dimname] )
            # ===== add attributes =====
            if self.datts == {}:
                if self.xarray_flag:
                    for key in list( self.values[dimname].attrs.keys() ):
                        if dimname in list(self.datts_add.keys()):
                            if key in (self.datts_add[dimname].keys()):
                                dimvar.setncattr( key, self.datts_add[dimname][key] )
                                self.datts_final[dimname][key] = self.datts_add[dimname][key]
                            else:
                                dimvar.setncattr( key, self.values[dimname].attrs[key] )
                                self.datts_final[dimname][key] = self.values[dimname].attrs[key]
                        else:
                            dimvar.setncattr( key, self.values[dimname].attrs[key] )
                            self.datts_final[dimname][key] = self.values[dimname].attrs[key]
                        if self.verbose:
                            print( '-- ' + key + ', ' + self.datts_final[dimname][key] )
            else:
                if self.verbose & self.xarray_flag:
                    print( 'Dimension attributes of xarray ' + \
                           'will be overwritten by "datts" input' )
                for key in list(self.datts[dimname].keys()):
                    dimvar.setncattr( key, self.datts[dimname][key] )
                    self.datts_final[dimname][key] = self.datts[dimname][key]
                    if self.verbose:
                        print( '-- ' + key + ', ' + self.datts_final[dimname][key] )
            # Add addtional attributes
            if dimname in list(self.datts_add.keys()):
                for key in list(self.datts_add[dimname].keys()):
                    if key not in list(self.datts_final[dimname].keys()):
                        dimvar.setncattr( key, self.datts_add[dimname][key] )
                        self.datts_final[dimname][key] = self.datts_add[dimname][key]
            # ===== END add attributes =====

            if dimname == 'time':
                dimvar.setncattr( 'units', self.tunits )
                dimvar.setncattr( 'calendar', self.tcalendar )
                if self.verbose:
                    print( '-- ' + 'units' + ', ' + self.tunits )
                    print( '-- ' + 'calendar' + ', ' + self.tcalendar )
        # ===== END Create dimensions =====
        
        
        # ===== Create Variables =====
        var_dims = {}
        if self.verbose:
            print( '===== Variables =====' )
        for varname in self.vnames:
            self.vatts_final[varname] = {}
            # ===== Get dimension info for variable =====
            if self.verbose:
                print( '- ' + varname )
            if self.vdims == None:
                if self.xarray_flag:
                    if type(self.dimindices) == str:
                        tmp_dimindices = self.dimindices
                    elif type(self.dimindices) == dict:
                        tmp_dimindices = self.dimindices[varname]
                    else:
                        raise ValueError( '"dimindices" must be string or dictionary')
                    
                    if self.dimindices == '[:]':
                        var_dims[varname] = self.values[varname].dims
                    else:
                        var_dim_indices_array = \
                                tmp_dimindices.replace( ']' , '' ).replace( '[', '' ).split(',')
                        var_dims[varname] = "("

                        for ivd, vd in enumerate(var_dim_indices_array):
                            if ':' in vd:
                                var_dims[varname] += "'" + self.values[varname].dims[ivd] + "',"
                        var_dims[varname] = var_dims[varname][:-1] + ')'
                        var_dims[varname] = eval( var_dims[varname] )
                else:
                    raise ValueError( '"vdims" must be provided to specify ' + \
                                      'dimensions of variables to be saved' + \
                                      'unless xarray is passed')
            else:
                if type( self.vdims ) == tuple:
                    var_dims[varname] = self.vdims
                elif type( self.vdims ) == dict:
                    var_dims[varname] = self.vdims[varname]
                else:
                    raise ValueError( 'check "vdims"!!!' )
                tmp_dimindices = '[:]'
                
            if self.verbose:
                print( '-- dimensions:' + str(var_dims[varname]) )
            # ===== END Get dimension info for variable =====
            
            # ===== create Variable in NetCDF file =====
            # assign datatype from xarray
            if self.xarray_flag:
                self.datatype[varname] = np.dtype(self.values[varname].values.ravel()[0]).name
                
            var_temp = fid.createVariable( varname, self.datatype[varname], var_dims[varname],
                                           zlib=self.compress, complevel=self.complevel)
            if varname in list(self.values_add.keys()):
                str_tmp = 'var_temp[:] = self.values_add[varname]'
            else:
                str_tmp = 'var_temp[:] = self.values[varname]' + tmp_dimindices
            exec( str_tmp )
            
            # ===== add attributes =====
            if self.verbose:
                print ( '-- attributes: ')
            if self.vatts == {}:
                if self.xarray_flag:
                    for key in list( self.values[varname].attrs.keys() ):
                        if varname in list(self.vatts_add.keys()):
                            if key in list( self.vatts_add[varname].keys() ):
                                var_temp.setncattr( key, self.vatts_add[varname][key] )
                                self.vatts_final[varname][key] = self.vatts_add[varname][key]
                            else:
                                var_temp.setncattr( key, self.values[varname].attrs[key] )
                                self.vatts_final[varname][key] = self.values[varname].attrs[key]
                        else:
                            var_temp.setncattr( key, self.values[varname].attrs[key] )
                            self.vatts_final[varname][key] = self.values[varname].attrs[key]
                        if self.verbose:
                            print( '---- ' + key + ', ' + str(self.vatts_final[varname][key]) )
            else:
                if self.verbose & self.xarray_flag:
                    print( 'Dimension attributes of xarray ' + \
                           'will be overwritten by "vatts" input' )
                for key in self.vatts[varname]:
                    var_temp.setncattr( key, self.vatts[varname][key] )
                    self.vatts_final[dimname][key] = self.vatts[dimname][key]
                    if self.verbose:
                        print( '--- ' + key + ', ' + self.vatts[varname][key] )
            
            # Add additional attributes
            if varname in list(self.vatts_add.keys()):
                for key in list(self.vatts_add[varname].keys()):
                    if key not in list(self.vatts_final[varname].keys()):
                        var_temp.setncattr( key, self.vatts_add[varname][key] )
                        self.vatts_final[varname][key] = self.vatts_add[varname][key]
                    
            # ===== END add attributes =====
            
            # ===== END create Variable in NetCDF file =====
        self.var_dims = var_dims
        # ===== END Create Variables =====
        
        # ===== Save global attributes =====
        if self.xarray_flag:
            self.gatts = self.values.attrs
        
        if self.gatts != {}:
            self.gatts_final = {}
            if self.verbose:
                print( '===== Global attributes =====')
            for key in list( self.gatts.keys() ):
                if key in list( self.gatts_add.keys() ):
                    fid.setncattr( key, self.gatts_add[key] )
                    self.gatts_final[key] = self.gatts_add[key]
                else:
                    fid.setncattr( key, self.gatts[key] )
                    self.gatts_final[key] = self.gatts[key]
                if self.verbose:
                    print( '- ' + key + ', ' + self.gatts[key] )
        # add addtional global attributes
        for key in list(self.gatts_add.keys()):
            if key not in list(self.gatts_final.keys()):
                fid.setncattr( key, self.gatts_add[key] )
                self.gatts_final[key] = self.gatts_add[key]
        
        # ===== END Save global attributes =====
        
        
        fid.close()
        
    def construct_time_array(self):
        
        # ===== Retrieve time info =====
        if self.xarray_flag:
            if type( self.values['time'].values[0] ) == np.datetime64:
                self.time_year = self.values['time'].values.astype('datetime64[Y]').astype('int') + 1970
                self.time_month = self.values['time'].values.astype('datetime64[M]').astype('int') % 12 + 1
                self.time_day = ( self.values['time'].values.astype('datetime64[D]') - \
                                  self.values['time'].values.astype('datetime64[M]') + 1 ).astype('int')

            elif type( self.values['time'].values[0] ) in [cftime._cftime.DatetimeGregorian, 
                                                           cftime._cftime.DatetimeNoLeap]:
                self.time_year = np.zeros( len(self.values['time']) ).astype('int')
                self.time_month = np.zeros( len(self.values['time']) ).astype('int')
                self.time_day = np.zeros( len(self.values['time']) ).astype('int')

                for ti, timeraw in enumerate( self.values['time'].values ):
                    self.time_year[ti] = timeraw.year
                    self.time_month[ti] = timeraw.month
                    self.time_day[ti] = timeraw.day

            else:
                raise ValueError( 'Currently ' + str(type( self.values['time'].values[0] )) \
                                + ' is not supported! Check type of the time dimension' )
        
        else:
            if type( self.values['time'][0] ) == np.datetime64:
                self.time_year = self.values['time'].astype('datetime64[Y]').astype('int') + 1970
                self.time_month = self.values['time'].astype('datetime64[M]').astype('int') % 12 + 1
                self.time_day = ( self.values['time'].astype('datetime64[D]') - \
                                  self.values['time'].astype('datetime64[M]') + 1 ).astype('int')

            elif type( self.values['time'][0] ) in [cftime._cftime.DatetimeGregorian, 
                                                           cftime._cftime.DatetimeNoLeap]:
                self.time_year = np.zeros( len(self.values['time']) ).astype('int')
                self.time_month = np.zeros( len(self.values['time']) ).astype('int')
                self.time_day = np.zeros( len(self.values['time']) ).astype('int')

                for ti, timeraw in enumerate( self.values['time'] ):
                    self.time_year[ti] = timeraw.year
                    self.time_month[ti] = timeraw.month
                    self.time_day[ti] = timeraw.day

            else:
                raise ValueError( 'Currently ' + str(type( self.values['time'].values[0] )) \
                                + ' is not supported! Check type of the time dimension' )
        # ===== END Retrieve time info =====
        
        
        # ===== Calculate time array for NetCDF file save =====
        # Reference day
        Syear, Smonth, Sday = self.tunits.split(' ')[2].split('-')
        Sdate = datetime.date( int(Syear), int(Smonth), int(Sday) )
        
        # Number of days in this case
        if self.xarray_flag:
            self.time_array = np.zeros( len(self.values['time'].values) ).astype('f8')
        else:
            self.time_array = np.zeros( len(self.values['time']) ).astype('f8')
        
        if self.xarray_flag:
            for ti, timeraw in enumerate( self.values['time'].values ):
                self.time_array[ti] = ( datetime.date( self.time_year[ti], self.time_month[ti], 
                                                       self.time_day[ti] ) - Sdate ).days
        else:
            for ti, timeraw in enumerate( self.values['time'] ):
                self.time_array[ti] = ( datetime.date( self.time_year[ti], self.time_month[ti], 
                                                       self.time_day[ti] ) - Sdate ).days
        
        
        if self.verbose:
            print( '== Total length of time: ', len(self.time_array) )
            print( '---- from: ', self.time_year[0], self.time_month[0], 
                                  self.time_day[0], self.time_array[0] )
            print( '----   to: ', self.time_year[-1], self.time_month[-1], 
                                  self.time_day[-1], self.time_array[-1] )
            
        # ===== END Calculate time array for NetCDF file save =====
    
    
    
    # ===== Defining __call__ method =====
    def __call__(self):
        print( '=== filename ===')
        print( self.filename )
        print( '=== xarray_flag ===')
        print( self.xarray_flag )
        print( '=== vnames ===')
        print( self.vnames )
        print( '=== dnames ===')
        print( self.dnames )
        print( '=== vdims ===')
        print( self.vdims )
        print( '=== vatts ===')
        print( self.vatts )
        print( '=== vatts_add ===')
        print( self.vatts_add )
        print( '=== vatts_final ===')
        print( self.vatts_final )        
        print( '=== datts ===')
        print( self.datts )
        print( '=== datts_add ===')
        print( self.datts_add )
        print( '=== datts_final ===')
        print( self.datts_final )        
        print( '=== gatts ===')
        print( self.gatts )
        print( '=== gatts_add ===')
        print( self.gatts_add )
        print( '=== gatts_final ===')
        print( self.gatts_final )        
        print( '=== tunits ===')
        print( self.tunits )
        print( '=== tcalendar ===')
        print( self.tcalendar )
        print( '=== dimindices ===')
        print( self.dimindices )
        print( '=== noadd_date ===')
        print( self.noadd_date )
        print( '=== fileformat ===')
        print( self.fileformat )
        print( '=== compress ===')
        print( self.compress )
        print( '=== complevel ===')
        print( self.complevel )
        print( '=== datatype ===')
        print( self.datatype )
        print( '=== verbose ===')
        print( self.verbose )
        print( '=== values ===')
        for key in self.vnames + self.dnames:
            print( '---- name: ', key )
            if self.xarray_flag:
                print( '------   max:', np.max( self.values[key].values ) )
                print( '------   min:', np.min( self.values[key].values ) )
                print( '------ shape:', np.shape( self.values[key].values ) )
            else:
                print( '------   max:', np.max( self.values[key] ) )
                print( '------   min:', np.min( self.values[key] ) )
                print( '------ shape:', np.shape( self.values[key] ) )
               

def nccopy(filelist, ncformat='cdf5', use_ncks=True, debug_output=False):
    '''
    Available NetCDF file format:
    (1) 'classic'                  (NetCDF3 classic)
    (2) '64-bit offset'            (NetCDF3 64bit-offset)
    (3) 'cdf5'                     (NetCDF3 64bit-data)
    (4) 'netCDF-4'                 (NetCDF4 classic)
    (5) 'netCDF-4 classic model'   (NetCDF4)
    '''
    #subprocess.getoutput( 'echo "$USER"')
    #os.system( 'module load nco' )
    for file in filelist:
        print( 'processing: ', file, datetime.datetime.now() )
        #os.system( 'nccopy -k ' + ncformat + ' ' + file + ' ' + file + '2' )
        #os.system( 'mv ' + file + '2' + ' ' + file )
        if use_ncks:
            fid = open( 'tmp.sh', 'w' )
            fid.write( 'module load nco\n' )
            fid.write( 'ncks --fl_fmt=' + ncformat + ' ' + file + ' ' + file + '2' )
            fid.close()
            tmp1 = subprocess.getoutput( 'sh tmp.sh' )
        else:
            tmp1 = subprocess.getoutput( 'nccopy -k "' + ncformat + '" ' + file + ' ' + file + '2' )
        tmp2 = subprocess.getoutput( 'mv ' + file + '2' + ' ' + file )
        if debug_output:
            print( tmp1 )
            print( tmp2 )
            
            
            
            
    return  
