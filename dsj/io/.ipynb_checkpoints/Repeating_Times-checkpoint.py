'''
Repeating_Times.py
this code is designed for repeating times for CESM (CAM-chem) input
Single year file can be repeated for multi-years

MODIFICATION HISTORY:
    Duseong Jo, 16, MAR, 2021: VERSION 1.00
    - Initial version
    Duseong Jo, 18, MAR, 2021: VERSION 1.10
    - Change default NetCDF format for compatibility with CESM
'''

### Module import ###
import numpy as np
import xarray as xr
import datetime
from netCDF4 import Dataset
import subprocess


class Repeating_Times(object):
    '''
    NAME:
           Repeating_Times

    PURPOSE:
           Repeat single year (monthly) file for readability in CESM (CAM-chem)

    INPUTS:
           filename: NetCDF file name with a single year
           newfilename: New NetCDF file name with multiple years
           year_range: a list which has start and end years
           nc_file_format: NetCDF file format to be used in NetCDF4 library
    '''
    def __init__(self, filename, newfilename, Year_range=[2000,2020],
                 nc_file_format='NETCDF3_64BIT_DATA' ):

        ds_emis = xr.open_dataset( filename, decode_times=False )

        self.Year_all = np.arange( Year_range[0], Year_range[1]+1 )
        self.tunits = 'days since 1950-01-01 00:00:00'

        fid = Dataset( newfilename, 'w', format=nc_file_format )

        for di, dimname in enumerate( list(ds_emis.dims.keys()) ):
            if dimname == 'time':
                # unlimited time dimension
                fid.createDimension( dimname )
            else:
                fid.createDimension( dimname, len( ds_emis[dimname] ) )

            # Write dimension variables
            if dimname in list(ds_emis.coords):
                dimvar = fid.createVariable( dimname,
                                             np.dtype( ds_emis[dimname].values.flat[0] ).name,
                                             ds_emis[dimname].dims )
                if dimname == 'time':
                    date_array = []
                    time_array = []
                    Syear, Smonth, Sday = self.tunits.split(' ')[2].split('-')
                    Sdate = datetime.date( int(Syear), int(Smonth), int(Sday) )
                    day = 15
                    for year in self.Year_all:
                        for month in np.arange(12)+1:
                            date_array.append( year*10000 + month*100 + day )
                            time_array.append( ( datetime.date( year, month, day ) - Sdate ).days )
                    dimvar[:] = time_array
                    dimvar.units = self.tunits
                    dimvar.calendar = 'standard'

                    # Add Date
                    dimvar = fid.createVariable( 'date', 'int32', ('time',) )
                    dimvar[:] = date_array
                    dimvar.units = 'YYYYMMDD'
                    dimvar.long_name = 'Date'

                else:
                    dimvar[:] = ds_emis[dimname].values[:]
                    for key in list( ds_emis[dimname].attrs.keys() ):
                        dimvar.setncattr( key, ds_emis[dimname].attrs[key] )

        var_list = []
        var_shape = {}
        tile_shape = {}
        for dv in list(ds_emis.data_vars):
            if dv.lower() == 'date':
                continue

            if 'time' in ds_emis[dv].dims:
                var_list.append( dv )
                var_shape[dv] = []
                tile_shape[dv] = []
                for ii, dim in enumerate(ds_emis[dv].dims):
                    if dim == 'time':
                        var_shape[dv].append( ds_emis[dv].shape[ii] * len(self.Year_all) )
                        tile_shape[dv].append( len(self.Year_all) )
                    else:
                        var_shape[dv].append( ds_emis[dv].shape[ii] )
                        tile_shape[dv].append( 1 )

                var = fid.createVariable( dv, np.dtype( ds_emis[dv].values.flat[0] ).name,
                                          ds_emis[dv].dims )
                var[:] = np.tile( ds_emis[dv].values, tile_shape[dv] )


            else:
                var = fid.createVariable( dv, np.dtype( ds_emis[dv].values.flat[0] ).name,
                                          ds_emis[dv].dims )

            for key in list( ds_emis[dv].attrs.keys() ):
                var.setncattr( key, ds_emis[dv].attrs[key] )


        for key in list( ds_emis.attrs.keys() ):
            fid.setncattr( key, ds_emis.attrs[key] )

        fid.comment_repeat = '===== Below are created from repeating script ====='
        fid.repeated_from = filename
        fid.repeated_by = 'Repeating_Times.py'
        fid.repeated_time = str( datetime.datetime.now() )
        user_name = subprocess.getoutput( 'echo "$USER"')
        host_name = subprocess.getoutput( 'hostname -f' )
        fid.repeated_username = user_name + ' on ' + host_name

        fid.close()