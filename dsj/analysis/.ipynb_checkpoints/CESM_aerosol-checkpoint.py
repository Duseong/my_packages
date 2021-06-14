'''
CESM_aerosol.py
this code is designed for calculating aerosol concentrations from CESM (CAM-chem) output
Also calculates pressure, air density and air number density 

MODIFICATION HISTORY:
    Duseong Jo, 29, MAR, 2021: VERSION 1.00
    - Initial version
'''

### Module import ###
import numpy as np
import xarray as xr



class Convert_aerosol_ugm3(object):
    '''
    NAME:
           Convert_aerosol

    PURPOSE:
           Convert the unit of CESM aerosol output 
           from kg/kg to ug/m3

    INPUTS:
           var: xarray variable for aerosol field
           var_info: xarray variable for pressure and temperature fields
           verbose: verbose output? - True or False
           
    '''
    def __init__(self, var, var_info, verbose=True):
        
        if ( 'pressure' in list(var_info.keys()) ) & \
           ( 'airden' in list(var_info.keys()) ) & \
           ( 'airnum' in list(var_info.keys()) ):
            self.info_field = var_info
        else:
            self.info_field = Calculate_pressure(var_info, verbose).field
        
        self.field = var * self.info_field['airden'] * 1e9
        
        
        
class Calculate_pressure(object):
    '''
    NAME:
           Calculate_pressure

    PURPOSE:
           Calculates 3D pressure field using P0, PS, hyam, and hybm
           Air density [kg m-3] and Air number density are also calculated

    INPUTS:
           var_info: xarray variable for pressure and temperature fields
           verbose: verbose output? - True or False
    '''
    def __init__(self, var_info, verbose):
        
        if 'P0' in list(var_info.keys()):
            self.P0 = np.mean( var_info['P0'] )
        else:
            self.P0 = 100000     # Pa
            if verbose:
                print( 'P0 is not available in info variable! using 100000 Pa as P0..' )

        if not 'T' in list(var_info.keys()):
            raise ValueError( 'Temperature field (T) is not available!' )
        if not 'PS' in list(var_info.keys()):
            raise ValueError( 'Surface pressure (PS) is not available!' )
        if not 'hyam' in list(var_info.keys()):
            raise ValueError( 'hybrid A coefficient (hyam) is not available!' )
        if not 'hybm' in list(var_info.keys()):
            raise ValueError( 'hybrid B coefficient (hybm) is not available!' )

        self.field = {}
        # pressure in Pa
        self.field['pressure'] = np.copy( var_info['T'].values ) * 0.0
        # air density in kg m-3
        self.field['airden'] = np.copy( var_info['T'].values ) * 0.0
        # air number density in molecules m-3
        self.field['airnum'] = np.copy( var_info['T'].values ) * 0.0
        
        dimension = var_info['T'].dims

        str_pressure_index = '['
        str_hym_index = '['
        str_ps_index = '['
        for di, dim in enumerate(dimension):
            if dim == 'time':
                str_pressure_index += 'ti,'
                str_ps_index += 'ti,'
                str_hym_index += '0,'
            elif dim == 'lev':
                str_pressure_index += 'ki,'
                str_hym_index += 'ki'
            elif dim == dimension[-1]:
                str_pressure_index += ':]'
                str_hym_index += ']'
                str_ps_index += ':]'
            else:
                str_pressure_index += ':,'
                str_ps_index += ':,'
                        
        calc_str = ''
        if 'time' in dimension:
            calc_str += "for ti in np.arange( len(var_info['time']) ):\n"
        if 'lev' in dimension:
            calc_str += "    for ki in np.arange( len(var_info['lev']) ):\n"
        calc_str += "        self.field['pressure']" + str_pressure_index + ' = '
        calc_str += "var_info.hyam" + str_hym_index + "*self.P0 + var_info.hybm" + str_hym_index 
        calc_str += "*var_info.PS" + str_ps_index
        exec(calc_str)

        self.field['airden'] = self.field['pressure'] / 287.05 / var_info['T']
        # n_air = ( A_v * p ) / ( R * T )
        self.field['airnum'] = self.field['pressure'] * 6.022e23 / ( 8.314 * var_info['T'] )
        