'''
CESM_aerosol.py
this code is designed for calculating aerosol concentrations from CESM (CAM-chem) output
Also calculates pressure, air density and air number density 

MODIFICATION HISTORY:
    Duseong Jo, 29, MAR, 2021: VERSION 1.00
    - Initial version
    Duseong Jo, 24, JUN, 2021: VERSION 1.10
    - Add custom options for regional history files
    Duseong Jo, 24, JUN, 2021: VERSION 1.20
    - Add a custom option to constrain levels
    Duseong Jo, 24, JUN, 2021: VERSION 1.30
    - Consider a case for 1-d array with only 'time' dimension
    Duseong Jo, 12, JUL, 2021: VERSION 1.31
    - Minor bug fix for air number density field
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
           T_name: Variable name for temperature
           PS_name: Variable name for surface pressure
           levels: specify which levels should be calculated (e.g., ':', '0', '2', '0:5', etc.)
           verbose: verbose output? - True or False
           
    '''
    def __init__(self, var, var_info, T_name='T', PS_name='PS', levels=':', verbose=True):
        
        self.T_name = T_name
        self.PS_name = PS_name
        self.levels = levels
        self.var = var
        
        if ( 'pressure' in list(var_info.keys()) ) & \
           ( 'airden' in list(var_info.keys()) ) & \
           ( 'airnum' in list(var_info.keys()) ):
            self.info_field = var_info
            self.Calculate_pressure(var_info, verbose, calc_index_only=True)
        else:
            self.Calculate_pressure(var_info, verbose)
        
        calc_str = "self.field = var" + self.str_t_index + " * self.info_field['airden'] * 1e9"
        exec(calc_str)
        
        
        
#class Calculate_pressure(object):
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
    #def __init__(self, var_info, verbose):
    def Calculate_pressure(self, var_info, verbose, calc_index_only=False ):
        

        if calc_index_only:
            dimension = self.var.dims
        else:
            dimension = var_info[self.T_name].dims

        str_pressure_index = '['
        str_hym_index = '['
        str_ps_index = '['
        str_t_index = '['
        for di, dim in enumerate(dimension):
            if dim == dimension[-1]:
                str_pressure_index += ':]'
                if len(dimension) == 1:
                    str_hym_index += '0]'
                else:
                    str_hym_index += ']'
                str_ps_index += ':]'
                str_t_index += ':]'
            elif dim == 'time':
                str_pressure_index += 'ti,'
                str_ps_index += 'ti,'
                str_hym_index += '0,'
                str_t_index += ':,'
            elif dim == 'lev':
                str_pressure_index += 'kii,'
                str_hym_index += 'ki'
                str_t_index += self.levels + ','
            else:
                str_pressure_index += ':,'
                str_ps_index += ':,'
                str_t_index += ':,'
        self.str_t_index = str_t_index    
            

        if not calc_index_only:

            if 'P0' in list(var_info.keys()):
                self.P0 = np.mean( var_info['P0'] )
            else:
                self.P0 = 100000     # Pa
                if verbose:
                    print( 'P0 is not available in info variable! using 100000 Pa as P0..' )

            if not self.T_name in list(var_info.keys()):
                raise ValueError( 'Temperature field (' + self.T_name + ') is not available!' )
            if not self.PS_name in list(var_info.keys()):
                raise ValueError( 'Surface pressure (' + self.PS_name + ') is not available!' )
            if not 'hyam' in list(var_info.keys()):
                raise ValueError( 'hybrid A coefficient (hyam) is not available!' )
            if not 'hybm' in list(var_info.keys()):
                raise ValueError( 'hybrid B coefficient (hybm) is not available!' )

            self.info_field = {}      

            # pressure in Pa
            calc_str = "self.info_field['pressure'] = np.copy( var_info[self.T_name].values" + \
                       str_t_index + " ) * 0.0"
            exec( calc_str )
            # air density in kg m-3
            calc_str = "self.info_field['airden'] = np.copy( var_info[self.T_name].values" + \
                       str_t_index + " ) * 0.0"
            exec( calc_str )
            # air number density in molecules m-3
            calc_str = "self.info_field['airnum'] = np.copy( var_info[self.T_name].values" + \
                       str_t_index + " ) * 0.0"
            exec( calc_str )


            # Calculate pressure            
            calc_str = ''
            if 'time' in dimension:
                calc_str += "for ti in np.arange( len(var_info['time']) ):\n"
            if 'lev' in dimension:
                calc_str += "    for kii, ki in enumerate( np.arange( len(var_info['lev']) )[" + self.levels + "]):\n"

            calc_str += "        self.info_field['pressure']" + str_pressure_index + ' = '
            calc_str += "var_info.hyam" + str_hym_index + "*self.P0 + var_info.hybm" + str_hym_index 
            calc_str += "*var_info." + self.PS_name + str_ps_index
            exec(calc_str)


            # calculate air density
            calc_str = ''
            calc_str += "self.info_field['airden'] = self.info_field['pressure'] / 287.05 / var_info[self.T_name]" + str_t_index
            exec(calc_str)

            # calculate air number density
            # n_air = ( A_v * p ) / ( R * T )
            calc_str = ''
            calc_str += "self.info_field['airnum'] = self.info_field['pressure'] * 6.022e23 / ( 8.314 * var_info[self.T_name]" + \
                        str_t_index + " )"
            exec(calc_str)


        
        
