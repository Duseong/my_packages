'''
Species_Mapping.py
this code is designed for recalculating emissions for CAM-chem
because CAM-chem and emission inventory have different chemical species
(1) Species mapping (class sp_map)

MODIFICATION HISTORY:
    Duseong Jo, 5, MAR, 2021: VERSION 1.00
    - Initial version
    Duseong Jo, 15, MAR, 2021: VERSION 2.00
    - Update for aircraft emisssions
    Duseong Jo, 16, MAR, 2021: VERSION 2.10
    - Bug fix for additonal variables
    Duseong Jo, 17, MAR, 2021: VERSION 2.20
    - Update for adding "date" variable
    Duseong Jo, 18, MAR, 2021: VERSION 2.30
    - Minor bug fix in printing lines for user_nl_cam
    Duseong Jo, 18, MAR, 2021: VERSION 2.40
    - Change default NetCDF format for compatibility with CESM
    Duseong Jo, 19, MAR, 2021: VERSION 2.50
    - Change default vertical dimension name for compatibility with CESM
    Duseong Jo, 23, MAR, 2021: VERSION 2.60
    - Add an additional error diagnostic
    Duseong Jo, 08, APR, 2021: VERSION 2.70
    - Add an additional keyword for reading file only
    Duseong Jo, 16, APR, 2021: VERSION 2.80
    - Add an additional attribute for the mapping equation
    Duseong Jo, 22, APR, 2021: VERSION 2.90
    - Add an additional keyword for custom date option
    Duseong Jo, 19, JUL, 2021: VERSION 3.00
    - Consider a case that some species has missing source categories when lumping
    Duseong Jo, 24, SEP, 2021: VERSION 3.10
    - Add an additional unit to deal with HTAPv3 emissions (ton/month)
    Duseong Jo, 24, SEP, 2021: VERSION 4.00
    - Add an additional keyword to add date_array variable manually
    Duseong Jo, 22, NOV, 2021: VERSION 4.10
    - Add if statements to check different grid area variables
'''

### Module import ###
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import glob
import cftime
import calendar
import datetime
import subprocess
from dsj.io.get_sp import get_sp


class sp_map(object):
    '''
    NAME:
           sp_map

    PURPOSE:
           Species mapping + unit conversion + vertical redistribution for CESM/CAM-chem

    INPUTS:
           mapping_file: input file for species mapping
           nc_file_format: NetCDF file format to be used in NetCDF4 library
           print_user_nl_cam: if True, print lines for the use of user_nl_cam
           keep_sector: if True, the script will keep the sectors
                        if False, the script will sum up all the sectors and make it one total
           read_file_only: if True, the script will read mapping file without calculation
           date_array: in case of adding date array manually
           day_of_each_month: if provided (integer), use it for day of each month, 
                               overwrite the day of month in the original file
           ignore_warning: if True, ignore warning messages
           verbose: Display detailed information on what is being done
    '''
    
    def __init__(self, mapping_file, nc_file_format='NETCDF3_64BIT_DATA', print_user_nl_cam=True, 
                 keep_sector=True, read_file_only=False, date_array=None, day_of_each_month=None, 
                 ignore_warning=False, verbose=False):
        
        self.mapping_file = mapping_file
        self.nc_file_format = nc_file_format
        self.print_user_nl_cam = print_user_nl_cam
        self.keep_sector = keep_sector
        self.date_array = date_array
        self.day_of_each_month = day_of_each_month
        self.ignore_warning = ignore_warning
        self.verbose = verbose
        self.conversion_factor = {}
        self.mw_Inventory = {}
        self.mass_per_particle = {}
        self.altitude_interfaces = {}
        self.altitude_intervals = {}
        
        if self.print_user_nl_cam:
            self.ext_frc_specifier = "ext_frc_specifier = "
            self.srf_emis_specifier = "srf_emis_specifier = "          
            
        # Read species mapping file
        self.read_mapping_file()
        
        # Remapping
        if not read_file_only:
            self.remapping()
        
        
    def read_mapping_file(self):
        # ========================================================================
        # ====================== Read species mapping file =======================
        # ========================================================================

        fid = open( self.mapping_file, 'r' )
        self.file_all = fid.readlines()
        fid.close()

        Reading_Part1 = True; Reading_Part2 = False; Reading_Part3 = False
        self.species_CESM = []
        self.mw_CESM = {}; self.mw_Inventory_force = {}; self.yield_type = {}
        self.N_equations = {}; self.eq_frac = {}; self.eq_species = {}
        self.eq_operators = {}; self.sectors = {}; self.dimensions = {}
        self.altitudes = {}; self.alt_fractions = {}; self.eq_attribute = {}
        self.diameter = {}; self.density = {}
        yield_types_list = ['mass','molar']
        
        Start_Reading_Part2 = True
        
        for ii, line in enumerate(self.file_all):
            if 'file locations end' in line.lower():
                Reading_Part1 = False
                Reading_Part2 = True
                continue
            elif 'mappaing equations end' in line.lower():
                Reading_Part2 = False
                Reading_Part3 = True    

            if line[0] == '#':
                continue

            # Reading Part 1
            if Reading_Part1:
                exec( 'self.' + line )
                
            # Reading Part 2    
            if Reading_Part2:
                if Start_Reading_Part2:
                    Start_Reading_Part2 = False
                    self.Source_filelist = glob.glob( self.Source_fileformat )
                    self.Source_sp, self.Source_files = get_sp( self.Source_filelist )

                    if self.Source_sp == []:
                        raise ValueError( 'Check Source file format!')        
                
                tmp = line.split()
                sp_tmp = tmp[0].split('[')[0]
                self.species_CESM.append( sp_tmp )
                self.eq_attribute[sp_tmp] = line.replace('\n','')
                self.mw_CESM[sp_tmp] = float( tmp[0].split('[')[1][:-1] )

                yield_type_check = True
                for yi in yield_types_list:
                    if yi in tmp:
                        yield_type_check = False
                        self.yield_type[sp_tmp] = yi
                        yield_index = np.where( np.char.equal( tmp, yi ) )[0][0]
                        self.N_equations[sp_tmp] = ( yield_index - 1) / 2

                        self.eq_species[sp_tmp] = []; self.eq_frac[sp_tmp] = []; self.eq_operators[sp_tmp] = []
                        for ei in np.arange(self.N_equations[sp_tmp]):
                            if tmp[2+np.int(ei*2)].split('*')[-1] in self.Source_sp:
                                self.eq_species[sp_tmp].append( tmp[2+np.int(ei*2)].split('*')[-1] )
                                pos_sp_border = tmp[2+np.int(ei*2)].rfind('*')
                                exec( 'self.eq_frac[sp_tmp].append( ' + tmp[2+np.int(ei*2)][:pos_sp_border] + ')' )
                                if ei > 0:
                                    self.eq_operators[sp_tmp].append( tmp[1+np.int(ei*2)] )

                            else:
                                raise ValueError( 'Check Source file format and species name!' + \
                                                  'Error in reading: ', tmp[2+np.int(ei*2)].split('*')[-1] )
                        
                        
                        self.sectors[sp_tmp] = tmp[yield_index+1].replace('[','').replace(']','').split(',')
                        self.dimensions[sp_tmp] = tmp[yield_index+2].replace('{','').replace('}','')

                        if self.dimensions[sp_tmp] == 'vert':
                            self.altitudes[sp_tmp] = np.array( tmp[yield_index+3].replace('[','').replace(']','').split(',') ).astype('f')
                            self.alt_fractions[sp_tmp] = np.array( tmp[yield_index+4].replace('[','').replace(']','').split(',') ).astype('f')
                            add_n = 5
                        else:
                            add_n = 3

                        if sp_tmp[:4] == 'num_':
                            self.diameter[sp_tmp] = float( tmp[yield_index+add_n] )
                            self.density[sp_tmp] = float( tmp[yield_index+add_n+1] )

                if yield_type_check:
                    raise ValueError( 'No information for yield type! ' + \
                                      'It should be "mass" or "molar"\n' + \
                                      'Check your input file' )

            # Reading Part 3
            if Reading_Part3:
                tmp = line.split()
                self.mw_Inventory_force[tmp[0]] = float( tmp[2] )
    

    def remapping(self):

        self.ds_source = {}
        self.ds_source_time = {}
        self.Destination_files = {}
        self.nc_dims = {}
        if self.species_list.lower() == 'all':
            self.species_to_process = self.species_CESM
        else:
            self.species_to_process = self.species_list

        First_read = True
        for sp_cesm in self.species_to_process:
            if self.verbose:
                print( 'Processing: ', sp_cesm, datetime.datetime.now() )

            # Read source files
            for sp_src in self.eq_species[sp_cesm]:
                if sp_src not in self.ds_source.keys():
                    self.ds_source[sp_src] = xr.open_dataset( self.Source_files[sp_src],
                                                              decode_times=False)
                    self.ds_source_time[sp_src] = xr.open_dataset( self.Source_files[sp_src] )

                    # Check whether the grid is FV or SE
                    if First_read:
                        self.dim_var = {}
                        if 'ncol' in list( self.ds_source[sp_src].dims ): 
                            self.grid_type = 'SE'
                            if 'grid_area_rad2' in list( self.ds_source[sp_src].variables ):
                                self.dim_var['grid_area_rad2'] = self.ds_source[sp_src]['grid_area_rad2'].values
                            elif 'area' in list( self.ds_source[sp_src].variables ):
                                self.dim_var['grid_area_rad2'] = self.ds_source[sp_src]['area'].values
                                if 'units' in self.ds_source[sp_src]['area'].attrs:
                                    if self.ds_source[sp_src]['area'].units != 'radians^2':
                                        if not self.ignore_warning:
                                            print( 'Warning: Check units in area, it should be "radians^2"' )
                                else:
                                    if not self.ignore_warning:
                                        print( 'Warning: Unit attribute is not availabe for area, it is assumed to be radians^2' )
                            else:
                                raise ValueError( 'Grid area variable is not available!' )
                                
                        elif 'lat' in list( self.ds_source[sp_src].dims ):
                            self.grid_type = 'FV'
                            self.dim_var['lat'] = self.ds_source[sp_src]['lat'].values
                            self.dim_var['lon'] = self.ds_source[sp_src]['lon'].values
                        else:
                            raise ValueError( 'Check source dimensions!' )                    
                        First_read = False
                        

                        
                    
                if sp_src == self.eq_species[sp_cesm][0]:
                    self.nc_dims[sp_cesm] = list( self.ds_source[sp_src].dims )

                # get mw info if available or use values from the input file
                if sp_src in self.mw_Inventory_force.keys():
                    self.mw_Inventory[sp_src] = self.mw_Inventory_force[sp_src]
                else:
                    nc_mw_check = False
                    for dvar in self.ds_source[sp_src].data_vars:
                        for key in list( self.ds_source[sp_src][dvar].attrs.keys() ):
                            if key in ['molecular_weight', 'molecular_weights']:
                                self.mw_Inventory[sp_src] = self.ds_source[sp_src][dvar].attrs[key]
                                nc_mw_check = True
                    for key in list( self.ds_source[sp_src].attrs.keys() ):
                        if key in ['molecular_weight', 'molecular_weights']:
                            self.mw_Inventory[sp_src] = self.ds_source[sp_src].attrs[key]
                            nc_mw_check = True

                    if not nc_mw_check:
                        raise ValueError('You must provide molecular weights of species in source files' + \
                                         '\nif those are not available in NetCDF file attributes')

                    if type( self.mw_Inventory[sp_src] ) == str:
                        try:
                            self.mw_Inventory[sp_src] = np.copy( self.mw_Inventory[sp_src].replace('.f'),'' ).astype('f')
                        except:
                            self.mw_Inventory[sp_src] = np.copy( self.mw_Inventory[sp_src] ).astype('f')
                    elif type( self.mw_Inventory[sp_src] ) in [int, float, np.float32, np.float64]:
                        self.mw_Inventory[sp_src] = np.copy( self.mw_Inventory[sp_src] ).astype('f')
                    else:
                        raise ValueError( 'Check molecular weight variable in NetCDF file!')

                    if self.mw_Inventory[sp_src] < 1:
                        if not self.ignore_warning:
                            print( 'It looks like molecular weights in source NetCDF file ' + \
                                   'have the unit of kg/mole, converting it to g/mole')
                        self.mw_Inventory[sp_src] *= 1000.


            # ===== Open NetCDF file to write =====
            self.Destination_files[sp_cesm] = self.Destination_fileformat.replace('SPC',sp_cesm)
            fid = Dataset( self.Destination_files[sp_cesm], 'w', format=self.nc_file_format )

            # construct text for the use of user_nl_cam
            if self.print_user_nl_cam:
                bin_name = ''
                for abin in ['a1','a2','a3','a4']:
                    if '_' + abin in sp_cesm:
                        bin_name = abin

                sp_name = ''
                for nn in sp_cesm.split('_'):
                    if nn == 'num':
                        sp_name = nn
                        break
                if sp_name == '':
                    sp_name = sp_cesm.split('_')[0]

                if bin_name == '':
                    sp_bin = sp_name
                else:
                    sp_bin = sp_name + '_' + bin_name

                if self.dimensions[sp_cesm] in ['vert', 'aircraft']:
                    self.ext_frc_specifier += "'" + sp_bin + " -> " + \
                                              self.Destination_files[sp_cesm] + "',\n" + \
                                              "                    "
                elif self.dimensions[sp_cesm] == 'sfc':
                    self.srf_emis_specifier += "'" + sp_bin + " -> " + \
                                              self.Destination_files[sp_cesm] + "',\n" + \
                                              "                     "
                else:
                    raise ValueError( 'Check input file for dimensions! - sfc/vert/aircraft' )

            # ===== Create Dimensions =====
            for di, dimname in enumerate( self.nc_dims[sp_cesm] ):
                # level dimension will be replaced by "altitude" for CESM compatibility
                if dimname in ['level']:
                    continue
                    
                if dimname == 'time':
                    # unlimited time dimension
                    fid.createDimension( dimname )
                else:
                    fid.createDimension( dimname, len( self.ds_source[sp_src][dimname] ) )

                # Write dimension variables
                dimvar = fid.createVariable( dimname, 
                                             np.dtype( self.ds_source[sp_src][dimname].values.flat[0] ).name,
                                             self.ds_source[sp_src][dimname].dims )
                dimvar[:] = self.ds_source[sp_src][dimname].values
                for key in list( self.ds_source[sp_src][dimname].attrs.keys() ):
                    dimvar.setncattr( key, self.ds_source[sp_src][dimname].attrs[key] )


            # Add vertical dimension if applied
            if self.dimensions[sp_cesm] in ['vert','aircraft']:
                if self.dimensions[sp_cesm] == 'aircraft':
                    if 'level' in list( self.ds_source[sp_src].coords ):
                        vert_var_name = 'level'
                        self.altitudes[sp_cesm] = self.ds_source[sp_src][vert_var_name].values
                    else:
                        raise ValueError( 'Check vertical dimension of aircraft emissions!\n' + \
                                          'Currently only "level" is supported' )
                
                fid.createDimension( 'altitude', len( self.altitudes[sp_cesm] ) )
                fid.createDimension( 'altitude_int', len( self.altitudes[sp_cesm] ) + 1 )
                    
                # Write dimension variables for vertical dimension
                altvar = fid.createVariable( 'altitude', 'f8', ('altitude',) )
                altivar = fid.createVariable( 'altitude_int', 'f8', ('altitude_int',) )
                altvar[:] = self.altitudes[sp_cesm]

                # calculate altitude interfaces and intervals
                self.altitude_interfaces[sp_cesm] = np.zeros( len(self.altitudes[sp_cesm])+1 )
                self.altitude_intervals[sp_cesm] = np.zeros( len(self.altitudes[sp_cesm]) )
                for vi in np.arange( len(self.altitudes[sp_cesm])-1 ):
                    self.altitude_interfaces[sp_cesm][vi] = self.altitudes[sp_cesm][vi] - \
                                        (self.altitudes[sp_cesm][vi+1] - self.altitudes[sp_cesm][vi]) / 2.
                self.altitude_interfaces[sp_cesm][-1] = self.altitudes[sp_cesm][-1] + \
                                    (self.altitudes[sp_cesm][-1] - self.altitudes[sp_cesm][-2]) / 2.
                self.altitude_interfaces[sp_cesm][-2] = self.altitudes[sp_cesm][-2] + \
                                    (self.altitudes[sp_cesm][-2] - self.altitudes[sp_cesm][-3]) / 2.

                # calculate factor
                if self.Destination_unit in ['molecules cm-2 s-1', 'molecules/cm2/s']:
                    vert_factor = 1e5 # cm from km
                else:
                    raise ValueError( 'The destination unit is not currently not supported: ' + \
                                      self.Destination_unit )

                # calculate altitude intervals
                for vi in np.arange( len(self.altitudes[sp_cesm]) ):
                    self.altitude_intervals[sp_cesm][vi] = ( self.altitude_interfaces[sp_cesm][vi+1] - \
                                                             self.altitude_interfaces[sp_cesm][vi] ) * \
                                                            vert_factor

                altivar[:] = self.altitude_interfaces[sp_cesm]

                # set atrributes
                altvar.long_name = 'Altitude'
                altvar.units = 'km'
                altivar.long_name = 'Altitude interfaces'
                altivar.units = 'km'
            # ===== END Create Dimensions =====

            # ===== Create Variables =====
            self.vars_source = self.ds_source[sp_src].data_vars

            # Pass info variables
            for varname in self.vars_source:
                if varname.lower() in ['lon','lat','ncol','area','rrfac','date']:
                    if (varname == 'date') & (self.day_of_each_month != None):
                        continue
                    var = fid.createVariable( varname, 
                                              np.dtype( self.ds_source[sp_src][varname].values.flat[0] ).name,
                                              self.ds_source[sp_src][varname].dims )
                    var[:] = self.ds_source[sp_src][varname].values
                    for key in list( self.ds_source[sp_src][varname].attrs.keys() ):
                        var.setncattr( key, self.ds_source[sp_src][varname].attrs[key] )

            # add date variable if not available | or 
            if ('date' not in self.vars_source) | (self.day_of_each_month != None):
                # Retrieve time info
                if self.date_array == None:
                    if type( self.ds_source_time[sp_src]['time'].values[0] ) == np.datetime64:
                        self.time_year = self.ds_source_time[sp_src]['time'].values.astype('datetime64[Y]').astype('int') + 1970
                        self.time_month = self.ds_source_time[sp_src]['time'].values.astype('datetime64[M]').astype('int') % 12 + 1
                        self.time_day = ( self.ds_source_time[sp_src]['time'].values.astype('datetime64[D]') - \
                                          self.ds_source_time[sp_src]['time'].values.astype('datetime64[M]') + 1 ).astype('int')

                    elif type( self.ds_source_time[sp_src]['time'].values[0] ) in [cftime._cftime.DatetimeGregorian, 
                                                                                   cftime._cftime.DatetimeNoLeap]:
                        self.time_year = np.zeros( len(self.ds_source_time[sp_src]['time']) ).astype('int')
                        self.time_month = np.zeros( len(self.ds_source_time[sp_src]['time']) ).astype('int')
                        self.time_day = np.zeros( len(self.ds_source_time[sp_src]['time']) ).astype('int')

                        for ti, timeraw in enumerate( self.ds_source_time[sp_src]['time'].values ):
                            self.time_year[ti] = timeraw.year
                            self.time_month[ti] = timeraw.month
                            self.time_day[ti] = timeraw.day

                    else:
                        raise ValueError( 'Currently ' + str(type( self.ds_source_time[sp_src]['time'].values[0] )) \
                                        + ' is not supported! Check type of the time dimension' )
            
                    # Construct Date
                    self.date_array = []
                    for ti, ty in enumerate( self.time_year ):
                        if self.day_of_each_month == None:
                            self.date_array.append( self.time_year[ti]*10000 + self.time_month[ti]*100 + self.time_day[ti] )
                        else:
                            self.date_array.append( self.time_year[ti]*10000 + self.time_month[ti]*100 + self.day_of_each_month )

                else:
                    self.time_year = np.zeros( len(self.ds_source_time[sp_src]['time']) ).astype('int')
                    self.time_month = np.zeros( len(self.ds_source_time[sp_src]['time']) ).astype('int')
                    self.time_day = np.zeros( len(self.ds_source_time[sp_src]['time']) ).astype('int')

                    for dateii, dateval in enumerate(self.date_array):
                        self.time_year[dateii] = int( int(dateval) / 10000 )
                        self.time_month[dateii] = int( (int(dateval) - self.time_year[dateii] * 10000) / 100 )
                        self.time_day[dateii] = int(dateval) - self.time_year[dateii] * 10000 - self.time_month[dateii] * 100
                    
                # Add date variable to NetCDF file
                var = fid.createVariable( 'date', 'int32', ('time',) )
                var[:] = self.date_array
                var.units = 'YYYYMMDD'
                var.long_name = 'Date'



            # Calculation - Remappaing
            if self.dimensions[sp_cesm] == 'sfc':

                if self.keep_sector:
                    for sc in self.sectors[sp_cesm]:

                        # Calculation first
                        calc_str = 'self.var_value = '
                        for ii, sp_src in enumerate(self.eq_species[sp_cesm]):
                            # to consider missing sectors in some species
                            if sc in list( self.ds_source[sp_src].data_vars ):
                                sp_src_type = sp_src
                                self.convert_unit(sp_cesm, sp_src, sc)
                                calc_str += 'self.eq_frac["' + sp_cesm + '"][' + str(ii) + \
                                            '] * self.conversion_factor["' + sp_cesm + '"]["' + \
                                            sp_src + '"]'                         
                                calc_str += '* self.ds_source["' + sp_src + '"]["' + sc + '"].values '
                                if (self.N_equations[sp_cesm] > 1) & (ii + 1 < self.N_equations[sp_cesm]):
                                    calc_str += self.eq_operators[sp_cesm][ii]
                                
                        
                        if calc_str[-1] == '+':
                            calc_str = calc_str[:-1]
                        
                        exec(calc_str)

                        # Create Variable in NetCDF
                        var = fid.createVariable( sc,
                                                  np.dtype( self.ds_source[sp_src_type][sc].values.flat[0] ).name,
                                                  self.ds_source[sp_src_type][sc].dims )                
                        var[:] = self.var_value
                        var.standard_name = sp_cesm
                        var.long_name = 'Emissions of ' + sp_cesm + ' for ' + sc + ' sector'
                        if sp_cesm in self.diameter.keys():
                            var.units = '(particles cm-2 s-1)(molecules mole-1)(g kg-1)'
                        else:
                            var.units = self.Destination_unit
                        var.molecular_weight = self.mw_CESM[sp_cesm]
                        var.molecular_weight_units = 'g mole-1'

                        for key in list( self.ds_source[sp_src_type][sc].attrs.keys() ):
                            if key in ['standard_name', 'long_name', 'units',
                                       'molecular_weight', 'molecular_weight_units']:
                                continue
                            else:
                                var.setncattr( key, self.ds_source[sp_src_type][sc].attrs[key] )

                else:
                    # Calculation first
                    calc_str = 'self.var_value = '
                    str_sector = ''
                    for sc in self.sectors[sp_cesm]:

                        for ii, sp_src in enumerate(self.eq_species[sp_cesm]):
                            # to consider missing sectors in some species
                            if sc in list( self.ds_source[sp_src].data_vars ):
                                sp_src_type = sp_src
                                self.convert_unit(sp_cesm, sp_src, sc)
                                calc_str += 'self.eq_frac["' + sp_cesm + '"][' + str(ii) + \
                                            '] * self.conversion_factor["' + sp_cesm + '"]["' + \
                                            sp_src + '"]'                         
                                calc_str += '* self.ds_source["' + sp_src + '"]["' + sc + '"].values '
                                if (self.N_equations[sp_cesm] > 1) & (ii + 1 < self.N_equations[sp_cesm]):
                                    calc_str += self.eq_operators[sp_cesm][ii]

                        str_sector += sc
                        if not sc == self.sectors[sp_cesm][-1]:
                            calc_str += ' + '
                            str_sector += ' + '

                        # to consider missing sectors in some species
                        if calc_str[-2] == '+':
                            calc_str = calc_str[:-2]
                    exec(calc_str)

                    # Create Variable in NetCDF
                    var = fid.createVariable( 'emiss',
                                              np.dtype( self.ds_source[sp_src_type][sc].values.flat[0] ).name,
                                              self.ds_source[sp_src_type][sc].dims )                
                    var[:] = self.var_value
                    var.standard_name = sp_cesm
                    var.long_name = 'Emissions of ' + sp_cesm
                    if sp_cesm in self.diameter.keys():
                        var.units = '(particles cm-2 s-1)(molecules mole-1)(g kg-1)'
                    else:
                        var.units = self.Destination_unit
                    var.molecular_weight = self.mw_CESM[sp_cesm]
                    var.molecular_weight_units = 'g mole-1'
                    var.sectors = str_sector

            elif self.dimensions[sp_cesm] == 'vert':

                if self.keep_sector:
                    for sc in self.sectors[sp_cesm]:

                        # Calculation first for 2d dimension
                        calc_str = 'self.var_value_2d = '
                        for ii, sp_src in enumerate(self.eq_species[sp_cesm]):
                            
                            self.convert_unit(sp_cesm, sp_src, sc)
                            calc_str += 'self.eq_frac["' + sp_cesm + '"][' + str(ii) + \
                                        '] * self.conversion_factor["' + sp_cesm + '"]["' + \
                                        sp_src + '"]'                         
                            calc_str += '* self.ds_source["' + sp_src + '"]["' + sc + '"].values '
                            if (self.N_equations[sp_cesm] > 1) & (ii + 1 < self.N_equations[sp_cesm]):
                                calc_str += self.eq_operators[sp_cesm][ii]
                        exec(calc_str)

                        # === Distribute vertically ===
                        # Set array dimension first
                        src_dims = self.ds_source[sp_src][sc].dims

                        if 'lat' in src_dims:
                            vert_ind = np.where( np.char.equal( src_dims, 'lat' ) )[0][0]
                        elif 'ncol' in src_dims:
                            vert_ind = np.where( np.char.equal( src_dims, 'ncol' ) )[0][0]
                        else:
                            raise ValueError( 'Check dimension of source file!' + \
                                              ' It should have "lat/lon" or "ncol"' )
                        array_shape = []
                        array_shape_name = []
                        str_vert_index = '['
                        for di, dim in enumerate(src_dims):
                            if di == vert_ind:
                                array_shape_name.append( 'altitude' )
                                array_shape.append( len( self.altitudes[sp_cesm] ) )
                                str_vert_index += 'vi,'
                            array_shape_name.append( dim )
                            array_shape.append( self.ds_source[sp_src][dim].size )
                            if dim == src_dims[-1]:
                                str_vert_index += ':]'
                            else:
                                str_vert_index += ':,'

                        self.var_value_3d = np.zeros( array_shape )

                        # Redistribute 2d emissions to 3d
                        calc_3d_str = 'self.var_value_3d' + str_vert_index + ' = ' + \
                                      'self.var_value_2d * self.alt_fractions[sp_cesm][vi] / ' + \
                                      'self.altitude_intervals[sp_cesm][vi]'
                        for vi in np.arange( len(self.altitudes[sp_cesm]) ):
                            exec( calc_3d_str )


                        # Create Variable in NetCDF
                        var = fid.createVariable( sc,
                                                  np.dtype( self.ds_source[sp_src][sc].values.flat[0] ).name,
                                                  array_shape_name )
                        var[:] = self.var_value_3d
                        var.standard_name = sp_cesm
                        var.long_name = 'Emissions of ' + sp_cesm + ' for ' + sc + ' sector'
                        if sp_cesm in self.diameter.keys():
                            var.units = '(particles cm-3 s-1)(molecules mole-1)(g kg-1)'
                        else:
                            var.units = self.Destination_unit.replace('m-2','m-3')
                        var.molecular_weight = self.mw_CESM[sp_cesm]
                        var.molecular_weight_units = 'g mole-1'

                        for key in list( self.ds_source[sp_src][sc].attrs.keys() ):
                            if key in ['standard_name', 'long_name', 'units',
                                       'molecular_weight', 'molecular_weight_units']:
                                continue
                            else:
                                var.setncattr( key, self.ds_source[sp_src][sc].attrs[key] )

                else:

                    # Calculation first for 2d dimension
                    calc_str = 'self.var_value_2d = '
                    str_sector = ''
                    for sc in self.sectors[sp_cesm]:

                        for ii, sp_src in enumerate(self.eq_species[sp_cesm]):
                            self.convert_unit(sp_cesm, sp_src, sc)
                            calc_str += 'self.eq_frac["' + sp_cesm + '"][' + str(ii) + \
                                        '] * self.conversion_factor["' + sp_cesm + '"]["' + \
                                        sp_src + '"]'                         
                            calc_str += '* self.ds_source["' + sp_src + '"]["' + sc + '"].values '
                            if (self.N_equations[sp_cesm] > 1) & (ii + 1 < self.N_equations[sp_cesm]):
                                calc_str += self.eq_operators[sp_cesm][ii]

                        str_sector += sc
                        if not sc == self.sectors[sp_cesm][-1]:
                            calc_str += ' + '
                            str_sector += ' + '

                    exec(calc_str)

                    # === Distribute vertically ===
                    # Set array dimension first
                    src_dims = self.ds_source[sp_src][sc].dims

                    if 'lat' in src_dims:
                        vert_ind = np.where( np.char.equal( src_dims, 'lat' ) )[0][0]
                    elif 'ncol' in src_dims:
                        vert_ind = np.where( np.char.equal( src_dims, 'ncol' ) )[0][0]
                    else:
                        raise ValueError( 'Check dimension of source file!' + \
                                          ' It should have "lat/lon" or "ncol"' )
                    array_shape = []
                    array_shape_name = []
                    str_vert_index = '['
                    for di, dim in enumerate(src_dims):
                        if di == vert_ind:
                            array_shape_name.append( 'altitude' )
                            array_shape.append( len( self.altitudes[sp_cesm] ) )
                            str_vert_index += 'vi,'
                        array_shape_name.append( dim )
                        array_shape.append( self.ds_source[sp_src][dim].size )
                        if dim == src_dims[-1]:
                            str_vert_index += ':]'
                        else:
                            str_vert_index += ':,'

                    self.var_value_3d = np.zeros( array_shape )

                    # Redistribute 2d emissions to 3d
                    calc_3d_str = 'self.var_value_3d' + str_vert_index + ' = ' + \
                                  'self.var_value_2d * self.alt_fractions[sp_cesm][vi] / ' + \
                                  'self.altitude_intervals[sp_cesm][vi]'
                    for vi in np.arange( len(self.altitudes[sp_cesm]) ):
                        exec( calc_3d_str )


                    # Create Variable in NetCDF
                    var = fid.createVariable( 'emiss',
                                              np.dtype( self.ds_source[sp_src][sc].values.flat[0] ).name,
                                              array_shape_name )
                    var[:] = self.var_value_3d
                    var.standard_name = sp_cesm
                    var.long_name = 'Emissions of ' + sp_cesm
                    if sp_cesm in self.diameter.keys():
                        var.units = '(particles cm-3 s-1)(molecules mole-1)(g kg-1)'
                    else:
                        var.units = self.Destination_unit.replace('m-2','m-3')
                    var.molecular_weight = self.mw_CESM[sp_cesm]
                    var.molecular_weight_units = 'g mole-1'
                    var.sectors = str_sector

            elif self.dimensions[sp_cesm] == 'aircraft':

                if not self.ignore_warning:
                    print( 'Warning: some aircraft emissions have a unit of xx (c)m-3 s-1\n' + \
                           '         here we assume a unit of xx (c)m-2 s-1\n' + \
                           '         i.e. divide each layer emission by vertical depth')
                    
                if self.keep_sector:
                    for sc in self.sectors[sp_cesm]:

                        # Iteration vertically
                        calc_str = 'self.var_value = '
                        for ii, sp_src in enumerate(self.eq_species[sp_cesm]):
                            self.convert_unit(sp_cesm, sp_src, sc)
                            calc_str += 'self.eq_frac["' + sp_cesm + '"][' + str(ii) + \
                                        '] * self.conversion_factor["' + sp_cesm + '"]["' + \
                                        sp_src + '"]'                         
                            calc_str += '* self.ds_source["' + sp_src + '"]["' + sc + '"].values '
                            if (self.N_equations[sp_cesm] > 1) & (ii + 1 < self.N_equations[sp_cesm]):
                                calc_str += self.eq_operators[sp_cesm][ii]
                        exec(calc_str)

                        # === Divide by vertical depth ===
                        src_dims = self.ds_source[sp_src][sc].dims
                        
                        vert_dim_candidate = ['altitude','level']
                        vert_dim_check = True
                        for vdc in vert_dim_candidate:
                            if vdc in src_dims:
                                self.aircraft_vert_name = vdc
                                vert_dim_check = False
                                break
                        if vert_dim_check:
                            raise ValueError( 'Check vertical dimension of aircraft emissions!' + \
                                              ' It should have "altitude" or "level"' )
                        vert_ind = np.where( np.char.equal( src_dims, self.aircraft_vert_name ) )[0][0]

                        str_vert_index = '['
                        for di, dim in enumerate(src_dims):
                            if di == vert_ind:
                                str_vert_index += 'vi,'
                            elif dim == src_dims[-1]:
                                str_vert_index += ':]'
                            else:
                                str_vert_index += ':,'

                        self.var_value_final = np.zeros( self.var_value.shape )

                        # Divide by vertical depth
                        calc_3d_str = 'self.var_value_final' + str_vert_index + ' = ' + \
                                      'self.var_value' + str_vert_index +  ' / ' + \
                                      'self.altitude_intervals[sp_cesm][vi]'
                        for vi in np.arange( len(self.altitudes[sp_cesm]) ):
                            exec( calc_3d_str )


                        # Create Variable in NetCDF
                        if 'level' in src_dims:
                            src_dims_re = []
                            for sd in src_dims:
                                if sd == 'level':
                                    src_dims_re.append('altitude')
                                else:
                                    src_dims_re.append(sd)
                        else:
                            src_dims_re = src_dims
                        
                        var = fid.createVariable( sc,
                                                  np.dtype( self.ds_source[sp_src][sc].values.flat[0] ).name,
                                                  src_dims_re )
                        var[:] = self.var_value_final
                        var.standard_name = sp_cesm
                        var.long_name = 'Emissions of ' + sp_cesm + ' for ' + sc + ' sector'
                        if sp_cesm in self.diameter.keys():
                            var.units = '(particles cm-3 s-1)(molecules mole-1)(g kg-1)'
                        else:
                            var.units = self.Destination_unit.replace('m-2','m-3')
                        var.molecular_weight = self.mw_CESM[sp_cesm]
                        var.molecular_weight_units = 'g mole-1'

                        for key in list( self.ds_source[sp_src][sc].attrs.keys() ):
                            if key in ['standard_name', 'long_name', 'units',
                                       'molecular_weight', 'molecular_weight_units']:
                                continue
                            else:
                                var.setncattr( key, self.ds_source[sp_src][sc].attrs[key] )

                else:

                    # Interation vertically
                    calc_str = 'self.var_value = '
                    str_sector = ''
                    for sc in self.sectors[sp_cesm]:
                        
                        for ii, sp_src in enumerate(self.eq_species[sp_cesm]):
                            self.convert_unit(sp_cesm, sp_src, sc)
                            calc_str += 'self.eq_frac["' + sp_cesm + '"][' + str(ii) + \
                                        '] * self.conversion_factor["' + sp_cesm + '"]["' + \
                                        sp_src + '"]'                         
                            calc_str += '* self.ds_source["' + sp_src + '"]["' + sc + '"].values '
                            if (self.N_equations[sp_cesm] > 1) & (ii + 1 < self.N_equations[sp_cesm]):
                                calc_str += self.eq_operators[sp_cesm][ii]
                        
                        str_sector += sc
                        if not sc == self.sectors[sp_cesm][-1]:
                            calc_str += ' + '
                            str_sector += ' + '
                    
                    exec(calc_str)

                    # === Divide by vertical depth ===
                    src_dims = self.ds_source[sp_src][sc].dims

                    vert_dim_candidate = ['altitude','level']
                    vert_dim_check = True
                    for vdc in vert_dim_candidate:
                        if vdc in src_dims:
                            self.aircraft_vert_name = vdc
                            vert_dim_check = False
                            break
                    if vert_dim_check:
                        raise ValueError( 'Check vertical dimension of aircraft emissions!' + \
                                          ' It should have "altitude" or "level"' )
                    vert_ind = np.where( np.char.equal( src_dims, self.aircraft_vert_name ) )[0][0]

                    str_vert_index = '['
                    for di, dim in enumerate(src_dims):
                        if di == vert_ind:
                            str_vert_index += 'vi,'
                        elif dim == src_dims[-1]:
                            str_vert_index += ':]'
                        else:
                            str_vert_index += ':,'

                    self.var_value_final = np.zeros( self.var_value.shape )

                    # Divide by vertical depth
                    calc_3d_str = 'self.var_value_final' + str_vert_index + ' = ' + \
                                  'self.var_value' + str_vert_index +  ' / ' + \
                                  'self.altitude_intervals[sp_cesm][vi]'
                    for vi in np.arange( len(self.altitudes[sp_cesm]) ):
                        exec( calc_3d_str )
                        
                    # Create Variable in NetCDF
                    if 'level' in src_dims:
                        src_dims_re = []
                        for sd in src_dims:
                            if sd == 'level':
                                src_dims_re.append('altitude')
                            else:
                                src_dims_re.append(sd)
                    else:
                        src_dims_re = src_dims                    
                    
                    var = fid.createVariable( 'emiss',
                                              np.dtype( self.ds_source[sp_src][sc].values.flat[0] ).name,
                                              src_dims_re )
                    var[:] = self.var_value_final
                    var.standard_name = sp_cesm
                    var.long_name = 'Emissions of ' + sp_cesm 
                    if sp_cesm in self.diameter.keys():
                        var.units = '(particles cm-3 s-1)(molecules mole-1)(g kg-1)'
                    else:
                        var.units = self.Destination_unit.replace('m-2','m-3')
                    var.molecular_weight = self.mw_CESM[sp_cesm]
                    var.molecular_weight_units = 'g mole-1'
                    var.sectors = str_sector
                    
            else:
                raise ValueError( 'Check input file for dimensions! - sfc/vert/aircraft' )
            # ===== END Create Variables =====


            # ===== Set Global Attributes =====
            for key in list( self.ds_source[sp_src].attrs.keys() ):
                fid.setncattr( key, self.ds_source[sp_src].attrs[key] )
            fid.comment_remapping = '===== Below are created from remapping script ====='
            str_file_history = 'from file(s): '
            for sp_src in self.eq_species[sp_cesm]:
                if sp_src == self.eq_species[sp_cesm][-1]:
                    str_file_history += self.Source_files[sp_src]
                else:
                    str_file_history += self.Source_files[sp_src] + '; '    
            fid.remmaped_from_files = str_file_history
            fid.mapping_equation = self.eq_attribute[sp_cesm]
            fid.remapped_by = 'Species_Mapping.py tool'
            fid.remapping_time = str( datetime.datetime.now() )
            user_name = subprocess.getoutput( 'echo "$USER"')
            host_name = subprocess.getoutput( 'hostname -f' )
            fid.remmaping_username = user_name + ' on ' + host_name
            str_comment_from_Species_Mapping_script = ''
            
            if (self.N_equations[sp_cesm] > 1):
                str_comment_from_Species_Mapping_script += \
                    'Note that some global attribute information can be wrong as the script ' + \
                    'just copies global attributes from the last species file in the mapping equation.'
            if ( self.day_of_each_month != None ):
                if len(str_comment_from_Species_Mapping_script) > 5:
                    str_comment_from_Species_Mapping_script += '\n'
                str_comment_from_Species_Mapping_script += \
                    'day of each month of date array was fixed to 15th, ' + \
                    'which can be different from day of each month of time array. ' + \
                    'Time array is copied from original input files, ' + \
                    'and not foced by day_of_each_month keyword'
            fid.comment_from_Species_Mapping_script = str_comment_from_Species_Mapping_script
                
                
            # ===== END Set Global Attributes =====

            # Close NetCDF file
            fid.close()

        # print text for the use of user_nl_cam if user wants
        if self.print_user_nl_cam:
            print( 'Printing text for the use of user_nl_cam: ')
            date_now = datetime.datetime.now()
            YMD = str(date_now.year).zfill(4) + str(date_now.month).zfill(2) + \
                  str(date_now.day).zfill(2)
            print( 'Text is also saved in user_nl_cam_from_sp_map_c' + YMD )
            print( self.ext_frc_specifier )
            print( self.srf_emis_specifier )

            fid_nl = open( 'user_nl_cam_from_sp_map_c' + YMD, 'w' )
            fid_nl.write( self.ext_frc_specifier )
            fid_nl.write( self.srf_emis_specifier )
            fid_nl.close()            


    # convert units
    def convert_unit(self, sp_cesm, sp_src, sc):
        Avo = 6.022e23 # molecules / mole
        
        if sp_cesm not in self.conversion_factor.keys():
            self.conversion_factor[sp_cesm] = {}
        self.conversion_factor[sp_cesm][sp_src] = 1.
        
        conv_success = False
        if self.Source_unit in ['kg m-2 s-1', 'kg/m2/s']:
            if sp_cesm in self.diameter.keys():
                conv_success = True
                # kg particle-1
                self.mass_per_particle[sp_cesm] = self.density[sp_cesm] * (np.pi/6.) * \
                                                  self.diameter[sp_cesm]**(3)
                # unit: (particles/cm2/s)(molecules/mole)(g/kg)
                # MAM4 in CESM/CAM-chem has some weird unit for number
                self.conversion_factor[sp_cesm][sp_src] *= Avo * 1e3 * 1e-4 / \
                                                           self.mass_per_particle[sp_cesm] \
                                                          * self.mw_CESM[sp_cesm] \
                                                          / self.mw_Inventory[sp_src]
            else:
                if self.Destination_unit in ['molecules cm-2 s-1', 'molecules/cm2/s']:
                    conv_success = True
                    self.conversion_factor[sp_cesm][sp_src] *= \
                                                 Avo / (self.mw_Inventory[sp_src] / 1e3) / 1e4
                    if self.yield_type[sp_cesm] == 'mass':
                        self.conversion_factor[sp_cesm][sp_src] *= self.mw_Inventory[sp_src] / \
                                                                   self.mw_CESM[sp_cesm]
                        
        elif self.Source_unit in ['ton/month', 'ton month-1']:
            self.calc_area()
            self.conversion_factor[sp_cesm][sp_src] = np.zeros( self.ds_source[sp_src][sc].shape ) + 1
            
            if sp_cesm in self.diameter.keys():
                conv_success = True
                # kg particle-1
                self.mass_per_particle[sp_cesm] = self.density[sp_cesm] * (np.pi/6.) * \
                                                  self.diameter[sp_cesm]**(3)
                
                # unit: (particles/cm2/s)(molecules/mole)(g/kg)
                # MAM4 in CESM/CAM-chem has some weird unit for number
                for ti in np.arange( len( self.ds_source[sp_src]['time'] ) ):
                    DOM = calendar.monthrange( self.time_year[ti], self.time_month[ti] )[1]
                    self.conversion_factor[sp_cesm][sp_src][ti] *= Avo * 1e3 * 1e-4 / \
                                                               self.mass_per_particle[sp_cesm] \
                                                              * self.mw_CESM[sp_cesm] \
                                                              / self.mw_Inventory[sp_src] \
                                                              * 1000  / self.grid_area / ( DOM * 86400. )
                
            else:
                if self.Destination_unit in ['molecules cm-2 s-1', 'molecules/cm2/s']:
                    conv_success = True
                    
                    for ti in np.arange( len( self.ds_source[sp_src]['time'] ) ):
                        DOM = calendar.monthrange( self.time_year[ti], self.time_month[ti] )[1]
                        self.conversion_factor[sp_cesm][sp_src][ti] *= \
                                                     Avo / (self.mw_Inventory[sp_src] / 1e3) / 1e4 \
                                                     * 1000 / self.grid_area / (DOM * 86400.)
                    if self.yield_type[sp_cesm] == 'mass':
                        self.conversion_factor[sp_cesm][sp_src][ti] *= self.mw_Inventory[sp_src] / \
                                                                       self.mw_CESM[sp_cesm]


        if not conv_success:
            raise ValueError( 'The unit conversion from ', self.Source_unit, 
                              'to ', self.Destination_unit, ' is currently not supported')

    # calculate grid area
    def calc_area(self):
        
        Earth_rad = 6.371e6 # in m
        
        if self.grid_type == 'FV':
            if 'lat' not in self.dim_var.keys():
                raise ValueError( 'Check your dimension variables! ' + \
                                  'Latitude (lat) values are not available' )
            if 'lon' not in self.dim_var.keys():
                raise ValueError( 'Check your dimension variables! ' + \
                                  'Longitude (lon) values are not available' )              
            
            self.grid_area = np.zeros( ( len(self.dim_var['lat']), 
                                         len(self.dim_var['lon']) ) )
            
            for jj in np.arange( len(self.dim_var['lat']) ):
                if (jj == 0) & ( (self.dim_var['lat'][jj] + 90.) < 0.001 ):
                    self.dlat = self.dim_var['lat'][jj+1] - self.dim_var['lat'][jj]
                    self.dlon = self.dim_var['lon'][jj+1] - self.dim_var['lon'][jj]
                    sedge = ( self.dim_var['lat'][jj] ) * np.pi / 180.
                    nedge = ( self.dim_var['lat'][jj] + self.dlat / 2. ) * np.pi / 180.
                elif (jj == len(self.dim_var['lat'])-1) & ( (self.dim_var['lat'][jj] - 90.) < 0.001 ):
                    self.dlat = self.dim_var['lat'][jj] - self.dim_var['lat'][jj-1]
                    self.dlon = self.dim_var['lon'][jj] - self.dim_var['lon'][jj-1]
                    sedge = ( self.dim_var['lat'][jj] - self.dlat / 2. ) * np.pi / 180.
                    nedge = ( self.dim_var['lat'][jj] ) * np.pi / 180.
                else:
                    self.dlat = self.dim_var['lat'][jj+1] - self.dim_var['lat'][jj]
                    self.dlon = self.dim_var['lon'][jj+1] - self.dim_var['lon'][jj]                    
                    sedge = ( self.dim_var['lat'][jj] - self.dlat / 2. ) * np.pi / 180.
                    nedge = ( self.dim_var['lat'][jj] + self.dlat / 2. ) * np.pi / 180.
                    
                    
                self.grid_area[jj,:] = self.dlon * (np.pi/180.) * Earth_rad**(2) \
                                       * ( np.sin(nedge) - np.sin(sedge) )
            
        elif self.grid_type == 'SE':
            
            Earth_area = 4 * np.pi * Earth_rad**(2)
        
            self.grid_area = self.dim_var['grid_area_rad2'] / np.sum( self.dim_var['grid_area_rad2'] ) \
                            * Earth_area
            
    
    # ===== Defining __call__ method =====
    def __call__(self):
        print( '=== mapping file ===')
        print( self.mapping_file )
