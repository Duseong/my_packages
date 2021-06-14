'''
VisIt_code_generator.py
this code is designed for generating VisIt codes for command line interface option
(1) 2D map plotting (class Map_2D)
(2) Quick plotting using class Map_2D

MODIFICATION HISTORY:
    dsj, 30, DEC, 2020: VERSION 1.00
      - 2D map plotting code generation
    dsj, 04, JAN, 2021: VERSION 2.00
      - Add Quick plot class
'''

import numpy as np
import datetime
from PIL import ImageFont
from dsj.io.NetCDF_IO import Write_NC
import os
from dsj.plot.map import get_cbar_prop


class Map_2D(object):
    
    '''
    NAME:
           Map_2D

    PURPOSE:
           Code generation for 2D map plotting

    INPUTS:
           filename: a file name which has values for plotting
           field_name: variable name in the file for plotting
           winsize: VisIt window size for plotting
           slicing: in case the variable is a 3D array with (X,Y,Z)
                    Currently, only Intercept method is supported
           slice_index: layer number on Z-axis to be plotted (31 is the SFC for CESM/CAM-chem 32L)
           drawlines: draw coastlines, country boundaries, state boundaries, etc. 
                      can be provided as either number or string
                     (0) no - no lines at all
                     (1) coast - coastlines
                     (2) country - country boundaries
                     (3) country_lake - country boundaries + lake (default)
                     (4) state - state boundaries
                     (5) state_lake - states boundaries + lake
           lines_resolution: can be 10m (fine), 50m (intermediate), or 110m (coarse)
                             default=10m
           lon_range: longitude ranges for 2D map
           lat_range: latitude ranges for 2D map
           show_grid: show grid
           plot_position: position of 2D map plot in the window [left, right, bottom, top]
           plot_min: minimum value of the plot
           plot_max: maximum value of the plot
           color_scale: can be 'Linear' or 'Log'
           color_table: color table name
           scale_factor: scaling factor for variable
           colorbar: Colorbar toggle (True or False)
           colorbar_position: position of colorbar in the window [left, bottom]
           colorbar_orientation: orientation of the colorbar
                     - VerticalRight, VerticalLeft, HozirontalTop, HorizontalBottom
           colorbar_scale_x: a scale factor for horizontal length of the colorbar
           colorbar_scale_y: a scale factor for vertical length of the colorbar
           colorbar_nticks: number of ticks in the colorbar
           colorbar_ticks: locations of ticks
           colorbar_labels: labels of ticks
           number_format: number format used in colorbar values
                          e.g., %5.2f, %1.1e, %3i
                          https://en.wikipedia.org/wiki/Printf_format_string#Type_field           
           fort_family: font family for texts in the plot (Arial, Courier, Times)
           unit: unit for the plot
           title: plot title
           output_filename: image filename to be saved
           noadd_date: in case you don't want to add the creation date at the end of the filename
           
           smoothing/grid-preserve?
    '''
    
    def __init__(self, filename='', field_name='', winsize=[1600,900], slicing=False, slice_index=31,
                       drawlines='country', lines_resolution='10m', lon_range=[0,358.75], lat_range=[-90, 90], 
                       show_grid=False, plot_position=[0.10, 0.95, 0.25, 0.90], plot_min=0, plot_max=None,
                       color_scale='Linear', color_table='myct_cont', scale_factor=1., colorbar=True, 
                       colorbar_position=[0.12,0.13], colorbar_orientation='HorizontalBottom', colorbar_scale_x=3, 
                       colorbar_scale_y=0.7, colorbar_nticks=5, colorbar_ticks=[], colorbar_labels=[], 
                       number_format='%7.2f', font_family='Times', title='', unit='', 
                       output_filename='', noadd_date=False ):
        
        self.filename = filename
        self.field_name = field_name
        self.winsize = winsize
        self.slicing = slicing
        self.slice_index = slice_index
        #self.drawlines = drawlines -> checked below
        self.lines_resolution = lines_resolution
        self.lon_range = lon_range
        self.lat_range = lat_range
        self.show_grid = show_grid
        self.plot_min = plot_min
        self.plot_max = plot_max
        self.color_scale = color_scale
        self.color_table = color_table
        self.scale_factor = scale_factor
        self.colorbar = colorbar
        self.plot_position = plot_position
        self.colorbar = colorbar
        self.colorbar_position = colorbar_position
        self.colorbar_orientation = colorbar_orientation
        self.colorbar_scale_x = colorbar_scale_x
        self.colorbar_scale_y = colorbar_scale_y
        self.colorbar_nticks = colorbar_nticks
        self.colorbar_ticks = colorbar_ticks
        self.colorbar_labels = colorbar_labels
        self.number_format = number_format
        self.font_family = font_family
        self.title = title
        self.unit = unit
        self.output_filename = output_filename
        self.noadd_date = noadd_date
        
        
        # ========================
        # ===== Check Errors =====
        # ========================
        # --- Filename ---
        if self.filename == '':
            raise ValueError( '"filename" must be provided!' )
        # --- longitude & latitude ranges ---
        if (self.lon_range[0] < -180.) or (self.lon_range[1] > 360.):
            print( 'Check "lon_range!"')
            raise ValueError( 'current values: ', lon_range )
        if (self.lat_range[0] < -90. ) or (self.lat_range[1] > 90.):
            print( 'Check "lat_range!"')
            raise ValueError( 'current values: ', lat_range )
        if (self.plot_position[0] < 0.) or (self.plot_position[2] < 0.) or \
           (self.plot_position[1] > 1.) or (self.plot_position[3] > 1.) or \
           (self.plot_position[1] <= self.plot_position[0] ) or \
           (self.plot_position[3] <= self.plot_position[2] ):
            print( 'Check "plot_position!"')
            raise ValueError( 'current values: ', plot_position )
        if self.colorbar_orientation not in ['VerticalRight', 'VerticalLeft', 
                                             'HozirontalTop', 'HorizontalBottom']:
            print( 'Check "colorbar_orientation!"')
            raise ValueError( 'current values: ', colorbar_orientation )
        if (self.colorbar_ticks != []) & (len(self.colorbar_ticks) != self.colorbar_nticks):
            print( 'number of elements in colorbar_ticks must be the same as colorbar_nticks')
            raise ValueError( 'current values: ', colorbar_ticks )
        if (self.colorbar_labels != []) & (len(self.colorbar_labels) != self.colorbar_nticks):
            print( 'number of elements in colorbar_labels must be the same as colorbar_nticks')
            raise ValueError( 'current values: ', colorbar_labels )
        # ===== END Check Errors =====
        
        
        
        # =============================
        # ===== Ininital Setting ======
        # =============================
        # --- longitude range check ---
        if np.min( lon_range ) < 0:
            self.lon_range_flag = 0
        else:
            self.lon_range_flag = 1
        # --- convert drawlines number to string ---
        drawlines_conv_dict_num_str = {0:'no',
                                       1:'coast',
                                       2:'country',
                                       3:'country_lake',
                                       4:'state',
                                       5:'state_lake'}
        if type(drawlines) == int:
            self.drawlines = drawlines_conv_dict_num_str[drawlines]
        elif type(drawlines) == str:
            self.drawlines = drawlines
        else:
            print( 'Check drawlines!' )
            raise ValueError( 'It must be string or integer number' )
        if self.drawlines != 'no': # Skip below if drawlines="no"
            # --- drawlines file location ---
            shp_base_dir = '/glade/u/home/cdswk/python/Shape_files/NaturalEarth/'
            shp_files = { 'coast':'RES_physical/ne_RES_coastline.shp',
                          'country':'RES_cultural/ne_RES_admin_0_countries.shp',
                          'country_lake':'RES_cultural/ne_RES_admin_0_countries_lakes.shp',
                          'state':'RES_cultural/ne_RES_admin_1_states_provinces.shp',
                          'state_lake':'RES_cultural/ne_RES_admin_1_states_provinces_lakes.shp',
                          'bounding_box':'RES_physical/ne_RES_wgs84_bounding_box.shp' }       
            self.drawlines_file = shp_base_dir + shp_files[drawlines].replace('RES',lines_resolution)
            if (self.lon_range_flag):
                self.drawlines_file = self.drawlines_file.replace('.shp','_0_360_processed.shp')
            else:
                self.drawlines_file = self.drawlines_file
        # --- Construct output filename ---
        if self.output_filename == '':
            self.save_image = False
        else:
            self.save_image = True
            if not noadd_date:
                CTime = datetime.datetime.now()
                CDate = str(CTime.year).zfill(4) + str(CTime.month).zfill(2) + \
                        str(CTime.day).zfill(2)
                self.CDate = CDate
                self.output_filename += '_c' + CDate
            self.output_filename += '_'
        # ===== END Initial Setting =====
        
        self.code_gen()



    def code_gen(self):

        # ===== Check list =====
        print( "# ===== Check list =====")
        print( '# - Options -> Appearance -> uncheck use default system')
        # ===== END Check list =====
        
        # Buffer zone
        print( "" )
        print( "" )
        print("# ▼▼▼▼▼▼▼ COPY BELOW ▼▼▼▼▼▼▼")
        print( "" )
        print( "" )
        
        # ==============================================
        # ===== print code for VisIt CLI interface =====
        # ==============================================
        
        # ===== Basic plotting =====
        print( '# === Resize Window === ' )
        print( 'ResizeWindow(1,' + str(self.winsize[0]) + ',' + str(self.winsize[1]) + ')' )
        print( '# === Read files === ' )
        print( 'OpenDatabase("localhost:' + self.filename + '"' + ', 0)' )
        if self.drawlines != 'no':
            print( 'OpenDatabase("localhost:' + self.drawlines_file + '"' + ', 0)' )
        print( '# === Draw plots === ' )
        print( 'ActivateDatabase("localhost:' + self.filename + '")')
        print( 'DefineScalarExpression("' + self.field_name + '_scaled", "' + self.field_name + ' * ' \
                                          + str(self.scale_factor) + '")')
        print( 'AddPlot("Pseudocolor", "' + self.field_name + '_scaled' + '", 1, 0)' )
        if self.slicing:
            print( 'AddOperator("Slice", 0)' )
            print( 'SetActivePlots(0)' )
            print( 'SliceAtts = SliceAttributes()' )
            print( 'SliceAtts.originType = SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node' )
            print( 'SliceAtts.originPoint = (0, 0, 0)' )
            print( 'SliceAtts.originIntercept = ' + str(self.slice_index) )
            print( 'SliceAtts.originPercent = 0' )
            print( 'SliceAtts.originZone = 0' )
            print( 'SliceAtts.originNode = 0' )
            print( 'SliceAtts.normal = (0, 0, 1)' )
            print( 'SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi' )
            print( 'SliceAtts.upAxis = (0, 1, 0)' )
            print( 'SliceAtts.project2d = 1' )
            print( 'SliceAtts.interactive = 1' )
            print( 'SliceAtts.flip = 0' )
            print( 'SliceAtts.originZoneDomain = 0' )
            print( 'SliceAtts.originNodeDomain = 0' )
            print( 'SliceAtts.theta = 0' )
            print( 'SliceAtts.phi = 90' )
            print( 'SetOperatorOptions(SliceAtts, 0, 0)' )
        if self.drawlines != 'no':
            print( 'ActivateDatabase("localhost:' + self.drawlines_file + '")')
            print( 'AddPlot("Mesh", "polygon", 1, 0)' )
        # ===== END Basic plotting =====
        
        
        # ===== Adjustplot and legend, remove redundant info. ===== 
        print( '# === Plot Position ===' )
        print( 'View2DAtts = View2DAttributes()' )
        print( 'View2DAtts.viewportCoords = (' + str(self.plot_position)[1:-1] + ')' )
        print( 'View2DAtts.windowCoords = (' + str(self.lon_range[0]) + ',' + str(self.lon_range[1]) + ',' + \
                                               str(self.lat_range[0]) + ',' + str(self.lat_range[1]) + ' )' )
        print( 'SetView2D(View2DAtts)' )
        print( '# === Toggle information about user, database, and time ===' )
        print( 'AnnotationAtts = AnnotationAttributes()' )
        print( 'AnnotationAtts.userInfoFlag = 0' )
        print( 'AnnotationAtts.databaseInfoFlag = 0' )
        print( 'AnnotationAtts.timeInfoFlag = 0' )
        # --- Longitude & Latitude ---
        print( '# === Longitude & Latitude adjustment ===' )
        print ('# --- Longitude Title ---')
        print( 'AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.' + self.font_family )
        print( 'AnnotationAtts.axes2D.xAxis.title.font.scale = 2' )
        print( 'AnnotationAtts.axes2D.xAxis.title.font.bold = 1' )
        print( 'AnnotationAtts.axes2D.xAxis.title.font.italic = 1' )
        print( 'AnnotationAtts.axes2D.xAxis.title.userTitle = 1' ) ### Can be customized in the future 
        print( 'AnnotationAtts.axes2D.xAxis.title.userUnits = 1' ) ### if needed
        print( 'AnnotationAtts.axes2D.xAxis.title.title = "Longitude"' )
        print( 'AnnotationAtts.axes2D.xAxis.title.units = "" # Unit - º' ) # º
        print ('# --- Longitude Label ---')
        print( 'AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.' + self.font_family )
        print( 'AnnotationAtts.axes2D.xAxis.label.font.scale = 1.5' )
        print ('# --- Latitude Title ---')
        print( 'AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.' + self.font_family )
        print( 'AnnotationAtts.axes2D.yAxis.title.font.scale = 2' )
        print( 'AnnotationAtts.axes2D.yAxis.title.font.bold = 1' )
        print( 'AnnotationAtts.axes2D.yAxis.title.font.italic = 1' )
        print( 'AnnotationAtts.axes2D.yAxis.title.userTitle = 1' ) ### Can be customized in the future 
        print( 'AnnotationAtts.axes2D.yAxis.title.userUnits = 1' ) ### if needed
        print( 'AnnotationAtts.axes2D.yAxis.title.title = "Latitude"' )
        print( 'AnnotationAtts.axes2D.yAxis.title.units = "" # Unit - º' ) 
        print ('# --- Latitude Label ---')
        print( 'AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.' + self.font_family )
        print( 'AnnotationAtts.axes2D.yAxis.label.font.scale = 1.5' )
        if self.show_grid:
            print( '# --- Longitude and Latitude grids ---')
            print( 'AnnotationAtts.axes2D.xAxis.grid = 1')
            print( 'AnnotationAtts.axes2D.yAxis.grid = 1')
        
        if self.colorbar == False:
            print( 'AnnotationAtts.legendInfoFlag = 0' )
        
        print( 'SetAnnotationAttributes(AnnotationAtts)' )
        # ===== END Adjustplot and legend, remove redundant info. ===== 

        
        # ===== Plot color attributes =====
        print( '# === Plot color attributes ===')
        print( 'PseudocolorAtts = PseudocolorAttributes()' )
        print( 'PseudocolorAtts.scaling = PseudocolorAtts.' + self.color_scale )
        if self.plot_min != None:
            print( 'PseudocolorAtts.minFlag = 1' )
            print( 'PseudocolorAtts.min = ' + str(self.plot_min) )
        if self.plot_max != None:
            print( 'PseudocolorAtts.maxFlag = 1' )
            print( 'PseudocolorAtts.max = ' + str(self.plot_max) )
        print( 'PseudocolorAtts.colorTableName = "' + self.color_table + '"' )
        print( 'SetPlotOptions(PseudocolorAtts)' )
        # ===== END Plot color attributes =====
        
        
        # ===== Colorbar properties =====
        if self.colorbar:
            print( '# === Colorbar properties ===')
            print( 'plotName = GetPlotList().GetPlots(0).plotName' )
            print( 'legend = GetAnnotationObject(plotName)' )
            print( '# === Scale the legend ===' )
            print( 'legend.xScale = ' + str(self.colorbar_scale_x) )
            print( 'legend.yScale = ' + str(self.colorbar_scale_y) )
            print( '# === Draw Bounding box? ===' )
            print( 'legend.drawBoundingBox = 0' )
            print( '#legend.boundingBoxColor = (180,180,180,230)' )
            print( '# === Make it horizontal ===')
            print( 'legend.orientation = legend.' + str(self.colorbar_orientation) )
            print( '# === Moving the legend ===' )
            print( 'legend.managePosition = 0' )
            print( 'legend.position = (' + str(self.colorbar_position[0]) + ',' + \
                                           str(self.colorbar_position[1]) + ') # upper left point' )
            print( '# === text color ===' )
            print( '#InvertBackgroundColor() # For black background' )
            print( 'legend.useForegroundForTextColor = 0' )
            print( 'legend.textColor = (0, 0, 0, 0) # Black' )
            print( '# === Number format ===' )
            print( 'legend.numberFormat = "' + self.number_format + '"' )
            print( '# === Font ===' )
            print( 'legend.fontFamily = legend.' + self.font_family )
            print( 'legend.fontBold = 1' )
            print( 'legend.fontItalic = 1' )
            print( 'legend.fontHeight = 0.04 # Font size' )
            print( '# === Labels, min/max info., title, etc. ===' )
            print( 'legend.numTicks = ' + str(self.colorbar_nticks) )
            print( 'legend.drawMinMax = 0 # print min and max value' )
            print( 'legend.minMaxInclusive = 1' ) ### what is this?
            print( 'legend.drawTitle = 0 # Turning off the title' )           
            # Use user-supplied labels
            if self.colorbar_ticks != []:
                print ('# === User-supplied Labels ===')
                print( 'legend.controlTicks=0' )
                print( 'legend.suppliedValues = (' + str(self.colorbar_ticks)[1:-1] + ')' )
                if self.colorbar_labels != []:
                    print( 'legend.drawLabels = legend.Labels # Turn on labels - None, Values, Labels, Both' )
                    print( 'legend.suppliedLabels = (' + str(self.colorbar_labels)[1:-1] + ')' )
        # ===== END Colorbar properties =====
        
        
        # ===== Suppress colorbar properties of Mesh =====
        if self.drawlines != 'no':
            print( '# === Colorbar properties of Mesh ===')
            print( 'plotName = GetPlotList().GetPlots(1).plotName' )
            print( 'legend = GetAnnotationObject(plotName)' )
            print( 'legend.drawTitle = 0 # Turning off the title' )
        # ===== END Suppress colorbar properties of Mesh =====
        
        
        # ===== Add title =====
        if self.title != '':
            self.calc_title_position()
            print( '# === Add title ===' )
            print( 'TitleObj = CreateAnnotationObject("Text2D")' )
            print( 'TitleObj.textColor = (0, 0, 0, 255)' )
            print( 'TitleObj.height = 0.04' )
            print( 'TitleObj.position = (' + str(self.left_position_title) + ',' \
                                           + str(self.plot_position[3])  + ')' )
            print( 'TitleObj.text = "' + self.title + '"' )
            print( 'TitleObj.fontFamily = TitleObj.' + self.font_family )
            print( 'TitleObj.fontBold = 1' )
            print( 'TitleObj.fontItalic = 0' )
        # ===== END Add title =====

        
        # ===== Add Unit =====
        if self.unit != '':
            print( '# === Add unit ===' )
            if self.unit == 'ugm-3':
                print( 'UnitObj = CreateAnnotationObject("Image")' )
                print( 'UnitObj.transparencyColor = (255,255,255,255)' )
                print( 'UnitObj.useTransparencyColor = 1' )
                print( 'UnitObj.image = "/glade/work/cdswk/Visualization/Image_files/Unit_ugm3_Times_Italic.png"' )
                print( 'UnitObj.position = (0.915, 0.02)' )
            else:
                print( 'UnitObj = CreateAnnotationObject("Text2D")' )
                print( 'UnitObj.textColor = (0, 0, 0, 255)' )
                print( 'UnitObj.height = 0.03' )             
                print( 'UnitObj.position = (0.93, 0.05)' )
                print( 'UnitObj.text = "' + self.unit + '"' )
                print( 'UnitObj.fontFamily = UnitObj.' + self.font_family )
                print( 'UnitObj.fontBold = 0' )
                print( 'UnitObj.fontItalic = 1' )
        # ===== END Add title =====
        
        
        # ===== Draw plot =====
        print( 'DrawPlots()' )
        
        
        # ===== Save image to file =====
        if self.save_image:
            print( 'SWA = SaveWindowAttributes()' ) 
            print( 'SWA.format = SWA.PNG' )
            print( 'SWA.fileName = "' + self.output_filename + '"' )
            print( 'SWA.screenCapture = 1' )
            print( 'SetSaveWindowAttributes(SWA)' )
            print( 'name = SaveWindow()' )
        # ===== END Save image to file =====
        
        # ===== END print code for VisIt CLI interface =====
        
        
    # Calculate Title Position (Left position) for plot    
    def calc_title_position(self):
        
        font = ImageFont.truetype("/glade/u/home/cdswk/python/miniconda3/lib/python3.7/site-packages" + \
                                  "/matplotlib/mpl-data/fonts/ttf/STIXGeneralBol.ttf")
        length = font.getsize( self.title )[0]
        
        left_margin = self.plot_position[0]
        plot_half = ( self.plot_position[1] - self.plot_position[0] ) / 2.
        
        length_full = 397.
        string_half = length / length_full / 2.
        self.left_position_title = left_margin + plot_half - string_half

    # ===== Defining __call__ method =====
    def __call__(self):
        print( '=== filename ===')
        print( self.filename )
        
        print( '=== lon_ranges_flag ===')
        print( self.lon_ranges_flag )
        
        print( '=== drawlines_file ===')
        print( self.drawlines_file )
   
        print( '=== left_position_title ===')
        print( self.left_position_title )
        
        
        
class Quick_plot(object):
    '''
    NAME:
           Quick_plot
           
    PURPOSE:
           Quick code generation for 2D map plotting
           
    INPUTS:
           array: an array including values for 2D map plot
           dimindices: in case you want to extract a portion of xarray variable (e.g. [0,-1,:,:])
                       it can be a dictionary for different dimensions for different variables
           index: vertical index when 3D fields are provided
           name: variable name if you want to provide the variable name for plot title
           lon: longitude values for non-xarray case
           lat: latitude values for non-xarray case
           scalefactor: scale factor for plotting variable
    '''
    
    def __init__(self, array, dimindices='[:]', index=None, name='', lon=None, lat=None, scalefactor=1e9 ):
        
        # === Get keywords ===
        self.array = array
        self.dimindices = dimindices
        self.index = index
        self.name = name
        self.lon = lon
        self.lat = lat
        self.scalefactor = scalefactor
        # === END Get keywords ===
        
        
        # === Check if array is xarray variable first ===
        if type( self.array ) in [ xr.core.dataset.Dataset, 
                                   xr.core.dataarray.DataArray ]:
            self.xarray_flag = True
            if self.name == '':
                raise ValueError( 'name must be provided for xarray')
        else:
            self.xarray_flag = False
            if (self.lons==None):
                raise ValueError( 'lons must be provided for non-xarray' )
            if (self.lats==None):
                raise ValueError( 'lats must be provided for non-xarray' )            
        # === END Check if array is xarray variable first ===
        
        
        # === Setup dimension variables and check errors ===
        if not (self.xarray_flag):
            if np.ndim( self.array ) == 2:
                dnames = ['lat','lon']
            elif np.ndim( self.array ) == 3:
                if index == None:
                    raise ValueError( 'index must be provided when dimension of array is 3 for non-xarray' )
                else:
                    dnames = ['lev','lat','lon']
            else:
                raise ValueError( 'array dimension must be 3 or less for non-xarray' )
        # ===END Setup dimension variables and check errors ===      
        
        
        # === check name
        if self.name == '':
            self.name='TMP'
        
        
        # === Set filename ===
        base_dir = '/glade/scratch/cdswk/temp/'
        CTime = datetime.datetime.now()
        CDate = str(CTime.year).zfill(4) + str(CTime.month).zfill(2) + \
                str(CTime.day).zfill(2)
        self.cdate = CDate
        for number in np.arange(99):
            self.filename = base_dir + 'TMP_file_for_VisIt_2D_plot_' + name + '_c_' + CDate + \
                            '_' + str(number).zfill(2) + '.nc'
            if os.path.exists( self.filename ):
                continue
            else:
                break      
        # === END Set filename ===
        
        
        # === Save NetCDF file ===
        if self.xarray_flag:
            nctmp = Write_NC( vnames=[self.name], values=self.array, filename=self.filename, 
                              dimindices=self.dimindices, noadd_date=True )
        else:
            values_dict = {name:self.array}
            values_dict['lon'] = self.lon
            values_dict['lat'] = self.lat
            if np.ndim(self.array) == 3:
                self.lev = np.arange( len(array[:,0,0]) )
                values_dict['lev'] = self.lev
            nctmp = Write_NC( vnames=[self.name], dnames=dnames, vdims={name:(str(dnames[1:-1]))}, 
                              gatts={'comments':'created for temporary use for 2D ploting'},
                              values=values_dict, datatype='f8', filename=self.filename, noadd_date=True )
        # === END Save NetCDF file ===
    
        
        # === Setup colorbar ticks ===
        #if self.index==None:
        #    cbprop = get_cbar_prop( [self.array], min_set=0, Nticks_set=5 )
        #else:
        #    cbprop = get_cbar_prop( [self.array[self.index,:,:]], min_set=0, Nticks_set=5 )        
        # === END Setup colorbar ticks ===
    
    
        # === Generate VisIt code ===
        if self.index==None:
            VTMP = Map_2D( filename=self.filename, field_name=name, slicing=False, 
                           title=self.name, color_scale='Linear', output_filename=self.filename[:-3], 
                           scale_factor=self.scalefactor, noadd_date=True )
        else:
            VTMP = Map_2D( filename=self.filename, field_name=name, slicing=True, slice_index=self.index,
                           title=self.name, color_scale='Linear', output_filename=self.filename[:-3], 
                           scale_factor=self.scalefactor, noadd_date=True ) 
        # === END Generate VisIt code ===
    