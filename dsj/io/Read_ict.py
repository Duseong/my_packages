'''
Read_ict.py
To read an icartt (ict) format file

MODIFICATION HISTORY:
    Duseong Jo, 02, JUL, 2021: VERSION 1.00
    Duseong Jo, 06, JUL, 2021: VERSION 2.00
    - Now can read multiple files at once
'''

from collections import OrderedDict
import numpy as np


class Read_ict(object):
    '''
    NAME:
           Read_ict
    
    PURPOSE:
           Read an icartt format (ict) file
           
    INPUTS:
           filename: ict filename

    RESOURCES:
           https://cdn.earthdata.nasa.gov/conduit/upload/6158/ESDS-RFC-029v2.pdf
    '''

    def __init__( self, filename ):
                       
        self.filename = filename
        
        
        if type(filename) == str:
            ict_fid = open( self.filename, 'r' )
        else:
            ict_fid = open( self.filename[0], 'r' )
        
        
        # First line
        First_line = ict_fid.readline().rstrip("\n")
        self.header_lines = int( First_line.split(',')[0] ) - 1
        
        # header information
        self.header = OrderedDict()
        if type(filename) != str:
            self.header['Number_of_Files'] = len( self.filename )
        
        self.header['PI_Name'] = ict_fid.readline().rstrip("\n") 
        self.header['PI_Affiliation'] = ict_fid.readline().rstrip("\n") 
        self.header['Data_Source_Description'] = ict_fid.readline().rstrip("\n")
        self.header['Mission_Name'] = ict_fid.readline().rstrip("\n")
        self.header['File_Volume_Number'], self.header['Total_Number_of_File_Volumes'] = ict_fid.readline().rstrip("\n").split(',')
        self.header['Data_Collection_and_Revision_Dates'] = ict_fid.readline().rstrip("\n")
        self.header['Data_Interval_Code'] = int( ict_fid.readline().rstrip("\n") )
        self.header['Independent_Variable_Definition'] = ict_fid.readline().rstrip("\n")
        self.header['Number_of_Dependent_Variables'] = int( ict_fid.readline().rstrip("\n") )
        self.header['Scale_Factors_of_Dependent_Variables'] = ict_fid.readline().rstrip("\n")
        self.header['Missing_Data_Flags_of_Dependent_Variables'] = ict_fid.readline().rstrip("\n")
        
        self.header['Variables'] = OrderedDict()
        self.header['Variables']['Short_Name'] = []
        self.header['Variables']['Unit'] = OrderedDict()
        self.header['Variables']['Standard_Name'] = OrderedDict()
        self.header['Variables']['Long_Name'] = OrderedDict()
        
        # Independent Variable
        TMP_vars = self.header['Independent_Variable_Definition'].split(',')
        self.header['Variables']['Short_Name'].append( TMP_vars[0] )
        if len(TMP_vars) >= 2:
            self.header['Variables']['Unit'][TMP_vars[0]] = TMP_vars[1]
        if len(TMP_vars) >= 3: # Some researchers don't report standard name
            self.header['Variables']['Standard_Name'][TMP_vars[0]] = TMP_vars[2] 
        if len(TMP_vars) == 4: # Long name is optional
            self.header['Variables']['Long_Name'][TMP_vars[0]] = TMP_vars[3]         
        
        # Dependent Variables      
        for i in np.arange( self.header['Number_of_Dependent_Variables'] ):
            TMP_vars = ict_fid.readline().rstrip("\n").split(',')
            if (len(TMP_vars) > 4):
                raise ValueError( "Check ict file! - Error at Depedent Variable Definitions" )
            self.header['Variables']['Short_Name'].append( TMP_vars[0] )
            if len(TMP_vars) >= 2:
                self.header['Variables']['Unit'][TMP_vars[0]] = TMP_vars[1]
            if len(TMP_vars) >= 3: # Some researchers don't report standard name
                self.header['Variables']['Standard_Name'][TMP_vars[0]] = TMP_vars[2] 
            if len(TMP_vars) == 4: # Long name is optional
                self.header['Variables']['Long_Name'][TMP_vars[0]] = TMP_vars[3] 
            
        # Special comments
        self.header['Number_of_Special_Comment_Lines'] = int( ict_fid.readline().rstrip("\n") )
        self.header['Special_Comments'] = ''
        for i in np.arange( self.header['Number_of_Special_Comment_Lines'] ):
            if i + 1 == self.header['Number_of_Special_Comment_Lines']:
                self.header['Special_Comments'] += ict_fid.readline().rstrip("\n")
            else:
                self.header['Special_Comments'] += ict_fid.readline().rstrip("\n") + '\n'
        
        # Normal comments
        self.header['Number_of_Normal_Comment_Lines'] = int( ict_fid.readline().rstrip("\n") )
        self.header['Normal_Comments'] = OrderedDict()
        self.header['Normal_Comments']['All'] = ''
        
        for i in np.arange( self.header['Number_of_Normal_Comment_Lines']-1 ):
            TMP_str = ict_fid.readline().rstrip("\n")
            if ':' in TMP_str:
                self.header['Normal_Comments'][TMP_str.split(':')[0]] = TMP_str.split(':')[1]
            else:
                print( "Warning: Some Normal Comments don't have KEYWORD:VALUE structure in " + self.filename )
                
            if i + 1 == self.header['Number_of_Normal_Comment_Lines']:
                self.header['Normal_Comments']['All'] += TMP_str
            else:
                self.header['Normal_Comments']['All'] += TMP_str + '\n'
        
        ict_fid.close()
        
        
        if type(filename) == str:
            self.data = np.genfromtxt( self.filename, delimiter=',', dtype=None, deletechars='',
                                       skip_header=self.header_lines, names=True )
        else:
            TMPdict = OrderedDict()
            TMPdict_type = {}

            
            for i, filen in enumerate(self.filename):
                tmpdata = np.genfromtxt( filen, delimiter=',', dtype=None, deletechars='',
                                         skip_header=self.header_lines, names=True )
                
                if i == 0:
                    for name in tmpdata.dtype.names:
                        TMPdict[name] = tmpdata[name]
                        TMPdict_type[name] = tmpdata[name].dtype  
                else:
                    for name in tmpdata.dtype.names:
                        TMPdict[name] = np.concatenate( (TMPdict[name], tmpdata[name] ) )
            
            self.data = TMPdict
            
        
        
    def __call__(self):
        if 'Number_of_Files' in self.header.keys():
            print( 'Number of Files: ' + str(self.header['Number_of_Files']) )
        print( 'PI Name: ' + self.header['PI_Name'] )
        print( 'PI Affiliation: ' + self.header['PI_Affiliation'] )
        print( 'Data Source Description: ' + self.header['Data_Source_Description'] )
        print( 'Mission Name: ' + self.header['Mission_Name'] )
        print( 'File Volume Number: ' + self.header['File_Volume_Number'] )
        print( 'Total Number of File Volumes: ' + self.header['Total_Number_of_File_Volumes'] )
        print( 'Data Collection and Revision Dates: ' + self.header['Data_Collection_and_Revision_Dates'] )
        print( 'Data Interval Code: ' + str(self.header['Data_Interval_Code']) )
        print( 'Independent Variable Definition: ' + self.header['Independent_Variable_Definition'] )          
        print( 'Number of Dependent Variables: ' + str(self.header['Number_of_Dependent_Variables']) )
        print( 'Scale Factors of Dependent Variables: ' + self.header['Scale_Factors_of_Dependent_Variables'] )
        print( 'Missing Data Flags of Dependent Variables: ' + self.header['Missing_Data_Flags_of_Dependent_Variables'] )
        print( 'Variables ====================' )
        for i, sname in enumerate( self.header['Variables']['Short_Name'] ):
            print( i, sname )
            try:
                print( '-- Unit: ', self.header['Variables']['Unit'][sname] )
            except:
                #print( '-- Unit: N/A' )
                1

            try:
                print( '-- Standard Name: ', self.header['Variables']['Standard_Name'][sname] )
            except:
                #print( '-- Unit: N/A' )
                1

            try:
                print( '-- Long Name: ', self.header['Variables']['Long Name'][sname] )
            except:
                #print( '-- Long Name: N/A' )
                1

        # Special Comments
        
        if self.header['Special_Comments'] == '':
            print( 'No Special Comments')
        else:
            print( 'Special Comments =====================')
            print( self.header['Special_Comments'] )

        # Normal Comments
        print( 'Normal Comments ======================' )
        print( self.header['Normal_Comments']['All'] )
