'''
expand_date.py
This python code changes YYYYMMDDhh format to real time value

MODIFICATION HISTORY:
    dsj, 07, JUL, 2015: VERSION 1.00
    dsj, 17, MAY, 2019: VERSION 2.00
                      - updated for python 3.7
'''
import calendar



class Expand_date(object):
    '''
    NAME:
           Expand_date

    PURPOSE:
           Changes YYYYMMDDhh format to real time value

    INPUTS:
           string -> base string for file name lists
           starttime -> starttime in [ year, month, day, hour ]
           endtime -> endtime in [ year, month, day, hour ]
    '''
    
    def __init__(self, string, start=[], end=[], interval=None):


        if len(start) != 4:
            if len(start) == 1:
                start.append(0o1) 
            if len(start) == 2:
                start.append(0o1)
            if len(start) == 3:
                start.append(00)

        if len(end) != 4:
            if len(end) == 1:
                end.append(0o1)
            if len(end) == 2:
                end.append(0o1)
            if len(end) == 3:
                end.append(00)
                
        self.start = start
        self.end = end
        self.files = []

        # Error check for start & end time and interval
        if len(start) == 0:
            raise ValueError( 'Check start! \n' +
                              'start should be have' +
                              '[year,month,day,hour]' )
        if len(end) == 0:
            raise ValueError( 'Check end! \n' +
                              'end should be have' +
                              '[year,month,day,hour]' )        
        if interval == None:
            raise ValueError( 'Please put time interval value!' )


        # get start s
        Year = start[0]
        Month = start[1]
        Day = start[2]
        Hour = start[3]        

        # make file lists according to time interval (first step)
        TMP = string.replace( 'YYYY', str(Year).zfill(4) )
        TMP = TMP.replace( 'MM', str(Month).zfill(2) )
        TMP = TMP.replace( 'DD', str(Day).zfill(2) )
        TMP = TMP.replace( 'hh', str(Hour).zfill(2) )
        self.files.append( TMP )

            
        while Year != end[0] or Month != end[1] or \
        Day != end[2] or Hour != end[3]:
            
            if interval == 'Y':
                Year = Year+1
            elif interval == 'M':
                Month = Month+1
                if Month == 13:
                    Month = 1
                    Year = Year+1
            elif interval == 'D':
                Day = Day+1
                DOM = calendar.monthrange(Year,Month)[1]
                if Day > DOM:
                    Month = Month+1
                    Day = 1
                    if Month == 13:
                        Month = 1
                        Year = Year+1                
            elif interval == 'H':
                Hour = Hour+1
                if Hour == 24:
                    Hour = 0
                    Day = Day + 1
                    DOM = calendar.monthrange(Year,Month)[1]
                    if Day > DOM:
                        Month = Month+1
                        Day = 1
                        if Month == 13:
                            Month = 1
                            Year = Year+1 
                            
            # make file lists according to time interval
            TMP = string.replace( 'YYYY', str(Year).zfill(4) )
            TMP = TMP.replace( 'MM', str(Month).zfill(2) )
            TMP = TMP.replace( 'DD', str(Day).zfill(2) )
            TMP = TMP.replace( 'hh', str(Hour).zfill(2) )
            self.files.append( TMP )





'''
for python 2.7
class Expand_date(object):
'''
'''
    NAME:
           Expand_date

    PURPOSE:
           Changes YYYYMMDDhh format to real time value

    INPUTS:
           string -> base string for file name lists
           starttime -> starttime in [ year, month, day, hour ]
           endtime -> endtime in [ year, month, day, hour ]
'''
'''
    def __init__(self, string, start=[], end=[], interval=None):


        if len(start) != 4:
            if len(start) == 1:
                start.append(01)
            if len(start) == 2:
                start.append(01)
            if len(start) == 3:
                start.append(00)

        if len(end) != 4:
            if len(end) == 1:
                end.append(01)
            if len(end) == 2:
                end.append(01)
            if len(end) == 3:
                end.append(00)
                
        self.start = start
        self.end = end
        self.files = []

        # Error check for start & end time and interval
        if len(start) == 0:
            raise ValueError( 'Check start! \n' +
                              'start should be have' +
                              '[year,month,day,hour]' )
        if len(end) == 0:
            raise ValueError( 'Check end! \n' +
                              'end should be have' +
                              '[year,month,day,hour]' )        
        if interval == None:
            raise ValueError( 'Please put time interval value!' )


        # get start s
        Year = start[0]
        Month = start[1]
        Day = start[2]
        Hour = start[3]        

        # make file lists according to time interval (first step)
        TMP = string.replace( 'YYYY', str(Year).zfill(4) )
        TMP = TMP.replace( 'MM', str(Month).zfill(2) )
        TMP = TMP.replace( 'DD', str(Day).zfill(2) )
        TMP = TMP.replace( 'hh', str(Hour).zfill(2) )
        self.files.append( TMP )

            
        while Year != end[0] or Month != end[1] or \
        Day != end[2] or Hour != end[3]:
            
            if interval == 'Y':
                Year = Year+1
            elif interval == 'M':
                Month = Month+1
                if Month == 13:
                    Month = 1
                    Year = Year+1
            elif interval == 'D':
                Day = Day+1
                DOM = calendar.monthrange(Year,Month)[1]
                if Day > DOM:
                    Month = Month+1
                    Day = 1
                    if Month == 13:
                        Month = 1
                        Year = Year+1                
            elif interval == 'H':
                Hour = Hour+1
                if Hour == 24:
                    Hour = 0
                    Day = Day + 1
                    DOM = calendar.monthrange(Year,Month)[1]
                    if Day > DOM:
                        Month = Month+1
                        Day = 1
                        if Month == 13:
                            Month = 1
                            Year = Year+1 
                            
            # make file lists according to time interval
            TMP = string.replace( 'YYYY', str(Year).zfill(4) )
            TMP = TMP.replace( 'MM', str(Month).zfill(2) )
            TMP = TMP.replace( 'DD', str(Day).zfill(2) )
            TMP = TMP.replace( 'hh', str(Hour).zfill(2) )
            self.files.append( TMP )

'''
