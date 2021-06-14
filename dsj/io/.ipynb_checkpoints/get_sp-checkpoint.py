'''
get_sp.py
Extract species name from the file list by using the diff method

MODIFICATION HISTORY:
    Duseong Jo, 25, FEB, 2021: VERSION 1.00
    - Initial version
    Duseong Jo, 27, FEB, 2021: VERSION 2.00
    - Adding a option to ignore date
'''
import numpy as np


def get_sp(filelist, ignore_date=False):
    
    if ignore_date:
           
        filelist_f = []
        for fl in filelist:
            
            date_string = ''
            n_of_digit = 0
            digit_pos = []
            check_c = False
            check_underbar = False
            
            Right_Start_TMP = -1 
            while Right_Start_TMP > -1000:
                
                if (not check_underbar):
                    if fl[Right_Start_TMP].isdigit():
                        date_string = fl[Right_Start_TMP] + date_string
                        digit_pos.append( Right_Start_TMP )
                        n_of_digit += 1
                        if len(digit_pos) > 1.5:
                            for ii, dp in enumerate(digit_pos[:-1]):
                                if (dp - digit_pos[ii+1] -1) < 0.01 :
                                    continue
                                else:
                                    digit_pos = []
                                    n_of_digit = 0
                                    date_string = ''
                                    break

                    if n_of_digit > 0:
                        if fl[Right_Start_TMP] in ['c', 'C']:
                            date_string = fl[Right_Start_TMP] + date_string
                            check_c = True
                        if fl[Right_Start_TMP] == '_':
                            date_string = fl[Right_Start_TMP] + date_string
                            break
                
                Right_Start_TMP -= 1
                
            filelist_f.append( fl.replace(date_string,'') )
    
    else:
        filelist_f = filelist
        
        
        
    species = []
    species_files = {}
    for ii, fl in enumerate(filelist_f):
        Left_End = 1000
        Right_Start = -1000

        for ii2, fl2, in enumerate(filelist_f):
            if fl == fl2:
                continue
            else:
                strA = fl
                strB = fl2

                lenA = len(strA)
                lenB = len(strB)

                Left_End_TMP = 0
                while Left_End_TMP < 1000:
                    if strA[Left_End_TMP] == strB[Left_End_TMP]:
                        Left_End_TMP += 1
                        continue
                    else:
                        break

                Right_Start_TMP = -1
                while Right_Start_TMP > -1000:
                    
                    if strA[Right_Start_TMP] == strB[Right_Start_TMP]:
                        Right_Start_TMP -= 1
                        continue
                    else:
                        break
        
                Left_End = np.min([Left_End, Left_End_TMP])
                Right_Start = np.max([Right_Start, Right_Start_TMP])
        
        species.append( strA[Left_End:Right_Start+1] )
        species_files[species[ii]] = filelist[ii]

    return species, species_files