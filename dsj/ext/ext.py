'''
ext.py
Collection of some useful funcitons
(1) RMA (Reduced Major Axis)

MODIFICATION HISTORY:
    dsj, 16, MAY, 2015: VERSION 1.00
'''

import numpy as np


def rma(X, Y):
    """
    This code is obtained from
    http://pydoc.net/Python/oceans/0.2.1/oceans.ff_tools.teaching/
    Accessed May 16, 2015
    
    Calculate a "MODEL-2" least squares fit.
    
    The SLOPE of the line is determined by calculating the GEOMETRIC MEAN
    of the slopes from the regression of Y-on-X and X-on-Y.
    
    The equation of the line is:     y = mx + b.
    
    This line is called the GEOMETRIC MEAN or the REDUCED MAJOR AXIS.
    
    See Ricker (1973) Linear regressions in Fishery Research, J. Fish.
    Res. Board Can. 30: 409-434, for the derivation of the geometric
    mean regression.
    
    Since no statistical treatment exists for the estimation of the
    asymmetrical uncertainty limits for the geometric mean slope,
    I have used the symmetrical limits for a model I regression
    following Ricker's (1973) treatment.  For ease of computation,
    equations from Bevington and Robinson (1992) "Data Reduction and
    Error Analysis for the Physical Sciences, 2nd Ed."  pp: 104, and
    108-109, were used to calculate the symmetrical limits: sm and sb.
    
    Data are input and output as follows:
    
    m, b, r, sm, sb = lsqfitgm(X,Y)
    X    =    x data (vector)
    Y    =    y data (vector)
    m    =    slope
    b    =    y-intercept
    r    =    correlation coefficient
    sm   =    standard deviation of the slope
    sb   =    standard deviation of the y-intercept
    
    Note that the equation passes through the centroid:  (x-mean, y-mean)
        
    """
    
    X, Y = map(np.asanyarray, (X, Y))
    
    # Determine slope of Y-on-X regression.
    my = lsqfity(X, Y)[0]
    
    # Determine slope of X-on-Y regression.
    mx = lsqfitx(X, Y)[0]
    
    # Calculate geometric mean slope.
    m = np.sqrt(my * mx)
    
    if (my < 0) and (mx < 0):
        m = -m
        
    # Determine the size of the vector.
    n = len(X)
    
    # Calculate sums and means.
    Sx = np.sum(X)
    Sy = np.sum(Y)
    xbar = Sx / n
    ybar = Sy / n
    
    # Calculate geometric mean intercept.
    b = ybar - m * xbar
    
    # Calculate more sums.
    #Sxy = np.sum(X * Y)  # FIXME: Assigned but never used.
    Sx2 = np.sum(X ** 2)
    #Sy2 = np.sum(Y ** 2)  # FIXME: Assigned but never used.
    
    # Calculate re-used expressions.
    #num = n * Sxy - Sx * Sy  # FIXME: Assigned but never used.
    den = n * Sx2 - Sx ** 2
    
    # Calculate r, sm, sb and s2.    
    r = np.sqrt(my / mx)
    
    if (my < 0) and (mx < 0):
        r = -r
        
    diff = Y - b - m * X
    
    s2 = np.sum(diff * diff) / (n - 2)
    sm = np.sqrt(n * s2 / den)
    sb = np.sqrt(Sx2 * s2 / den)

    return {'S':m, 'Y':b, 'R':r, 'std_of_S':sm, 'std_of_Y':sb}


def lsqfity(X, Y):
    """
    Calculate a "MODEL-1" least squares fit.

    The line is fit by MINIMIZING the residuals in Y only.

    The equation of the line is:     Y = my * X + by.

    Equations are from Bevington & Robinson (1992)
    Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
    pp: 104, 108-109, 199.

    Data are input and output as follows:

    my, by, ry, smy, sby = lsqfity(X,Y)
    X     =    x data (vector)
    Y     =    y data (vector)
    my    =    slope
    by    =    y-intercept
    ry    =    correlation coefficient
    smy   =    standard deviation of the slope
    sby   =    standard deviation of the y-intercept

    """

    X, Y = map(np.asanyarray, (X, Y))

    # Determine the size of the vector.
    n = len(X)

    # Calculate the sums.

    Sx = np.sum(X)
    Sy = np.sum(Y)
    Sx2 = np.sum(X ** 2)
    Sxy = np.sum(X * Y)
    Sy2 = np.sum(Y ** 2)

    # Calculate re-used expressions.
    num = n * Sxy - Sx * Sy
    den = n * Sx2 - Sx ** 2

    # Calculate my, by, ry, s2, smy and sby.
    my = num / den
    by = (Sx2 * Sy - Sx * Sxy) / den
    ry = num / (np.sqrt(den) * np.sqrt(n * Sy2 - Sy ** 2))

    diff = Y - by - my * X

    s2 = np.sum(diff * diff) / (n - 2)
    smy = np.sqrt(n * s2 / den)
    sby = np.sqrt(Sx2 * s2 / den)

    return my, by, ry, smy, sby



def lsqfitx(X, Y):
    """Calculate a "MODEL-1" least squares fit.

    The line is fit by MINIMIZING the residuals in X only.

    The equation of the line is:     Y = mx * X + bx.

    Equations are modified from those in Bevington & Robinson (1992)
    Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
    pp: 104, 108-109, 199.

    Data are input and output as follows:

    mx, bx, rx, smx, sbx = lsqfitx(X, Y)
    X      =    x data (vector)
    Y      =    y data (vector)
    mx     =    slope
    bx     =    y-intercept
    rx     =    correlation coefficient
    smx    =    standard deviation of the slope
    sbx    =    standard deviation of the y-intercept

    """

    X, Y = map(np.asanyarray, (X, Y))

    # Determine the size of the vector.
    n = len(X)

    # Calculate the sums.
    Sx = np.sum(X)
    Sy = np.sum(Y)
    Sx2 = np.sum(X ** 2)
    Sxy = np.sum(X * Y)
    Sy2 = np.sum(Y ** 2)

    # Calculate re-used expressions.
    num = n * Sxy - Sy * Sx
    den = n * Sy2 - Sy ** 2

    # Calculate m, a, rx, s2, sm, and sb.
    mxi = num / den
    a = (Sy2 * Sx - Sy * Sxy) / den
    rx = num / (np.sqrt(den) * np.sqrt(n * Sx2 - Sx ** 2))

    diff = X - a - mxi * Y

    s2 = np.sum(diff * diff) / (n - 2)
    sm = np.sqrt(n * s2 / den)
    sa = np.sqrt(Sy2 * s2 / den)

    # Transpose coefficients
    mx = 1 / mxi
    bx = -a / mxi

    smx = mx * sm / mxi
    sbx = np.abs(sa / mxi)

    return mx, bx, rx, smx, sbx
