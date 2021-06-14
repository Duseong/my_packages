'''
bilinear.py
This code is made for bilinear interpolation

MODIFICATION HISTORY:
    dsj, 12, AUG, 2019: VERSION 1.00

'''
from bisect import bisect_left


class BilinearInterpolation(object):
    ''' Bilinear interpolation. 
        Usage:

            table = BilinearInterpolation(
            x_index=(54.458333, 54.5), 
            y_index=(17.041667, 17.083333), 
            values=((31.945, 31.866), (31.993, 31.911)) )

            print(table(54.4786674627, 17.0470721369))
    '''



    def __init__(self, x_index, y_index, values):
        self.x_index = x_index
        self.y_index = y_index
        self.values = values

    def __call__(self, x, y):
        # local lookups
        x_index, y_index, values = self.x_index, self.y_index, self.values

        i = bisect_left(x_index, x) - 1
        j = bisect_left(y_index, y) - 1

        x1, x2 = x_index[i:i + 2]
        y1, y2 = y_index[j:j + 2]
        z11, z12 = values[j][i:i + 2]
        z21, z22 = values[j + 1][i:i + 2]

        return (z11 * (x2 - x) * (y2 - y) +
                z21 * (x - x1) * (y2 - y) +
                z12 * (x2 - x) * (y - y1) +
                z22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1))
