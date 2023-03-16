from numpy import *
from scipy.interpolate import interp1d

from settings import kwargs as kwargsv
from extract import *

myextv = Extractor(**kwargsv)
myextv.set_fitsfile('/home/hols/vega2021/thar/datafiles/NEO_20220903_191404_th0.fits')

from settingsmoon import kwargs as kwargsm
kwargsm['OFFSET_LAMBDA']=0
myextm = Extractor(**kwargsm)
myextm.set_fitsfile('../lune_raw/NEO_20200202_173811_th0.fits')

def hom(d, l, x):
    return x + l*(x-x.shape[0]/2) + d

def optimize(d, l):

    def fitfunction(d, l, x, m, vv):
        xx = hom(d, l, x)
        try:
            return sum( m * vv(xx) )
        except:
            return 0


    d, l = meshgrid(d,l)
    d = d.ravel()
    l = l.ravel()

    res = {}
    for o in kwargsm['ORDERS']:
        v = myextv.voie1[o]
        v = where(isnan(v), 0, v)
        m = myextm.voie1[o]
        m = where(isnan(m), 0, m)
        n = arange(v.shape[0])
        vv = interp1d(n, v)
        N = 500
        x = n[N:-N]
        m = m[N:-N]

        dist=[fitfunction(dd, ll, x, m, vv) for dd, ll in zip(d,l)]
        i = argmax(dist)

        res[o] = (d[i], l[i])

        print('on ', o, res[o])
    return res
