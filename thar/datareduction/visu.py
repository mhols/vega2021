
#verifie flux sur 3eme voie
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import pickle
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import math
import string
import sys
import scipy
import os
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from scipy import exp
from scipy.signal import savgol_filter
from scipy.signal import correlate
from scipy.optimize import curve_fit
import extract
from settings import kwargs as refkwargs
import regal

#
clum=3e5

def get_ext(f_thar):
    try:
        myext = store.get(f_thar)
        print('retrieving precomputed object for ',  f_thar)
    except:
        myext = extract.Extractor(**kwargs)
        myext.set_fitsfile(f_thar)
        store.store(f_thar, myext)
    return myext
    
store = regal.Store()
myext = get_ext('/Users/boehm/Desktop/vega2021/thar/51Peg_raw/2020_0912/NEO_20200912_190502_th0.fits')

try:
    settingsmodule = os.environ['SETTINGSMODULE']
except:
    raise Exception('wrong SETTINGSMODULE')

try:
    exec('from %s import *'%(settingsmodule,))
except:
    raise Exception('could not import {settingsmodule}')






#lecture du fichier ThAr reduit avec le DRS de Neo-Narval
#a=pyfits.open('NEO_20200202_173811_th2.fits')
#a=pyfits.open('NEO_20200202_173811_thB.fits')

#urversion
#ThAr moon
a=pyfits.open('../51Peg_raw/2020_0912/HOBO_NEO_20200912_190502_th1.fits')

#ThAr Referenznacht
#a=pyfits.open('../reffiles/HOBO_NEO_20220903_191404_th1.fits')

#moon
#a=pyfits.open('../lune_raw/HOBO_NEO_20200202_182715_st1.fits')

wave1=a[1].data['wavelength_1']
intens1=a[1].data['flux_1']
wave2=a[1].data['wavelength_2']
intens2=a[1].data['flux_2']
ordernumber=a[1].data['true_order_number']

"""
np.savetxt('hobolambda1improved.dat',np.array(wave1),fmt=['%.7f'])
np.savetxt('hobolambda2improved.dat',np.array(wave2),fmt=['%.7f'])
"""

b=pyfits.open('../lune_res/HOBO_NEO_20200202_173811_th1.fits')
wave1b=b[1].data['wavelength_1']
intens1b=b[1].data['flux_1']
wave2b=b[1].data['wavelength_2']
intens2b=b[1].data['flux_2']
#ordernumber=a[1].data['true_order_number']


"""
b=pyfits.open('../lune_res/HOBO_NEO_20200202_173811_th1.fits')
wave1=a[1].data['wavelength_1']
intens1=a[1].data['flux_1']
wave2=a[1].data['wavelength_2']
intens2=a[1].data['flux_2']
ordernumber=a[1].data['true_order_number']
"""

exclusion = np.loadtxt(EXCLUSION)
    

"""
plt.figure(figsize=(16,6))
plt.plot(wave1,intens1,"b")
plt.plot(wave2,intens2,"r")
"""
plt.figure(figsize=(16,6))


#ORDERS=(33,33)
for o in myext.ORDERS:

    I=np.where(ordernumber == o)
    GI = myext.beams[o].I
    
    II, = np.where(exclusion[:,0] == o)
    
    wav1=wave1[I]
    flu1=intens1[I]
    
    wav1GI=wav1[GI]
    flu1GI=flu1[GI]
    
    for i in II:
            indexx, = np.where((wav1 >= exclusion[i,1]) & (wav1 <= exclusion[i,2]))
#            print(wav1[indexx])
            xe=[exclusion[i,1],exclusion[i,2]]
            ye=[-10.,-10.]
            if (o % 2) == 0:
                plt.plot(xe,ye,"r")
            else:
                plt.plot(xe,ye,"b")
#           flu1[indexx]=0.
    
    plt.ylim(-20,1000)
    if (o % 2) == 0:
#        plt.plot(wav1,flu1,'r')
        plt.plot(wav1GI,flu1GI,'r')
    else:
#        plt.plot(wav1,flu1,'b')
        plt.plot(wav1GI,flu1GI,'b')
#    print(o,min(wav1),(min(wav1)+max(wav1))/2, max(wav1))

plt.show()

store.save()
