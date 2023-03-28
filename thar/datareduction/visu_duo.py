
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
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from scipy import exp
from scipy.signal import savgol_filter
from scipy.signal import correlate
from scipy.optimize import curve_fit
from settings import *

#
clum=3e5






#lecture du fichier ThAr reduit avec le DRS de Neo-Narval
#a=pyfits.open('NEO_20200202_173811_th2.fits')
#a=pyfits.open('NEO_20200202_173811_thB.fits')

#urversion
#a=pyfits.open('../lune_res/HOBO_NEO_20200202_173811_th1.fits')
a=pyfits.open('../reffiles/HOBO_NEO_20220903_191404_th1.fits')
wave1=a[1].data['wavelength_1']
intens1=a[1].data['flux_1']
wave2=a[1].data['wavelength_2']
intens2=a[1].data['flux_2']
ordernumber=a[1].data['true_order_number']


b=pyfits.open('../51Peg_raw/2020_0912/HOBO_NEO_20200912_190502_th1.fits')
wave1b=b[1].data['wavelength_1']
intens1b=b[1].data['flux_1']
wave2b=b[1].data['wavelength_2']
intens2b=b[1].data['flux_2']
#ordernumber=a[1].data['true_order_number']
##wavehobo=np.loadtxt("../reffiles/hobo.txt")
#np.savetxt('../reffiles/hobo2.txt',np.array(wave1),fmt=['%.7f'])

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
plt.title("20220903: r/b, 20200202: y/g")


ORDERS=(33,33)
for o in ORDERS:

    I=np.where(ordernumber == o)
    
    II, = np.where(exclusion[:,0] == o)
    
    wav1=wave1[I]
    flu1=intens1[I]
    
    wav1b=wave1b[I]
#    waveh=wavehobo[I]

    flu1b=intens1b[I]
    
    
    for i in II:
            indexx, = np.where((wav1 >= exclusion[i,1]) & (wav1 <= exclusion[i,2]))
            print(wav1[indexx])
            xe=[exclusion[i,1],exclusion[i,2]]
            ye=[-10.,-10.]
            if (o % 2) == 0:
                plt.plot(xe,ye,"r")
            else:
                plt.plot(xe,ye,"b")
#           flu1[indexx]=0.
    
    plt.ylim(-20,1000)
    """
    if (o % 2) == 0:
        plt.plot(wav1,flu1,'r')
        plt.plot(wav1b,flu1b,'y')
    else:
        plt.plot(wav1,flu1,'b')
        plt.plot(wav1b,flu1b,'g')
    """

    plt.plot(flu1,'r')
    plt.plot(flu1b,'b')
  
#    print(o,min(wav1),(min(wav1)+max(wav1))/2, max(wav1))

plt.show()

