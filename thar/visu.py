
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
#
clum=3e5




#lecture du fichier ThAr reduit avec le DRS de Neo-Narval
#a=pyfits.open('NEO_20200202_173811_th2.fits')
#a=pyfits.open('NEO_20200202_173811_thB.fits')

#urversion
a=pyfits.open('NEO_20220219_173134_mo2.fits')
wave1=a[1].data['Wavelength1']
intens1=a[1].data['Beam1']
wave2=a[1].data['Wavelength2']
intens2=a[1].data['Beam2']
wave3=a[1].data['Wavelength3']
intens3=a[1].data['Beam3']

b=pyfits.open('NEO_20220219_173048_th2.fits')
wave1th=b[1].data['Wavelength1']
intens1th=b[1].data['Beam1']
wave2th=b[1].data['Wavelength2']
intens2th=b[1].data['Beam2']
wave3th=b[1].data['Wavelength3']
intens3th=b[1].data['Beam3']


plt.figure(figsize=(16,6))
plt.plot(wave1th,intens1th,"y")
plt.plot(wave2th,intens2th,"r")
plt.plot(wave3th,intens3th,"b")

plt.figure(figsize=(16,6))
plt.plot(wave1,intens1,"y")
plt.plot(wave2,intens2,"r")
plt.plot(wave3,intens3,"b")

plt.show()
