import matplotlib.pyplot as plt
import numpy as np
import pickle
import astropy.io.fits as pyfits
import os
import sys
from scipy.ndimage import minimum_filter1d, generic_filter
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline, bisplrep, bisplev
from scipy.signal import convolve2d
from numpy.polynomial import Polynomial
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
import pandas as pd
from extract import *

ORDERS = range(21, 22)
#NROWS = 7208
REF_SPECTRUM = '../reffiles/thar_spec_MM201006.dat'

hol1 = 1405.39
hol2 = 3712.15
ref1 = 2462.76
ref2 = 6200.49


pixlam54 = [1405.30, 4158.59],[2351.04, 4198.32],[2409.47, 4200.67],[2666.87, 4210.92],[3712.19, 4250.31]
pixlam22 = [666.02,10133.56],[935.86,10163.378],[2192.09,10308.54],[2584.03,10332.74],[4106.9, 10470.1]
toto = [666.02,10133.56],[2584.03,10332.74],[4106.9, 10470.1]
tot = np.array(toto)
lam=np.zeros(tot.size, dtype = float)
pix=np.zeros(tot.size, dtype = float)
lam=tot[:,0]
pix=tot[:,1]
ppix = np.poly1d(np.polyfit(lam, pix, 2))
plt.plot(ppix(xx),vh1)
plt.plot(refwave,refintens)
plt.show()


"""

hol1 = 644.20
hol2 = 3769.86
art1 = 1088.85
art2 = 7340.84
"""

with open(REF_SPECTRUM, 'r') as f:
    lines = f.readlines()

num_lines=len(lines)
npt=num_lines-2
linecount = 0
refwave = np.zeros (npt, dtype =  float)
refintens = np.zeros (npt, dtype =  float)

for l in range (0,npt):
    linefr = lines[2+l]
    sfr = linefr.split()
    refwave[l] = float(sfr[0])
    refintens[l] = float(sfr[1])

n = (hol1-hol2)/(ref1-ref2)
b = hol1 - n*ref1

m = (ref1-ref2)/(hol1-hol2)
bs = ref1 - m*hol1


def save(myext):
    pickle.dump(myext, open('pickledump', 'wb'))

def reload():
    myext = pickle.load(open('pickledump', 'rb'))
    return myext

def restart(**kwargs):

    myext=Extractor('../datafiles/NEO_20220903_190703_fla.fits',
        DATADIR='../datafiles',**kwargs)
    return myext

artfile = "/Users/boehm/Desktop/vega2021/thar/datafiles/NEO_20220219_173048_th2.fits"
#artfile = "/Users/boehm/Desktop/vega2021/thar/datafiles/NEO_20230125_175423_th2.fits"
a=pyfits.open(artfile)

wave1=a[1].data['Wavelength1']
intens1=a[1].data['Beam1']
wave2=a[1].data['Wavelength2']
intens2=a[1].data['Beam2']
wave3=a[1].data['Wavelength3']
intens3=a[1].data['Beam3']

orderlim=a[2].data['Orderlimit']

#myext=restart()
myext = restart(VOIE_METHOD='SUM_DIVIDE_CENTRALROW')
myext.set_fitsfile('../datafiles/NEO_20220903_191404_th0.fits')


waveart = []
voie1hols = []
voie2hols = []

for o in ORDERS:
#   for i in range(NROWS):
#   waveart.append(o)
    selectedorder = o
    mult = selectedorder - 21
    va=intens1[mult*7800:(mult+1)*7800]
    wa=wave1[mult*7800:(mult+1)*7800]
    myext.beams[selectedorder]
    
    vh1 = myext.voie1[selectedorder]
    vh2 = myext.voie2[selectedorder]
    
    indexx=np.where((wa[0] < refwave) & (refwave < wa[-1]))
    waveext=refwave[indexx]
    intensext=refintens[indexx]
    """
    plt.figure(figsize=(16,6))
    plt.plot(wa,va)
    tit = str(o)
    plt.title(tit)
    plt.plot(waveext,intensext,"r")
    plt.show()
    """


###da = np.loadtxt('art.dat')
###va = da[:,1]

    xh=np.arange(vh1.size)



# n * art1 = hol1
# n * art2 = hol2


    xh_r= m*xh +bs
    """
    xa=np.arange(va.size)
    xa_original=xa
    xa=xa*n+b
    """
    

    xr=np.arange(intensext.size)
    xr_original=xr
    xr=xr*n+b
    

  


    p1 = np.poly1d(np.polyfit(xr_original, waveext, 5))


#plt.plot(wa,va/va[1000:3000].max())
#plt.plot(p1(xh_a),vh1/vh1[1000:3000].max())


#plt.show()

#ici on est dans la grille hols, cad 4208 pix, avec 0-100 croix
    wh=p1(xh_r)
#    plt.plot(xh,wh)

#ici on est dans la grille art, cad 7800 pix, excluant la croix
#plt.plot(xa,wa)

#plt.show()
    waveart.append(wh)
    voie1hols.append(vh1)
    voie2hols.append(vh2)
#np.savetxt('vspan.dat', np.column_stack((time, tmp)), fmt=['%.5f', '%.7f'])
#le fichier artlambda.dat contient la longueur d'onde de 39 ordres, un après l'autre,
#en commencant par l'ordre IR 21. Chaque ordre a des lambda croissants.
# ce fichier contient donc 39*4208 = 164112 lignes

tmp=[]
tmpvh1=[]
tmpvh2=[]
for o in ORDERS:
    t=o-ORDERS[0]
    for i in range(NROWS):
#        print(t,i,waveart[t][i])
        tmp.append(waveart[t][i])
        tmpvh1.append(voie1hols[t][i])
        tmpvh2.append(voie2hols[t][i])
        

np.savetxt('toto.dat',np.array(tmp),fmt=['%.7f'])
np.savetxt('voie1hols.dat',np.array(tmpvh1),fmt=['%.7f'])
np.savetxt('voie2hols.dat',np.array(tmpvh2),fmt=['%.7f'])

#lecture par
# waveart = np.loadtxt('artlambda.dat')

"""
whols=np.loadtxt('artlambda.dat')
vhols=np.loadtxt('voie1hols.dat')
order = 54
torder = order-ORDERS[0]
wa=wave1[torder*7800:(torder+1)*7800]
va=intens1[torder*7800:(torder+1)*7800]
wh=whols[torder*4208:(torder+1)*4208]
vh=vhols[torder*4208:(torder+1)*4208]
va=va/va[1000:3000].max()
vh=vh/vh[1000:3000].max()
plt.plot(wa,va)
plt.plot(wh,vh)
plt.show()
"""

"""
myext.set_fitsfile('/Users/boehm/Desktop/extract/fitsfiles/NEO_20220903_200017_st0.fits')
for o in ORDERS:
    plt.figure(figsize=(16,6))
    plt.plot(myext.voie1[o])
    plt.plot(myext.voie2[o],"r")
plt.show()



# fichier reduit par le DRS:
artfile = "/Users/boehm/Desktop/NEO_20220903_195942_st1.fits"
a=pyfits.open(artfile)

wave1=a[1].data['Wavelength1']
intens1=a[1].data['Beam1']
wave2=a[1].data['Wavelength2']
intens2=a[1].data['Beam2']
wave3=a[1].data['Wavelength3']
intens3=a[1].data['Beam3']

for o in ORDERS:
    selectedorder = o
    mult = selectedorder - 21
    va1=intens1[mult*7800:(mult+1)*7800]
    wa1=wave1[mult*7800:(mult+1)*7800]
    va2=intens2[mult*7800:(mult+1)*7800]
    wa2=wave2[mult*7800:(mult+1)*7800]
    plt.figure(figsize=(16,6))
    plt.plot(wa1,va1,"b")
    plt.plot(wa2,va2,"r")
plt.show()

for o in ORDERS:
    plt.figure(1, figsize=(15, 5))
    
    plt.subplot(211)
    plt.plot(myext.voie1[o])
    plt.plot(myext.voie2[o],"r")

    plt.subplot(212)
    selectedorder = o
    mult = selectedorder - 21
    va1=intens1[mult*7800:(mult+1)*7800]
    wa1=wave1[mult*7800:(mult+1)*7800]
    va2=intens2[mult*7800:(mult+1)*7800]
    wa2=wave2[mult*7800:(mult+1)*7800]
    plt.plot(wa1,va1,"b")
    plt.plot(wa2,va2,"r")
    plt.show()
"""
