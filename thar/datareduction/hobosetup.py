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
from settings import *

#ORDERS = range(21, 60)
#NROWS = 7208

#pointage de deux raies du meme ordre, bien eloigne, une fois sur la reference (nuit  Vega) puis dans
#la nuit a reduire
#on a fait tourner visu_duo_pixel.py pour plotter les deux spectres et mesurer en pixel position deux
#raies extremes, afin de mesurer decalage et dilatation eventuelle
#dans un premier temps on considere que voie1 et  voie2 se comportent de la meme maniere differentielle
hoboref1 = 219.1
hoboref2 = 3984.7
hobo1 = 234.5
hobo2 = 4004.7

# longueur d'onde raie 1 tombe sur hoboref1, la meme tombe sur hobo1 dans le nouveau spectre
# idem pour la longueur d'onde 2. On a donc l'egalite: lambda1(m * hobo1 + b)  = lambda1(hoboref1)
#lambda1(m * hobo2 + b)  = lambda1(hoboref2)
# donc (hoboref2-hoboref1)/(hobo2-hobo1) = m
#


m = (hoboref2-hoboref1)/(hobo2-hobo1)
b = hoboref2 - m*hobo2



myext=Extractor(VOIE_METHOD='SUM_DIVIDE_CENTRALROW')

# fichier de reference, "nuit Vega", a partir de artlambdacorrect...dat

HOBOREFFILE = "/Users/boehm/Desktop/vega2021/thar/reffiles/HOBO_NEO_20220903_191404_th1.fits"

myext.set_fitsfile(HOBOREFFILE)

wave = []
waver = []
waveref=np.zeros(NROWS, dtype = float)
pixelref=np.arange(NROWS)
pixel=np.arange(NROWS)

for o in ORDERS:
    waveref = get_lambda(o)
    
# in pixeln:
# Wellenlaenge an pixel (m * hoboref1 + b) = Wellenlaenge an pixel hobo1 im neuen file
# m * hoboref2 + b = hobo2

#position fractionnaire de pixel sur la grille pixelref (donc waveref)

    pixelref = m * pixel +b

    p1 = np.poly1d(np.polyfit(pixelref, waveref, 2))


#plt.plot(wa,va/va[1000:3000].max())
#plt.plot(p1(xh_a),vh1/vh1[1000:3000].max())


#plt.show()

#ici on est dans la grille hols, cad 4208 pix, avec 0-100 croix
    wh=p1(pixelref)
#    plt.plot(xh,wh)

#ici on est dans la grille art, cad 7800 pix, excluant la croix
#plt.plot(xa,wa)

#plt.show()
    waver.append(waveref)
    wave.append(wh)
    
tmp=[]
tmpr=[]
for o in ORDERS:
    t=o-ORDERS[0]
    for i in range(NROWS):
        tmp.append(wave[t][i])
        tmpr.append(waver[t][i])
        

np.savetxt('../reffiles/hobo.txt',np.array(tmp),fmt=['%.7f'])

