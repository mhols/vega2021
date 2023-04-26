#compatible python 3, il a ete converti avec 2to3 et modification pickle
#ce code doit comparer deux spectres:
# 1) ThAr de reference UVES, calibre en longueur d'onde au mieux
#       voir http://astronomy.swin.edu.au/~mmurphy/thar/index.html
#  ici comparaison non pas au spectre même, mais à la liste de raies correspondante!
# 2) ThAr reduit par le DRS de Neo-Narval
#le but est de quantifier la precision du polynome de dispersion du DRS, ordre par ordre,
#region par region, le mieux serait raie par raie
#il faut lire les deux spectres, puis faire une fenetre glissante qu'on crosscorrelle avec la ref

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

yyy = [  807.74891654,   919.95541645,  1103.22732311,  1110.02040113,
        1352.5883779 ,  1441.04816395,  1645.55264777,  1932.6389426 ,
        2216.95243925,  2826.81228378,  3560.38575644,  9177.05836349,
       29495.92857131, 45006.63382431, 52219.44439841, 40922.512117  ,
       30367.71916259, 12240.50001547,  4005.2079407 ,  3140.22603902,
        2608.12668338,  2266.77052957,  1990.28250983,  1721.8064991 ,
        1478.99667981,  1365.73615686,  1118.18180981,  1100.31646158,
        1068.1224808 ,  1036.92762198,   937.73875075,   814.36162358,
         785.0888952 ]

xxx = np.arange(len(yyy))

def gauss(x, A, mu, sigma, y_offset):
# A, sigma, mu, y_offset = p
   return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y_offset



iymin=np.argmin(yyy)
iymax=np.argmax(yyy)
ixmin=np.argmin(xxx)
ixmax=np.argmax(xxx)
rA=(yyy[iymax]-yyy[iymin])
rsigma=0.2*(xxx[ixmax]-xxx[ixmin])
rmu=0.
ry_offset=yyy[iymin]
#        print("ry_offset: ",y_offset)
#            gy=gauss(reflambda[indrefextract],A,mu,sigma,y_offset)
rinit_vals = [rA, rmu,rsigma,ry_offset]  # for [amp, cen, wid]
#        print ("rinit_vals fit EWref:", rinit_vals)

rgopt, rcovar = curve_fit(gauss,xxx-mean(xxx),yyy,p0=rinit_vals)


