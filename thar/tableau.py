#c

import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import pickle
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import math
import string
import sys
import os
import json
import scipy
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from scipy.special import erf
from numpy import exp
from scipy.signal import savgol_filter
from scipy.signal import correlate
from scipy.optimize import curve_fit
#

n=5
o=7
crossmax=6

res=[]


count = 0
for i in range (n+1):
    for j in range (o+1):
        count+=1
        if (i+j <= crossmax) | (i==0) | (j==0):
            print(count,i,j)
            res.append((i,j))
        
print(res)
    

