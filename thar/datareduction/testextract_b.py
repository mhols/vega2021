import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
import astropy.io.fits as pyfits
import os
import sys
#import utils
from scipy.optimize import curve_fit


from extract import *

def save(myext):
    pickle.dump(myext, open('pickledump', 'wb'))

def reload():
    myext = pickle.load(open('pickledump', 'rb'))
    return myext

def restart(**kwargs):

    myext=Extractor(**kwargs)

    return myext

#myext = pickle.load(open('pickledump', 'rb'))

"""
myext = restart()
i=51
b=myext.beams[i]
plt.figure(figsize=(9.5,9.5))
plt.imshow(myext.masterflat,vmin=0,vmax=200)
for i in ORDERS:
    b=myext.beams[i]
    plt.plot(b.y,b.x,'k.')
#    plt.plot(b._yy,b._xx,'k.')
plt.show()
"""

"""
myext = restart(VOIE_METHOD='SUM_DIVIDE_BLAZE')
"""
