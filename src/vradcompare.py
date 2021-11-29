'''
Created on Jul 9, 2014

@author: hols
'''

from spectralutil import *
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectral

PROGDIR = os.path.dirname(__file__)
DATADIR = os.path.join(PROGDIR, '../data')

DATAFILE = os.path.join(DATADIR, 'filematrix_nn.dat')
time, velocity, intensity, list_time, list_inte = load_data(DATAFILE)

vrad_c = vrad_corr(velocity, intensity)
vrad_m = vrad( 0.3, velocity, intensity)
vrad_t = vrad_translat( intensity )

vspan = vspan((0.5,0.7), (0.15,0.3), time, velocity, intensity)
vrad_c -= vrad_c.mean()
vrad_m -= vrad_m.mean()
vrad_t -= vrad_t.mean()


#plt.plot( db[:,0], vsb )
#plt.plot( tim, vsh )



plt.figure()
plt.plot(vrad_c, vrad_m, 'o')

plt.figure()
plt.plot(vrad_c, vrad_t, 'o')

plt.figure()
plt.plot(vrad_m, vrad_t, 'o')

plt.figure()
plt.plot(vrad_c, vspan, 'ro')

plt.figure()
plt.plot(vrad_m, vspan, 'ro')

plt.figure()
plt.plot(vrad_t, vspan, 'ro')



plt.show()