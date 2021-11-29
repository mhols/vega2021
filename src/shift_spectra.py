'''
Created on Jul 10, 2014

@author: hols
'''

from spectralutil import *
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectral
from scipy import interpolate

PROGDIR = os.path.dirname(__file__)
DATADIR = os.path.join(PROGDIR, '../data')

DATAFILE = os.path.join(DATADIR, 'filematrix.dat')

nval      = 201                # number of data points

coltime = 0 # colum of time values
colspec = 1
colval  = colspec + nval
colvul  = colval + nval

data = np.loadtxt(DATAFILE)

print data.shape
time = data[:,coltime:colspec].ravel() - 2.45e6

velocity = data[0, colspec:colval]      # velocities of bins
intens   = data[:, colval:colvul]       # intensities

n,d = intens.shape
v=np.arange(d)
tck = interpolate.splrep(v, intens[0,:], s=0)


freq = 1./0.76

for i  in xrange(n):
    row = data[i,:]
    delta = 0.005*np.sin( row[0] * 2 * np.pi * freq )
    data[i,colval:colvul] = interpolate.splev(v+delta, tck, der=0)

print data.shape

np.savetxt('filematrix_shifted.dat', data)
