#compatible python 3, il a ete converti avec 2to3 et modification pickle
#ce code doit comparer deux spectres:
#lit le fichier FP de Neo-Narval.
# la sortie est un fichier .json donnant pour chaque ordre le flux par pixel et la
# longueur d'onde preliminaire par pixel, et ceci pour la voie 1, 2, et 3!

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
#
name1 = "FP_voie_1.dat"
name2 = "FP_voie_2.dat"
name3 = "FP_voie_3.dat"


wave2Dname = 'wave2D_Matthias.dat'
wave2D = np.loadtxt(wave2Dname)

#lecture du fichier FP reduit avec le DRS de Neo-Narval

a=pyfits.open('NEO_20220517_201653_mo2.fits')

# *****************************************
# hier nimmt man Matthias polynom


# *****************************************
wave1=a[1].data['Wavelength1']
wave1=wave2D
intens1=a[1].data['Beam1']

wave2=a[1].data['Wavelength2']
wave2=wave2D
intens2=a[1].data['Beam2']

wave3=a[1].data['Wavelength3']
wave3=wave2D
intens3=a[1].data['Beam3']


orderlim=a[2].data['Orderlimit']

ordernum = zeros (orderlim.size, dtype = float)
centrallam = zeros (orderlim.size, dtype = float)
tordernum = zeros (orderlim.size, dtype = float)

res1=[]
res2=[]
res3=[]


#for j in range (20,21):
#for j in range (31,orderlim.size-2):
for j in range (2,orderlim.size-2):

    ordernum[j] = j
# vrai numero d'ordre
    tordernum[j] = ordernum[j]+21
    print("j: ",j, " ordernum[j]", ordernum[j], " true ordernum[j]: ",tordernum[j])

    wext1=wave1[orderlim[j]:orderlim[j+1]]
    iext1=intens1[orderlim[j]:orderlim[j+1]]
    wext2=wave2[orderlim[j]:orderlim[j+1]]
    iext2=intens2[orderlim[j]:orderlim[j+1]]
    wext3=wave3[orderlim[j]:orderlim[j+1]]
    iext3=intens3[orderlim[j]:orderlim[j+1]]
  
#    centrallam[j]=(wave1[orderlim[j+1]]+wave1[orderlim[j]])/2.
    centrallam[j]=(np.max(wext3)-np.min(wext3))/2.+np.min(wext3)
#    plt.plot(wext3,iext3)
#    plt.show()
    """
    toto=np.logical_and(iext3[1:-1]-iext3[:-2]>0., iext3[2:]-iext3[1:-1]<0)
    I,=np.where(toto==True)
    #maximum jedes FP peaks
    Msnip=I+1
    lsnip=np.zeros((I.shape[0],9))
    for i in range(I.shape[0]):
        lsnip[i,:]=iext3[Msnip[i]-4:Msnip[i]+5]
    fitting=np.polyfit(np.arange(-4,5),lsnip.T,2)
#    fitmax=fitting[]
    """
    
    """
#    ym1=iext3[:-2]
#    y0=iext3[1:-1]
#    yp1=iext3[2:]
# wir versuchen pm 2 pix:
    ym1=iext3[:-4]
    y0=iext3[2:-2]
    yp1=iext3[4:]
    
#    deltamax=(ym1-yp1)/(2*(ym1-2*y0+yp1))
    deltamax=(ym1-yp1)/(ym1-2*y0+yp1)
#    deltamax[I]
    truemax=I+2+deltamax[I]
    plt.plot(truemax[1:]-truemax[:-1],'o')
    plt.show()
    plt.plot(iext3,'o')
    plt.vlines(truemax,0,3000,'r')
    plt.show()
    """

    res1.append({'true_order_number': tordernum[j],
        'wave1': wext1.tolist(),
        'flux1': iext1.tolist()
        })
        
    res2.append({'true_order_number': tordernum[j],
        'wave2': wext2.tolist(),
        'flux2': iext2.tolist()
        })
    res3.append({'true_order_number': tordernum[j],
        'wave3': wext3.tolist(),
        'flux3': iext3.tolist()
        })

with open(name1+".json", 'w') as outfile:
        json.dump(res1, outfile, indent=2)
with open(name2+".json", 'w') as outfile:
        json.dump(res2, outfile, indent=2)
with open(name3+".json", 'w') as outfile:
        json.dump(res3, outfile, indent=2)


"""
toto=np.logical_and(iext3[1:-1]-iext3[:-2]>0., iext3[2:]-iext3[1:-1]<0)
I,=np.where(toto==True)
ym1=iext3[:-2]
y0=iext3[1:-1]
yp1=iext3[2:]
deltamax=(ym1-yp1)/(2*(ym1-2*y0+yp1))
deltamax[I]
truemax=I+1+deltamax[I+1]
plt.plot(truemax[1:]-truemax[:-1],'o')

plt.vlines(I+1,0,3000)
plt.plot(iext3,'o')
"""
