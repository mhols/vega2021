#compatible python 3, il a ete converti avec 2to3 et modification pickle
#ce code doit comparer deux spectres:
#lit le fichier FP de Neo-Narval.
# la sortie est un fichier .json donnant pour chaque ordre le flux par pixel et la
# longueur d'onde preliminaire par pixel, et ceci pour la voie 1, 2, et 3!

import astropy.io.fits as pyfits
import numpy as np
from pylab import *
import sys
import os
import json
import scipy

def fp_write(fitsfile,wave2Dname,outfile1,outfile2,outfile3):
    """
    fitsfile: fichier input Neo-Narval FP de reference
    wave2Dname: comme arturo.dat, output lambda d'un calcul thar2D
    outfile1 (2,3): output json file contenant snippets FP
    """
    wave2D = np.loadtxt(wave2Dname)

    a=pyfits.open(fitsfile)

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

    with open(outfile1, 'w') as outfile:
        json.dump(res1, outfile, indent=2)
    with open(outfile2, 'w') as outfile:
        json.dump(res2, outfile, indent=2)
    with open(outfile3, 'w') as outfile:
        json.dump(res3, outfile, indent=2)


if __name__ == "__main__":
    fitsfile = "./datafiles/NEO_20220517_201653_mo2.fits"
    wave2dname = "./reffiles/wave2D_Matthias.dat"
    outfile1 = "NEO_20220517_201653_mo2_voie1.json"
    outfile2 = "NEO_20220517_201653_mo2_voie2.json"
    outfile3 = "NEO_20220517_201653_mo2_voie3.json"
    
    fp_write(fitsfile,wave2dname,outfile1,outfile2,outfile3)


