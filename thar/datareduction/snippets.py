#ce code genere les snippets utilises pour la calibration 2D
#ce code doit comparer deux spectres:
# 1) ThAr de reference UVES, calibre en longueur d'onde au mieux
#       voir http://astronomy.swin.edu.au/~mmurphy/thar/index.html
#  ici comparaison non pas au spectre même, mais à la liste de raies correspondante!
# 2) ThAr reduit par hobo sur Neo-Narval
# au depart on s'est servi de artlambda.dat (genere par arthols.py) de la même nuit,
# afin d'avoir une grille en lambda approximative
#


import astropy.io.fits as pyfits
import numpy as np
import sys
import os
import json
import scipy
from scipy.optimize import curve_fit
from extract import *
#

def print(*x):
    pass

#NROWS=4208
#ORDERS=range(21,58)
#clum=3e5

REF_SPECTRUM = './reffiles/thar_spec_MM201006.dat'
REF_ATLASLINES = './reffiles/thar_UVES_MM090311.dat'
#seuil en ADU
SEUIL = 2000.
#vrange in km/s
VRANGE= 9.
def snippets(nvoie,order):
   

    """
    #lecture du fichier ascii du spectre ThAr de reference (UVES)
    #ce spectre contient une premiere colonne lambda, une seconde intensite,...
    #213358 lignes

    #refname = 'thar_spec_MM201006.dat'
    #thar
    #213358 4
    #3033.0867480159  103.99224526863826    0.36633407377209    0.36633407377209     97.132933537024
    """
    #nombre de lignes du catalogue de raies
    
    num_lines=len(open(REF_SPECTRUM,'r').readlines())
    npt=num_lines-2
    linecount = 0
    refwave = np.zeros (npt, dtype =  float)
    refintens = np.zeros (npt, dtype =  float)


    fr = open(REF_SPECTRUM, 'r')
    linefr = fr.readline()
    linefr = fr.readline()
    for l in range (0,npt):
            linefr = fr.readline()
            sfr = linefr.split()
            refwave[l] = float(sfr[0])
            refintens[l] = float(sfr[1])
    fr.close


    #lecture du fichier ascii catalogue de raies thar_UVES_MM090311.dat
    #33322.4046  3000.109256  0.450 Ar II  W
    #33319.9704  3000.328442 -0.007 XX 0   P
    #atlasname = 'thar_UVES_MM090311.dat'
    num_alines=len(open(REF_ATLASLINES,'r').readlines())
    npt=num_alines
    atlasline = np.zeros (npt, dtype =  float)
    ar = open(atlasname, 'r')
    for l in range (0,npt):
            linear = ar.readline()
            sfr = linear.split()
            atlasline[l] = float(sfr[1])
    fr.close
     
    o = order
    voie=voie+'nvoie'
    myext = restart()
    flux = myext.voie
    lam  = myext.get_lambda(o)

    res = []
    o = order
    
    # ICI ON DOIT LIRE exlude.dat et filter les atlasline d'office
    # cad ne pas prendre les raies tombant dans des endroits interdits
    
    
    
    
    
    
    
    #####
    
    #selectionner les raies de reference dans l'intervalle spectral
    indexx=np.where(lam[0] <atlasline) & (atlasline < lam[lam.size-1])
    atlasext=atlasline[indexx]
    

    # on ne veut choisir que les raies de l'atlas qui se retrouvent dans la zone atlasext[k] +/- vrange, et qui ont un flux max au dessus du seuil. Ca reduit la liste. En plus, il faut le faire sur les deux voies et appliquer le mini.
        
    latlasext=atlasext*(1.-VRANGE/CLUM)
    ratlasext=atlasext*(1.+VRANGE/CLUM)
        
    numlines=int(np.array(indexx).size)
    maxi = np.zeros (numlines, dtype =  float)
    maxir = np.zeros (numlines, dtype =  float)
        
    #ici on selectionne des snippets autour de chaque raie du catalogue
    #de reference. Il faut que la grille en longueur d'onde soit des le
    #depart suffisamment bonne pour que la raie observee tombe dans le snippet
    
    for k in range (numlines):
    # snippets autour des raies de references, extraits dans spectre obs
        indext =  np.where((lam <ratlasext[k]) & (lam >latlasext[k]))
        wav=lam[indext]
        int=flux[indext]
        
    # snippets autour des raies de references, extraits dans spectre de ref
        indextr =  np.where((refwave<ratlasext[k]) & (refwave>latlasext[k]))
        wavr=refwave[indextr]
        intr=refintens[indextr]
        
        
    # selectionner que les raies du ThAr observes au dessus du seuil.
    # pour chaque raie k on determine le maximum de flux

        maxi[k] = max(int)
        maxir[k] = max(intr)
     
    # ici on filtre sur les maxima
    indseuil= np.where((maxi > SEUIL) & (maxir > SEUIL))
        
    #number of selected reference lines
    nslines=int(np.array(indseuil).size)
    # liste des raies au dessus du seuil, avec bords rouges et bleus
    satlasext=atlasext[indseuil]
    slatlasext=latlasext[indseuil]
    sratlasext=ratlasext[indseuil]
        
    
    for k in range (nslines):
            linecount += 1
            indextr = np.where((wavr <sratlasext[k]) & (wavr>slatlasext[k]))
            wr=wavr[indextr]
            ir=intr[indextr]
        
            indext =  np.where((wav<max(wr)) & (wav>min(wr)))
            wavext=wav[indext]
            intext=int[indext]
            
            print(k,satlasext(k),wavext,intext)
   
   
"""
                res1.append({
                    'ThAr_line_number': linecount,
                    'ref_lambda': satlasext[k],
                    'true_order_number': tordernum[j],
                    'pixels_extract': pixel1.tolist(),
                    'flux_values_extract': yyy.tolist()
                })

"""

if __name__ == "__main__":
    print ('es geht los')
    myext = restart()
    myext.set_fitsfile('./datafiles/NEO_20220903_191404_th0.fits')

    snippets =  snippets(myext, 44)
    print ('snippets', snippets)
    
