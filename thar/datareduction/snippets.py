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
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import json
import scipy
from scipy.optimize import curve_fit
import extract
#

#NROWS=4208
#ORDERS=range(21,58)
CLUM=3e5

REF_SPECTRUM = '../reffiles/thar_spec_MM201006.dat'
REF_ATLASLINES = '../reffiles/thar_UVES_MM090311.dat'
EXCLUSION = '../reffiles/excluded.dat'
IVANLINES = '../reffiles/ivan.txt'
#seuil en ADU
SEUIL = 2000.
SEUILR = 800.
#vrange in km/s
VRANGE= 9.
    
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



def snippets(extractor,nvoie,order):
   

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
    



    #lecture du fichier ascii catalogue de raies thar_UVES_MM090311.dat
    #33322.4046  3000.109256  0.450 Ar II  W
    #33319.9704  3000.328442 -0.007 XX 0   P
    #atlasname = 'thar_UVES_MM090311.dat'
    with open(REF_ATLASLINES, 'r') as f:
        alines = f.readlines()

    num_alines=len(alines)
    npt=num_alines
    atlasline = np.zeros (npt, dtype =  float)
    for l in range (0,npt):
        linear = alines[l]
        sfr = linear.split()
        atlasline[l] = float(sfr[1])
    
    
    with open(IVANLINES, 'r') as kk:
        ilines = kk.readlines()
    num_ilines=len(ilines)
    npti=num_ilines
    iline = np.zeros (npti, dtype =  float)
    for l in range (0,npti):
        iline[l] = float(ilines[l])

     
    o = order
    if nvoie == 1:
        lam,flux = myext.get_lambda_intens1(o)
    elif nvoie == 2:
        lam,flux = myext.get_lambda_intens2(o)
    elif nvoie == 3:
        lam,flux = myext.get_lambda_intens3(o)
    else:
        raise Exception('no such voie')
        
    
    res = []

    
    
    #####
    
    #selectionner les raies de reference dans l'intervalle spectral
    #lam[-1] deniere valeur de lambda
    
    indexx=np.where((lam[0] <atlasline) & (atlasline < lam[-1]))
    atlasext=atlasline[indexx]
    
    indexi=np.where((lam[0] <iline) & (iline < lam[-1]))
    ivanext=iline[indexi]
   
       
    # ICI ON DOIT LIRE exlude.dat et filter les atlasline d'office
    # cad ne pas prendre les raies tombant dans des endroits interdits
    exclusion = np.loadtxt(EXCLUSION)
    
    #I = []
    I, = np.where(exclusion[:,0] == o)
    #print(atlasext)
    
    goodlines = []
    for l in atlasext:
        for i in I:
            if l >= exclusion[i,1] and l <= exclusion[i,2]:
                break
        goodlines.append(l)
        
    atlasext=np.array(goodlines)
    #(atlasext)
    
    plt.figure(figsize=(16,6))
    plt.plot(lam,flux)
    plt.plot(refwave,refintens)
    for ll in atlasext:
        plt.vlines(ll,0.,20000.,'y')
    for kk in ivanext:
        plt.vlines(kk,10000.,20000.,'r')
    
    # on ne veut choisir que les raies de l'atlas qui se retrouvent dans la zone atlasext[k] +/- vrange, et qui ont un flux max au dessus du seuil. Ca reduit la liste.
        
    latlasext=atlasext*(1.-VRANGE/CLUM)
    ratlasext=atlasext*(1.+VRANGE/CLUM)
        
    numlines=atlasext.shape[0]
    maxi = np.zeros (numlines, dtype =  float)
    maxir = np.zeros (numlines, dtype =  float)
        
    #ici on selectionne des snippets autour de chaque raie du catalogue
    #de reference. Il faut que la grille en longueur d'onde soit des le
    #depart suffisamment bonne pour que la raie observee tombe dans le snippet
    
   
    snip = []
    for l,c,r in zip(latlasext,atlasext,ratlasext):
        indext =  np.where((lam > l) & (lam < r))
        wave=lam[indext]
        inte=flux[indext]
        
        indextr = np.where((refwave > l) & (refwave < r))
        waver = refwave[indextr]
        inter = refintens[indextr]
        
  
       
      
    # selectionner que les raies du ThAr observes au dessus du seuil.
    # pour chaque raie k on determine le maximum de flux
        distmax = 2.
        goodsnippet = True
        goodsnippet = goodsnippet and (np.max(inte) - np.min(inte)) >= SEUIL
        goodsnippet = goodsnippet and (np.max(inter) - np.min(inter)) >= SEUILR
        goodsnippet = goodsnippet and (np.argmax(inte)>= distmax) and (np.argmax(inte) <= inte.shape[0]-distmax)
        if goodsnippet:
            snip.append({"o":o,"refwave":c ,"wave":wave,"inte":inte})
            
            plt.vlines(c,0.,10000.,'b')
            plt.plot(wave,inte,"r")
            
    plt.show()
    
    return snip
  


if __name__ == "__main__":
    print ('es geht los')
    myext = extract.Extractor(DATADIR='./../datafiles',VOIE_METHOD='SUM_DIVIDE_CENTRALROW')
    myext.set_fitsfile('../datafiles/NEO_20220903_191404_th0.fits')

    snip =  snippets(myext,1, 44)
    print ('snip', snip, len(snip))



"""
    #plotten im ipython
o=54
lam,flux = myext.get_lambda_intens1(o)
plt.plot(lam,flux,"b")
plt.plot(refwave,refintens,"r")
plt.show()
    
    
for o in extract.ORDERS:
   lam,flux = myext.get_lambda_intens1(o)
   plt.figure(figsize=(16,6))
   plt.xlim(lam[0],lam[-1])
   plt.ylim(0.,100000)
   plt.plot(lam,flux,"b")
   
   
   plt.plot(refwave,refintens,"r")
   plt.show()

how to read the snippets:

order = np.array(snip[ll]['o'])
lamc = np.zeros(len(snip))
pixc = np.zeros(len(snip))
for ll in range(len(snip)):
    lamc[ll]=np.array(snip[ll]['refwave'])
    pixc[ll]=np.array(np.median(snip[ll]['pixel']))
plt.figure(figsize=(16,6))
plt.plot(order*lamc,pixc)
plt.show()



plt.figure(figsize=(16,6))
for order in extract.ORDERS:
    snip =  snippets(myext,1, order)
    lamc = np.zeros(len(snip))
    pixc = np.zeros(len(snip))
    for ll in range(len(snip)):
        order = np.array(snip[ll]['o'])
        lamc[ll]=np.array(snip[ll]['refwave'])
        pixc[ll]=np.array(np.median(snip[ll]['pixel']))
    plt.plot(order*lamc,pixc,'o')
    plt.plot(order*lamc,pixc)
plt.show()
"""
