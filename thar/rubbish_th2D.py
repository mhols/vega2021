#compatible python 3, il a ete converti avec 2to3 et modification pickle
#ce code doit comparer deux spectres:
# 1) ThAr de reference UVES, calibre en longueur d'onde au mieux
#       voir http://astronomy.swin.edu.au/~mmurphy/thar/index.html
#  ici comparaison non pas au spectre même, mais à la liste de raies correspondante!
# 2) ThAr reduit par le DRS de Neo-Narval
#le but est de quantifier la precision du polynome de dispersion du DRS, ordre par ordre,
#region par region, le mieux serait raie par raie
#il faut lire les deux spectres, puis faire une fenetre glissante qu'on crosscorrelle avec la ref
# la sortie est un fichier .json pour effectuer un polyfit 2D sur les raies

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
f=open('ThAr_Line_Cal_v2.pkl','rb')
b=pickle.load(f,encoding='bytes')
# liste de raies utilises par Arturo
thararturo = np.array(b[b'Lines Atlas'])
clum=3e5

#lecture du fichier ascii du spectre ThAr de reference (UVES)
#ce spectre contient une premiere colonne lambda, une seconde intensite,...
#213358 lignes

wave2Dname = 'wave2D_Matthias.dat'
wave2D = np.loadtxt(wave2Dname)

refname = 'thar_spec_MM201006.dat'
#thar
#213358 4
#3033.0867480159  103.99224526863826    0.36633407377209    0.36633407377209     97.132933537024

#nombre de lignes du catalogue de raies
num_lines=len(open(refname,'r').readlines())
npt=num_lines-2
linecount = 0
refwave = zeros (npt, dtype =  float)
refintens = zeros (npt, dtype =  float)


fr = open(refname, 'r')
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
atlasname = 'thar_UVES_MM090311.dat'
num_alines=len(open(atlasname,'r').readlines())
npt=num_alines
atlasline = zeros (npt, dtype =  float)
ar = open(atlasname, 'r')
for l in range (0,npt):
        linear = ar.readline()
        sfr = linear.split()
        atlasline[l] = float(sfr[1])
fr.close



#lecture du fichier ThAr reduit avec le DRS de Neo-Narval

a=pyfits.open('NEO_20220517_201608_th2.fits')

# *****************************************
# hier nimmt man Matthias polynom
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

def gauss(x, A, mu, sigma, y_offset):
# A, sigma, mu, y_offset = p
   return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y_offset

"""
def _gauss(x, A, mu, sigma, y_offset):
    return y_offset + 0.5 * np.sqrt(np.pi) * (erf((x+0.5-mu)/(np.sqrt(2)*sigma)) - erf((x-0.5-mu)/(np.sqrt(2)*sigma)))
"""

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
    
def plot_xcorr(x, y):
    "Plot cross-correlation (full) between two signals."
    N = max(len(x), len(y))
    n = min(len(x), len(y))

    if N == len(y):
        lags = np.arange(-N + 1, n)
    else:
        lags = np.arange(-n + 1, N)
    c = correlate(x / np.std(x), y / np.std(y), 'full')

    #plt.plot(lags, c / n)
    
    return lags,c,n


#outtableau1 = open("ThAr2D_voie1.dat","w")
#outtableau1 = open("ThAr2D_voie1.dat","w")
#outtableau2 = open("ThAr2D_voie2.dat","w")


#montrer = 'montremoi1'

montrer = 'non'

# orderlim va de 0 a 152100, la 41eme valeur est a 0. Il y a donc 39 ordres en tout
# in range 38,39 veut dire que tout ce qui est en dessous de 39, a partir de 0, donc 38
# uniquement


ordernum = zeros (orderlim.size, dtype = float)
centrallam = zeros (orderlim.size, dtype = float)
tordernum = zeros (orderlim.size, dtype = float)

res1=[]
res2=[]


for j in range (32,34):
#for j in range (31,orderlim.size-2):
#for j in range (2,orderlim.size-2):

    ordernum[j] = j
# vrai numero d'ordre
    tordernum[j] = ordernum[j]+21
    print("j: ",j, " ordernum[j]", ordernum[j], " true ordernum[j]: ",tordernum[j])
#ici on ne selectionne que les raies de l'atlas dans l'ordre j, index.size est ce nombre de raies du catalogue UVES
    indexx=where((wave1[orderlim[j]]<atlasline) & (atlasline < wave1[orderlim[j+1]-1]))
    atlasext=atlasline[indexx]
    wext1=wave1[orderlim[j]:orderlim[j+1]-1]
    iext1=intens1[orderlim[j]:orderlim[j+1]-1]
    wext2=wave2[orderlim[j]:orderlim[j+1]-1]
    iext2=intens2[orderlim[j]:orderlim[j+1]-1]
    wext3=wave3[orderlim[j]:orderlim[j+1]-1]
    iext3=intens3[orderlim[j]:orderlim[j+1]-1]

    
    
    indexxa=where((wave1[orderlim[j]]<thararturo) & (thararturo < wave1[orderlim[j+1]-1]))
    thararturoext=thararturo[indexxa]
    centrallam[j]=(wave1[orderlim[j+1]-1]+wave1[orderlim[j]])/2.
#
#  ICI on pourrait sortir les valeurs min et max de thararturo pour chaque ordre!
#
#orderlim[j+1]-orderlim[j] vaut 7800, il y a donc 7800 pixels (=Arturos) par ordre
    
    
# seuil en ADU pour qu'on puisse selectionner une raie de l'atlas. Il faut qu'elle depasse de ce seuil le zero. Ensuite, il faut donner une zone +/- autour de la raie du catalogue, pour detecter la valeur maxi du flux tu ThAr obs.
# le vrange est donc cette largeur en km/s
    seuil = 50.
    vrange = 9.
# on ne veut choisir que les raies de l'atlas qui se retrouvent dans la zone atlasext[k] +/- vrange, et qui ont un flux max au dessus du seuil. Ca reduit la liste. En plus, il faut le faire sur les deux voies et appliquer le mini.
    
    latlasext=atlasext*(1.-vrange/clum)
    ratlasext=atlasext*(1.+vrange/clum)
    
    numlines=int(np.array(indexx).size)
    numlinesa=int(np.array(indexxa).size)
    maxi = zeros (numlines, dtype =  float)
    maxir = zeros (numlines, dtype =  float)
    intline1 = zeros (numlines, dtype =  float)
    intline2 = zeros (numlines, dtype =  float)
    intline3 = zeros (numlines, dtype =  float)
    
# selectionner que les raies du ThAr observes au dessus du seuil.
    for k in range (numlines):
        indext1 =  where((wext1<ratlasext[k]) & (wext1>latlasext[k]))
        w1=wext1[indext1]
        i1=iext1[indext1]
    
        indext2 =  where((wext2<ratlasext[k]) & (wext2>latlasext[k]))
        w2=wext2[indext2]
        i2=iext2[indext2]
        
        indext3 =  where((wext3<ratlasext[k]) & (wext3>latlasext[k]))
        w3=wext3[indext3]
        i3=iext3[indext3]
        
        indextr =  where((refwave<ratlasext[k]) & (refwave>latlasext[k]))
        wr=refwave[indextr]
        ir=refintens[indextr]
  
        maxi[k]=min(max(i1),max(i2))
        maxir[k] = max(ir)
        
        
    indseuil= where((maxi > seuil) & (maxir > seuil))
    #number of selected reference lines
    nslines=int(np.array(indseuil).size)
    satlasext=atlasext[indseuil]
    slatlasext=latlasext[indseuil]
    sratlasext=ratlasext[indseuil]
    
# ici boucle sur toutes les raies selectionnees, correlation sur une petite zone autour de chaque raie entre l'intensite des deux voies individuellement et l'intensite du spectre de ref.
    
    corr1 = zeros (nslines, dtype =  float)
    corr2 = zeros (nslines, dtype =  float)
    corr12 = zeros (nslines, dtype =  float)
    corr1ms = zeros (nslines, dtype =  float)
    corr2ms = zeros (nslines, dtype =  float)
    corr12ms = zeros (nslines, dtype =  float)
    resol1 = zeros (nslines, dtype =  float)
    resol2 = zeros (nslines, dtype =  float)
    gaussposi1 = zeros (nslines, dtype =  float)
    gaussposi2 = zeros (nslines, dtype =  float)
    gaussposi3 = zeros (nslines, dtype =  float)
    pixelposi1 = zeros (nslines, dtype =  float)
    pixelposi2 = zeros (nslines, dtype =  float)
    pixelposi3 = zeros (nslines, dtype =  float)
    erreurposi1 = zeros (nslines, dtype =  float)
    erreurposi2 = zeros (nslines, dtype =  float)
    erreurposi3 = zeros (nslines, dtype =  float)
    diff1 = zeros (nslines, dtype =  float)
    diff2 = zeros (nslines, dtype =  float)
    diff1ms= zeros (nslines, dtype =  float)
    diff2ms = zeros (nslines, dtype =  float)
    


    amp1 = zeros (nslines, dtype =  float)
    amp2 = zeros (nslines, dtype =  float)
    amp3 = zeros (nslines, dtype =  float)
    
    for k in range (nslines):
# attention, si on teste pour seulement 2 raies, le programme se plante plus tard!
#    for k in range (5,6):

# on incremente le compteur general de raies spectrales de 1
# ici on pretend que voies 1 et 2 sont suffisamment proches pour que le nombre
# de raies est identique
        print("linecount: ",linecount)
        linecount += 1
        indextr = where((refwave <sratlasext[k]) & (refwave>slatlasext[k]))
        wr=refwave[indextr]
        ir=refintens[indextr]
       
        indext1 =  where((wext1<max(wr)) & (wext1>min(wr)))
        w1=wext1[indext1]
        i1=iext1[indext1]
        pix1=np.array(indext1)
        pixel1=pix1[0]
    
        indext2 =  where((wext2<max(wr)) & (wext2>min(wr)))
        w2=wext2[indext2]
        i2=iext2[indext2]
        pix2=np.array(indext2)
        pixel2=pix2[0]
        
        indext3 =  where((wext3<max(wr)) & (wext3>min(wr)))
        w3=wext3[indext3]
        i3=iext3[indext3]
        pix3=np.array(indext3)
        pixel3=pix3[0]


##############################resolution voie 1#############################################
            
        xxx=w1
        yyy=i1
        print("satlasext[k]", satlasext[k])
#        print("xxx:", xxx)
#        print("yyy:", yyy)
#        print("pixel1:",pixel1)
#        print("mean(xxx):",mean(xxx))
        indrefextractmin=np.argmin(yyy)
        indrefextractmax=np.argmax(yyy)
        indrefmin=np.argmin(xxx)
        indrefmax=np.argmax(xxx)
        rA=(yyy[indrefextractmax]-yyy[indrefextractmin])
        rsigma=0.2*(xxx[indrefmax]-xxx[indrefmin])
#        rmu=xx[indrefextractmin]
        rmu=0.
        ry_offset=yyy[indrefextractmin]
#        print("ry_offset: ",y_offset)
#            gy=gauss(reflambda[indrefextract],A,mu,sigma,y_offset)
        rinit_vals = [rA, rmu,rsigma,ry_offset]  # for [amp, cen, wid]
        rsigmap=0.2*(np.argmax(pixel1)-np.argmin(pixel1))
        rinit_valsp = [rA, rmu, rsigmap, ry_offset]
        print ("rinit_vals fit EWref:", rinit_vals)
        errorbar=np.sqrt(yyy)
        try:
            rgopt, rcovar = curve_fit(gauss,xxx-mean(xxx),yyy,p0=rinit_vals,sigma=errorbar)
            rgoptp, rcovarp = curve_fit(gauss,pixel1-mean(pixel1),yyy,p0=rinit_valsp,sigma=errorbar)
# la matrice de covariance rcovarp (en pixels) contient sur sa diagonale les erreurs associes
# aux differents parametres. L'erreur sur la position du fit, donc relatif a rgoptp[1]
# s'exprime alors
            errs=np.diag(rcovarp)
            erreurposi1[k] = errs[1]
#            print("rmu:", gopt[1]," rA:",gopt[0]," rSigma:",gopt[2]," rY_offset:",gopt[3])
        #Return co-effs for fit and covariance
#            plt.figure()
#            plt.plot(xxx,yyy)
#            plt.plot(xxx,gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]),"r")
#            plt.show()
            FWHM =  rgopt[2] * 2 * sqrt(2*np.log(2))
            resol = mean(xxx)/FWHM
            resol1[k]=resol
            amp1[k]=rgopt[0]
            gaussposi1[k] = rgopt[1]+mean(xxx)
            pixelposi1[k] = rgoptp[1]+mean(pixel1)
            print("gaussposi1[k]:", gaussposi1[k])
            print("pixelposi1[k]:", pixelposi1[k])
            name1 = "ThAr2D_voie_1.dat"
            name2 = "ThAr2D_voie_2.dat"

# ici on ecrit le tableau de sortie (i, lam_i, true_ordre_i, posi_i,err_posi_i)
#            outtableau1.write("%4d %20.10e %2d %20.10e %20.10e" % (linecount, gaussposi1[k], tordernum[j], pixelposi1[k], erreurposi1[k])+str(pixel1)+str(yyy)+"\n")
            res1.append({
                'ThAr_line_number': linecount,
                'ref_lambda': satlasext[k],
                'meas_lambda': gaussposi1[k],
                'true_order_number': tordernum[j],
                'meas_pixel_position': pixelposi1[k],
                'fit_error_pixel_position': erreurposi1[k],
                'pixels_extract': pixel1.tolist(),
                'flux_values_extract': yyy.tolist()
            })




# difference satlasext - gaussposi
            diff1[k] = satlasext[k]-gaussposi1[k]
            diff1ms[k] = diff1[k]/satlasext[k]*clum*1000.
            print("Resolution:", resol)
            print("gauss:",gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]))
            print("gopt",rgopt)
        except (RuntimeError,TypeError):
            print("Ca n'a pas marche")
#[gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
# EW= A*(sqrt(2pi) sigma) si le continu est a 1!! Donc renormalisation AVANT
# la division par gopt[3] est une normalisation locale


##############################resolution voie 2#############################################
            
        xxx=w2
        yyy=i2
        print("satlasext[k]", satlasext[k])
#        print("xxx:", xxx)
#        print("yyy:", yyy)
#        print("pixel2:",pixel2)
#        print("mean(xxx):",mean(xxx))
        indrefextractmin=np.argmin(yyy)
        indrefextractmax=np.argmax(yyy)
        indrefmin=np.argmin(xxx)
        indrefmax=np.argmax(xxx)
        rA=(yyy[indrefextractmax]-yyy[indrefextractmin])
        rsigma=0.2*(xxx[indrefmax]-xxx[indrefmin])
#        rmu=xx[indrefextractmin]
        rmu=0.
        ry_offset=yyy[indrefextractmin]
#        print("ry_offset: ",y_offset)
#            gy=gauss(reflambda[indrefextract],A,mu,sigma,y_offset)
        rinit_vals = [rA, rmu,rsigma,ry_offset]  # for [amp, cen, wid]
        rsigmap=0.2*(np.argmax(pixel2)-np.argmin(pixel2))
        rinit_valsp = [rA, rmu, rsigmap, ry_offset]
        print ("rinit_vals fit EWref:", rinit_vals)
        errorbar = np.sqrt(yyy)
        try:
            rgopt, rcovar = curve_fit(gauss,xxx-mean(xxx),yyy,p0=rinit_vals,sigma=errorbar)
            rgoptp, rcovarp = curve_fit(gauss,pixel2-mean(pixel2),yyy,p0=rinit_valsp,sigma=errorbar)
            errs=np.diag(rcovarp)
            erreurposi2[k] = errs[1]

            
#            print("rmu:", gopt[1]," rA:",gopt[0]," rSigma:",gopt[2]," rY_offset:",gopt[3])
        #Return co-effs for fit and covariance
#            plt.figure()
#            plt.plot(xxx,yyy)
#            plt.plot(xxx,gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]),"r")
#            plt.show()
            FWHM =  rgopt[2] * 2 * sqrt(2*np.log(2))
            resol = mean(xxx)/FWHM
            resol2[k]=resol
            amp2[k]=rgopt[0]
            gaussposi2[k] = rgopt[1]+mean(xxx)
            pixelposi2[k] = rgoptp[1]+mean(pixel2)
            print("k: ",k)
            print("gaussposi1[k]:", gaussposi2[k])
            print("pixelposi1[k]:", pixelposi2[k])
            # ici on ecrit le tableau de sortie (i, lam_i, true_ordre_i, posi_i,err_posi_i)
#        outtableau2.write("%15.7f %20.10e %20.10e %20.10e %20.10e\n" % (linecount, gaussposi2[k], tordernum[j], pixelposi2[k], erreurposi2[k])+str(pixel2)+str(yyy)+"\n")
            res2.append({
                'ThAr_line_number': linecount,
                'ref_lambda': satlasext[k],
                'meas_lambda': gaussposi2[k],
                'true_order_number': tordernum[j],
                'meas_pixel_position': pixelposi2[k],
                'fit_error_pixel_position': erreurposi2[k],
                'pixels_extract': pixel2.tolist(),
                'flux_values_extract': yyy.tolist()
            })
# difference satlasext - gaussposi
            diff2[k] = satlasext[k]-gaussposi2[k]
            diff2ms[k] = diff2[k]/satlasext[k]*clum*1000.
#            print("Resolution:", resol)
#            print("gauss:",gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]))
#            print("gopt",rgopt)
        except (RuntimeError,TypeError):
            print("Ca n'a pas marche")
#[gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
# EW= A*(sqrt(2pi) sigma) si le continu est a 1!! Donc renormalisation AVANT
# la division par gopt[3] est une normalisation locale

##############################resolution voie 3#############################################
            
        xxx=w3
        yyy=i3
        print("satlasext[k]", satlasext[k])
#        print("xxx:", xxx)
#        print("yyy:", yyy)
#        print("pixel2:",pixel2)
#        print("mean(xxx):",mean(xxx))
        indrefextractmin=np.argmin(yyy)
        indrefextractmax=np.argmax(yyy)
        indrefmin=np.argmin(xxx)
        indrefmax=np.argmax(xxx)
        rA=(yyy[indrefextractmax]-yyy[indrefextractmin])
        rsigma=0.2*(xxx[indrefmax]-xxx[indrefmin])
#        rmu=xx[indrefextractmin]
        rmu=0.
        ry_offset=yyy[indrefextractmin]
#        print("ry_offset: ",y_offset)
#            gy=gauss(reflambda[indrefextract],A,mu,sigma,y_offset)
        rinit_vals = [rA, rmu,rsigma,ry_offset]  # for [amp, cen, wid]
        rsigmap=0.2*(np.argmax(pixel2)-np.argmin(pixel2))
        rinit_valsp = [rA, rmu, rsigmap, ry_offset]
        print ("rinit_vals fit EWref:", rinit_vals)
        errorbar = np.sqrt(yyy)
        try:
            rgopt, rcovar = curve_fit(gauss,xxx-mean(xxx),yyy,p0=rinit_vals,sigma=errorbar)
            rgoptp, rcovarp = curve_fit(gauss,pixel3-mean(pixel3),yyy,p0=rinit_valsp,sigma=errorbar)
            errs=np.diag(rcovarp)
            erreurposi3[k] = errs[1]

            
#            print("rmu:", gopt[1]," rA:",gopt[0]," rSigma:",gopt[2]," rY_offset:",gopt[3])
        #Return co-effs for fit and covariance
#            plt.figure()
#            plt.plot(xxx,yyy)
#            plt.plot(xxx,gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]),"r")
#            plt.show()
            FWHM =  rgopt[2] * 2 * sqrt(2*np.log(2))
            resol = mean(xxx)/FWHM
            resol2[k]=resol
            amp2[k]=rgopt[0]
            gaussposi3[k] = rgopt[1]+mean(xxx)
            pixelposi3[k] = rgoptp[1]+mean(pixel3)
            print("k: ",k)
            print("gaussposi3[k]:", gaussposi2[k])
            print("pixelposi3[k]:", pixelposi2[k])
            # ici on ecrit le tableau de sortie (i, lam_i, true_ordre_i, posi_i,err_posi_i)
#        outtableau2.write("%15.7f %20.10e %20.10e %20.10e %20.10e\n" % (linecount, gaussposi2[k], tordernum[j], pixelposi2[k], erreurposi2[k])+str(pixel2)+str(yyy)+"\n")
            res3.append({
                'ThAr_line_number': linecount,
                'ref_lambda': satlasext[k],
                'meas_lambda': gaussposi3[k],
                'true_order_number': tordernum[j],
                'meas_pixel_position': pixelposi3[k],
                'fit_error_pixel_position': erreurposi3[k],
                'pixels_extract': pixel3.tolist(),
                'flux_values_extract': yyy.tolist()
            })
# difference satlasext - gaussposi
            diff2[k] = satlasext[k]-gaussposi2[k]
            diff2ms[k] = diff2[k]/satlasext[k]*clum*1000.
#            print("Resolution:", resol)
#            print("gauss:",gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]))
#            print("gopt",rgopt)
        except (RuntimeError,TypeError):
            print("Ca n'a pas marche")



#outtableau1.close()
#outtableau2.close()

with open(name1+".json", 'w') as outfile:
        json.dump(res1, outfile, indent=2)
with open(name2+".json", 'w') as outfile:
        json.dump(res2, outfile, indent=2)



#                'name1': name1 + '.json',
#                'description': "reference file voie 1",
#with open("ThAr2D_voie_1.dat.json", "r") as read_file:
#    data = json.load(read_file)
#data[7]['true_order_number']
#plt.plot(data[77]['pixels_extract'],data[77]['flux_values_extract'])
#for line in data:
#plt.plot(line['pixels_extract'],line['flux_values_extract'])
